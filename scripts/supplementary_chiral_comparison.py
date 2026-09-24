"""Compare CHIRAL with current Figure 2 using identical probes and scoring.

Use the saved 240-probe panels and sample coordinates in results/fig1/all_tissues.
Both methods use lib.circadian.align: reflection/rotation maximizing the phase
resultant R, followed by mean absolute circular phase error in hours. Collection
time is used for scoring, and was also used upstream to select Figure 2 panels;
this is a comparison conditional on those panels, not label-free discovery.

Requires the CHIRAL R package. Use --rescore-only to reuse saved raw CHIRAL
phases after checking their input fingerprints.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.circadian import LABELS, align, read_matrix
from lib.common import write_versions, run_r

RAW = ROOT / "data/raw/GSE54650/GSE54650_series_matrix.txt.gz"
FIG1 = ROOT / "results/fig1/all_tissues"
OUT = ROOT / "results/chiral_comparison"
SEED = 17
TWO_PI = 2 * np.pi

# No collection times enter the R fitting process. Preserve unaligned phases
# at numerical precision so the common Python scorer handles both methods.
CHIRAL_CODE = r'''
suppressMessages(library(CHIRAL))
a <- commandArgs(TRUE)
mat <- as.matrix(read.csv(a[1], row.names = 1, check.names = FALSE))
set.seed(as.integer(a[3]))
out <- CHIRAL(mat, clockgenes = rownames(mat), pbar = FALSE)
if (length(out$phi) != ncol(mat) || any(!is.finite(out$phi))) {
    stop("CHIRAL did not return finite phases for every sample")
}
write.csv(data.frame(sample = colnames(mat), raw_phase_rad = out$phi,
                     iterations = out$iter), a[2], row.names = FALSE)
cat("CHIRAL", as.character(packageVersion("CHIRAL")),
    "R", as.character(getRversion()), "\n")
'''


def score_phase(phase, ct):
    """Use exactly the Figure 2 alignment and metrics for either method."""
    phase, ct = np.asarray(phase, float), np.asarray(ct, float)
    if phase.ndim != 1 or phase.shape != ct.shape or len(phase) < 3:
        raise ValueError("Phase and CT must be matched vectors of at least three samples")
    if not np.isfinite(phase).all() or not np.isfinite(ct).all():
        raise ValueError("Phase and CT must be finite")
    R, negerr, aligned, reflected, residual = align(np.cos(phase), np.sin(phase), ct)
    return {
        "R": float(R),
        "mean_abs_err_h": float(-negerr),
        "alignment_reflected": bool(reflected),
    }, np.mod(np.angle(aligned), TWO_PI), residual * 24 / TWO_PI


def matched_panel(expr, features, samples):
    """Extract exact Figure 2 probes and verify the saved CLR projection."""
    features = features.sort_values("acs_position")
    ids, names = features.ID_REF.tolist(), samples["sample"].tolist()
    if not ids or len(set(ids)) != len(ids) or len(set(names)) != len(names):
        raise ValueError("Expected nonempty, unique probe and sample IDs")
    if not np.array_equal(features.acs_position, np.arange(len(ids))):
        raise ValueError("Feature order must contain consecutive positions")
    if not np.allclose(features.acs_angle, TWO_PI * np.arange(len(ids)) / len(ids)):
        raise ValueError("Feature angles do not match the saved circular order")
    ct_from_names = np.array([float(name.rsplit("_CT", 1)[1]) for name in names])
    if not np.array_equal(ct_from_names, samples.ct):
        raise ValueError("Sample IDs and collection times do not match")
    values = expr.loc[ids, names].to_numpy(float)
    if not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError("Expected finite positive expression for every selected probe")
    log_values = np.log(values).T
    clr = log_values - log_values.mean(axis=1, keepdims=True)
    centered = clr - clr.mean(axis=0, keepdims=True)
    angles = TWO_PI * np.arange(len(ids)) / len(ids)
    reconstructed = centered @ np.array([np.cos(angles), np.sin(angles)]).T * np.sqrt(2 / len(ids))
    saved = samples[["acs_u", "acs_v"]].to_numpy(float)
    if not np.allclose(reconstructed, saved, rtol=1e-10, atol=1e-10):
        raise ValueError("Shared expression/probe panel does not reproduce saved clrcycle coordinates")
    matrix = pd.DataFrame(np.log2(values), index=ids, columns=names)
    matrix.index.name = "ID_REF"
    return matrix, float(np.max(np.abs(reconstructed - saved)))


def run_chiral(matrix, tissue_dir, rescore_only=False):
    """Fit CHIRAL or load phases whose saved input matches this exact panel."""
    tissue_dir.mkdir(parents=True, exist_ok=True)
    matrix_text = matrix.to_csv()
    fingerprint = hashlib.sha256(matrix_text.encode()).hexdigest()
    mat_path = tissue_dir / "log2_matrix.csv"
    phase_path = tissue_dir / "chiral_raw_phases.csv"
    meta_path = tissue_dir / "chiral_fit.json"
    if rescore_only:
        metadata = json.loads(meta_path.read_text())
        if metadata["input_sha256"] != fingerprint or metadata["seed"] != SEED:
            raise ValueError(f"Saved CHIRAL phases have stale inputs/settings: {tissue_dir}")
        if metadata["phase_sha256"] != hashlib.sha256(phase_path.read_bytes()).hexdigest():
            raise ValueError(f"Saved CHIRAL phase file has changed: {phase_path}")
    else:
        mat_path.write_text(matrix_text)
        done = run_r(CHIRAL_CODE, mat_path, phase_path, SEED, capture_output=True, text=True)
        metadata = {
            "input_sha256": fingerprint,
            "phase_sha256": hashlib.sha256(phase_path.read_bytes()).hexdigest(),
            "seed": SEED,
            "software": done.stdout.strip().splitlines()[-1],
        }
        meta_path.write_text(json.dumps(metadata, indent=2) + "\n")
    phases = pd.read_csv(phase_path)
    if phases["sample"].duplicated().any() or set(phases["sample"]) != set(matrix.columns):
        raise ValueError(f"CHIRAL sample IDs do not match the input: {phase_path}")
    phases = phases.set_index("sample").loc[matrix.columns]
    if phases.iterations.nunique() != 1 or not np.isfinite(phases.raw_phase_rad).all():
        raise ValueError(f"Invalid CHIRAL result: {phase_path}")
    return phases.raw_phase_rad.to_numpy(), int(phases.iterations.iloc[0]), metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rescore-only", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    expr = read_matrix(RAW)
    features = pd.read_csv(FIG1 / "all_tissue_feature_orders.csv", dtype={"ID_REF": str})
    samples = pd.read_csv(FIG1 / "all_tissue_sample_coordinates.csv")
    saved_metrics = pd.read_csv(FIG1 / "all_tissue_metrics.csv").set_index("tissue")
    rows, phase_tables, audits = [], [], {}
    for prefix, label in LABELS.items():
        ft = features.loc[features.tissue.eq(prefix)].copy()
        st = samples.loc[samples.tissue.eq(prefix)].copy()
        matrix, reconstruction_error = matched_panel(expr, ft, st)
        ct = st.ct.to_numpy(float)
        clr_phase = np.arctan2(st.acs_v, st.acs_u).to_numpy()
        clr, clr_aligned, clr_error = score_phase(clr_phase, ct)
        reference = saved_metrics.loc[prefix]
        if len(ft) != reference.n_selected_probes or len(st) != reference.n_samples:
            raise ValueError(f"Panel/sample counts differ from Figure 2 for {label}")
        if not np.allclose(
            [clr["R"], clr["mean_abs_err_h"]],
            [reference.phase_agreement_r, reference.mean_abs_phase_error_hours],
            rtol=1e-10, atol=1e-10,
        ):
            raise ValueError(f"Common scoring does not reproduce Figure 2 metrics for {label}")
        tissue_dir = OUT / "matrices" / label.replace(" ", "_")
        print(f"Fitting/scoring {label}: {len(ft)} matched probes, {len(st)} samples", flush=True)
        chiral_phase, iterations, metadata = run_chiral(matrix, tissue_dir, args.rescore_only)
        chiral, chiral_aligned, chiral_error = score_phase(chiral_phase, ct)
        pd.DataFrame({"sample": matrix.columns, "CT": ct}).to_csv(tissue_dir / "ct.csv", index=False)
        rows.append({
            "tissue": label, "n_probes": len(ft), "n_samples": len(st),
            "chiral_R": chiral["R"], "chiral_mean_abs_err_h": chiral["mean_abs_err_h"],
            "chiral_iterations": iterations,
            "clrcycle_R": clr["R"], "clrcycle_mean_abs_err_h": clr["mean_abs_err_h"],
            "delta_R_clrcycle_minus_chiral": clr["R"] - chiral["R"],
            "delta_error_h_clrcycle_minus_chiral": clr["mean_abs_err_h"] - chiral["mean_abs_err_h"],
        })
        phase_tables.append(pd.DataFrame({
            "tissue": label, "sample": matrix.columns, "ct": ct, "ct_mod24": ct % 24,
            "chiral_raw_phase_rad": chiral_phase, "clrcycle_raw_phase_rad": clr_phase,
            "chiral_aligned_phase_rad": chiral_aligned, "clrcycle_aligned_phase_rad": clr_aligned,
            "chiral_phase_residual_hours": chiral_error,
            "clrcycle_phase_residual_hours": clr_error,
        }))
        audits[label] = {
            **metadata, "max_coordinate_reconstruction_error": reconstruction_error,
            "chiral_alignment_reflected": chiral["alignment_reflected"],
            "clrcycle_alignment_reflected": clr["alignment_reflected"],
        }
        print(f"  CHIRAL R={chiral['R']:.4f}, error={chiral['mean_abs_err_h']:.3f} h; "
              f"clrcycle R={clr['R']:.4f}, error={clr['mean_abs_err_h']:.3f} h", flush=True)

    table = pd.DataFrame(rows)
    summary = []
    for method, other in (("chiral", "clrcycle"), ("clrcycle", "chiral")):
        summary.append({
            "method": "CHIRAL" if method == "chiral" else method,
            "mean_R": table[f"{method}_R"].mean(),
            "mean_abs_err_h": table[f"{method}_mean_abs_err_h"].mean(),
            "tissues_with_higher_R": int((table[f"{method}_R"] > table[f"{other}_R"] + 1e-12).sum()),
            "tissues_with_lower_error": int((table[f"{method}_mean_abs_err_h"] < table[f"{other}_mean_abs_err_h"] - 1e-12).sum()),
        })
    table.to_csv(OUT / "chiral_vs_clrcycle.csv", index=False)
    pd.concat(phase_tables, ignore_index=True).to_csv(OUT / "sample_phases.csv", index=False)
    pd.DataFrame(summary).to_csv(OUT / "summary.csv", index=False)
    audit = {
        "panel_source": "results/fig1/all_tissues/all_tissue_feature_orders.csv",
        "coordinate_source": "results/fig1/all_tissues/all_tissue_sample_coordinates.csv",
        "alignment": "lib.circadian.align; choose reflection maximizing phase resultant, then circular-mean rotation",
        "R": "abs(mean(exp(1j * (aligned_phase - reference_phase))))",
        "error": "mean(abs(angle(exp(1j * (aligned_phase - reference_phase))))) * 24 / (2*pi)",
        "aggregation": "Unweighted arithmetic mean across the 12 tissues",
        "scope": "Single CHIRAL seed 17 versus saved Figure 2 fits; no repeated-initialization robustness claim",
        "selection": "Exact Figure 2 probes, selected upstream using collection-time harmonic scores; no reselection or gene-symbol deduplication",
        "preprocessing": "Same positive expression entries: log2 for CHIRAL, native centered log-ratio transform for clrcycle",
        "tissues": audits,
    }
    (OUT / "comparison_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    write_versions(OUT / "software_versions.json")
    versions = json.loads((OUT / "software_versions.json").read_text())
    fit_software = {entry["software"] for entry in audits.values()}
    versions["chiral_fit_software"] = sorted(fit_software)
    if len(fit_software) == 1:
        tokens = next(iter(fit_software)).split()
        versions["CHIRAL"], versions["r"] = tokens[1], tokens[3]
    (OUT / "software_versions.json").write_text(json.dumps(versions, indent=2) + "\n")
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
