"""Probe-level GSE54650 circadian analysis workflow."""

from __future__ import annotations
import gzip, re
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "vendor/clrcycle"))
import clrcycle
from lib.common import paper_ellipse

LABELS = {
    "Adr": "Adrenal",
    "Aor": "Aorta",
    "Bstm": "Brainstem",
    "BFat": "Brown fat",
    "Cer": "Cerebellum",
    "Hrt": "Heart",
    "Hyp": "Hypothalamus",
    "Kid": "Kidney",
    "Liv": "Liver",
    "Lun": "Lung",
    "Mus": "Muscle",
    "WFat": "White fat",
}


def read_matrix(path: Path):
    """Read the expression table and sample titles from a GEO matrix file."""
    titles = begin = end = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as h:
        for i, line in enumerate(h):
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                begin = i
            elif line.startswith("!series_matrix_table_end"):
                end = i
                break
    x = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        skiprows=begin + 1,
        nrows=end - begin - 2,
        dtype={"ID_REF": str},
    )
    x = x.rename(columns={x.columns[0]: "ID_REF"}).set_index("ID_REF")
    x.columns = titles
    return x.astype(float)


def read_annotation(path: Path):
    """Read probe-to-gene-symbol mappings from the GEO platform annotation."""
    begin = end = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as h:
        for i, line in enumerate(h):
            if line.startswith("!platform_table_begin"):
                begin = i
            elif line.startswith("!platform_table_end"):
                end = i
                break
    x = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        skiprows=begin + 1,
        nrows=end - begin - 2,
        dtype=str,
        low_memory=False,
    )
    x = x.rename(columns={"ID": "ID_REF", "Gene symbol": "gene_symbol"})
    x["gene_symbol"] = x.gene_symbol.fillna("")
    return x[["ID_REF", "gene_symbol"]].drop_duplicates("ID_REF")


def clr(x):
    """Apply a zero-safe centered log-ratio transform by sample."""
    if np.any(x <= 0):
        x = x + x[x > 0].min() * 1e-3
    z = np.log(x)
    return z - z.mean(1, keepdims=True)


def harmonic(logx, ct):
    """Compute first-harmonic fit, amplitude, and peak time per feature."""
    phase = 2 * np.pi * (ct % 24) / 24
    design = np.c_[np.cos(phase), np.sin(phase)]
    y = logx - logx.mean(0, keepdims=True)
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    fitted = design @ beta
    ss = np.square(fitted).sum(0)
    total = np.square(y).sum(0)
    return (
        np.divide(ss, total, out=np.zeros_like(ss), where=total > 0),
        np.sqrt(np.square(beta).sum(0)),
        np.mod(np.arctan2(beta[1], beta[0]), 2 * np.pi) * 24 / (2 * np.pi),
    )


def obj(cov, order):
    """Evaluate first-Fourier variance for one circular feature order."""
    n = len(order)
    t = 2 * np.pi * np.arange(n) / n
    c = np.sqrt(2 / n) * np.cos(t)
    s = np.sqrt(2 / n) * np.sin(t)
    q = cov[np.ix_(order, order)]
    return float(c @ q @ c + s @ q @ s)


def optimize(cov, seed):
    """Initialize spectrally and improve the order with random pair swaps."""
    rng = np.random.default_rng(seed)
    val, vec = np.linalg.eigh(cov)
    xy = vec[:, -2:] * np.sqrt(np.maximum(val[-2:], 0))
    order = np.argsort(np.arctan2(xy[:, 1], xy[:, 0]))
    best = obj(cov, order)
    order = order.copy()
    for _ in range(4000):
        i, j = rng.choice(len(order), 2, replace=False)
        candidate = order.copy()
        candidate[i], candidate[j] = candidate[j], candidate[i]
        score = obj(cov, candidate)
        if score > best:
            order, best = candidate, score
    return order, best


def align(u, v, ct):
    """Orient projected coordinates to circadian time for interpretation."""
    z = u + 1j * v
    truth = 2 * np.pi * (ct % 24) / 24
    best = None
    for reflected in (False, True):
        q = np.conj(z) if reflected else z
        delta = np.angle(np.mean(np.exp(1j * (np.angle(q) - truth))))
        a = q * np.exp(-1j * delta)
        residual = np.angle(np.exp(1j * (np.angle(a) - truth)))
        R = np.abs(np.mean(np.exp(1j * residual)))
        err = np.mean(abs(residual)) * 24 / (2 * np.pi)
        candidate = (R, -err, a, reflected, residual)
        if best is None or candidate > best:
            best = candidate
    return best


def prefixes(columns):
    """Return tissue prefixes represented by two-cycle sample series."""
    out = []
    for x in columns:
        m = re.match(r"([A-Za-z]+)_CT\d+$", x)
        if m and m.group(1) not in out:
            out.append(m.group(1))
    return out


def analyze(expr, annot, prefix, seed):
    """Run selection, clrcycle fitting, and post-fit alignment for one tissue."""
    cols = sorted(
        [x for x in expr if x.startswith(prefix + "_CT")],
        key=lambda x: int(re.search(r"CT(\d+)", x).group(1)),
    )
    tissue = expr[cols].T
    ct = np.array([int(re.search(r"CT(\d+)", x).group(1)) for x in tissue.index], float)
    values = tissue.to_numpy()
    tissue = tissue.loc[:, np.isfinite(values).all(0) & (values.min(0) > 0)]
    r2, amp, peak = harmonic(np.log(tissue.to_numpy()), ct)
    scores = pd.DataFrame(
        {
            "ID_REF": tissue.columns.astype(str),
            "harmonic_r2": r2,
            "harmonic_amplitude": amp,
            "peak_ct": peak,
            "mean_expression": tissue.mean(0).to_numpy(),
        }
    ).merge(annot, on="ID_REF", how="left")
    scores["gene_symbol"] = scores.gene_symbol.fillna("")
    scores = scores[scores.gene_symbol.str.len() > 0].sort_values(
        ["harmonic_r2", "harmonic_amplitude"], ascending=False
    )
    selected = scores.head(240)
    ids = selected.ID_REF.tolist()
    result = clrcycle.fit(
        tissue[ids], max_features=None, seed=seed, n_swaps=4000
    )
    coordinates = result.coordinates.set_index("sample").loc[tissue.index]
    u = coordinates.clrcycle_cosine.to_numpy()
    v = coordinates.clrcycle_sine.to_numpy()
    R, negerr, a, reflected, residual = align(u, v, ct)
    n = len(result.feature_order)
    samples = pd.DataFrame(
        {
            "tissue": prefix,
            "sample": tissue.index,
            "ct": ct,
            "ct_mod24": ct % 24,
            "acs_u": u,
            "acs_v": v,
            "acs_u_aligned": a.real,
            "acs_v_aligned": a.imag,
            "acs_angle_aligned": np.mod(np.angle(a), 2 * np.pi),
            "acs_radius": abs(a),
            "phase_residual_hours": residual * 24 / (2 * np.pi),
        }
    )
    features = (
        selected.set_index("ID_REF")
        .loc[result.feature_order.feature]
        .reset_index()
    )
    features["tissue"] = prefix
    features["acs_position"] = result.feature_order.clrcycle_position.to_numpy()
    features["acs_angle"] = result.feature_order.clrcycle_angle.to_numpy()
    features = features[
        [
            "tissue",
            "ID_REF",
            "gene_symbol",
            "acs_position",
            "acs_angle",
            "harmonic_r2",
            "harmonic_amplitude",
            "peak_ct",
            "mean_expression",
        ]
    ]
    hi = samples.loc[samples.acs_radius.idxmax()]
    lo = samples.loc[samples.acs_radius.idxmin()]
    metrics = {
        "tissue": prefix,
        "label": LABELS[prefix],
        "n_samples": len(tissue),
        "n_candidate_positive_probes": tissue.shape[1],
        "n_selected_probes": n,
        "selection_mode": "rhythmic",
        "acs_objective": result.objective,
        "rho_clrcycle": result.rho,
        "phase_agreement_r": R,
        "mean_abs_phase_error_hours": -negerr,
        "alignment_reflected": reflected,
        "peak_radius_ct": hi.ct,
        "peak_radius_ct_mod24": hi.ct_mod24,
        "peak_radius": hi.acs_radius,
        "trough_radius_ct": lo.ct,
        "trough_radius_ct_mod24": lo.ct_mod24,
        "trough_radius": lo.acs_radius,
    }
    ellipse = paper_ellipse(result, ct, period=24, units="CT hours modulo 24")
    metrics.update(ellipse)
    weights = result.feature_weights.copy()
    weights["tissue"] = prefix
    weights["ct"] = weights["sample"].map(samples.set_index("sample").ct)
    weights["ct_mod24"] = weights["sample"].map(samples.set_index("sample").ct_mod24)
    weights["gene_symbol"] = weights["feature"].map(
        features.set_index("ID_REF").gene_symbol
    )
    return samples, features, scores, metrics, weights
