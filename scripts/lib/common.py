"""Shared parsing, circular-statistics, and publication-style helpers."""

from __future__ import annotations

import csv
import gzip
import json
import platform
import re
import sys
import subprocess
from tempfile import NamedTemporaryFile
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from .figure_style import font_rcparams
import numpy as np
import pandas as pd


def run_r(code, *args, executable="Rscript", **kwargs):
    """Run an R source file, avoiding launcher-dependent escaping of -e code."""
    with NamedTemporaryFile(mode="w", suffix=".R", encoding="utf-8") as script:
        script.write(code)
        script.flush()
        return subprocess.run(
            [executable, script.name, *(str(arg) for arg in args)], check=True, **kwargs
        )


SEED = 17
STAGE_NAMES = [
    "Menstrual",
    "Early prolif.",
    "Mid prolif.",
    "Late prolif.",
    "Early secr.",
    "Mid secr.",
    "Late secr.",
]
STAGE_COLORS = [
    "#7b5aa6",
    "#4f81bd",
    "#3aa39b",
    "#67a65b",
    "#e4c64f",
    "#f39a3c",
    "#d95b69",
]


def set_style() -> None:
    """Apply the compact Matplotlib style shared by manuscript figures."""
    mpl.rcParams.update(
        {
            **font_rcparams(),
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "legend.fontsize": 7,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "axes.linewidth": 0.7,
            "lines.linewidth": 1.2,
            "savefig.bbox": "tight",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def panel_label(axis, label: str, x: float = -0.12, y: float = 1.08) -> None:
    """Place a bold manuscript panel label on an axes."""
    axis.text(
        x, y, label, transform=axis.transAxes, fontsize=11, fontweight="bold", va="top"
    )


def save_figure(figure, root: Path, stem: str) -> None:
    """Save one figure as PDF, SVG, and high-resolution PNG."""
    out = root / "figures"
    out.mkdir(parents=True, exist_ok=True)
    figure.savefig(out / f"{stem}.pdf")
    figure.savefig(out / f"{stem}.svg")
    figure.savefig(out / f"{stem}.png", dpi=350)
    plt.close(figure)


def write_versions(
    path: Path,
    packages=(
        "clrcycle",
        "numpy",
        "pandas",
        "matplotlib",
        "scipy",
        "scikit-learn",
        "umap-learn",
    ),
) -> None:
    """Record Python, platform, and package versions as JSON."""
    rows = {"python": sys.version.split()[0], "platform": platform.platform()}
    for package in packages:
        try:
            rows[package] = version(package)
        except PackageNotFoundError:
            rows[package] = "not installed"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rows, indent=2) + "\n")


def parse_geo_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return samples-by-features values and aligned GEO sample metadata."""
    metadata: dict[str, list[str]] = {}
    rows: list[list[str]] = []
    in_table = False
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                rows.append(next(csv.reader([line], delimiter="\t")))
            elif line.startswith("!Sample_"):
                fields = next(csv.reader([line], delimiter="\t"))
                metadata.setdefault(fields[0], []).append(fields[1:])
    if rows:
        table = pd.DataFrame(rows[1:], columns=rows[0])
        table = table.rename(columns={table.columns[0]: "feature"}).set_index("feature")
        table = table.apply(pd.to_numeric)
    else:
        table = pd.DataFrame()
    sample_ids = table.columns.astype(str).tolist()
    if not sample_ids and "!Sample_geo_accession" in metadata:
        sample_ids = metadata["!Sample_geo_accession"][0]
    sample = pd.DataFrame(index=sample_ids)
    for key, blocks in metadata.items():
        if len(blocks) == 1 and len(blocks[0]) == len(sample):
            sample[key.removeprefix("!Sample_")] = blocks[0]
        elif key == "!Sample_characteristics_ch1":
            for block in blocks:
                if len(block) != len(sample):
                    continue
                name = block[0].split(":", 1)[0].strip().lower().replace(" ", "_")
                sample[name] = [value.split(":", 1)[-1].strip() for value in block]
    return table.T, sample


def parse_geo_metadata(paths: list[Path]) -> pd.DataFrame:
    """Combine sample metadata parsed from multiple GEO matrix files."""
    frames = []
    for path in paths:
        _, meta = parse_geo_matrix(path)
        if len(meta):
            frames.append(meta)
    return pd.concat(frames, axis=0)


def harmonic_score(values: np.ndarray, phase: np.ndarray) -> np.ndarray:
    """Score each feature by variance explained by its first harmonic."""
    design = np.c_[np.ones(len(phase)), np.cos(phase), np.sin(phase)]
    fitted = design @ np.linalg.lstsq(design, values, rcond=None)[0]
    residual = ((values - fitted) ** 2).sum(axis=0)
    total = ((values - values.mean(axis=0)) ** 2).sum(axis=0)
    score = np.zeros_like(total, dtype=float)
    valid = total > 0
    score[valid] = 1 - residual[valid] / total[valid]
    return score


def align_angles(
    predicted: np.ndarray, target: np.ndarray
) -> tuple[np.ndarray, float, float]:
    """Align circular angles by choosing a reflection and rotation."""
    candidates = []
    for direction in (1.0, -1.0):
        oriented = direction * predicted
        offset = np.angle(np.mean(np.exp(1j * (target - oriented))))
        aligned = np.mod(oriented + offset, 2 * np.pi)
        error = np.angle(np.exp(1j * (aligned - target)))
        candidates.append(
            (
                aligned,
                float(np.abs(np.mean(np.exp(1j * error)))),
                float(np.mean(np.abs(error))),
            )
        )
    return max(candidates, key=lambda item: item[1])


def phase_from_stage(stage: np.ndarray) -> np.ndarray:
    """Map integer endometrium stages to canonical circular phases."""
    return np.mod((np.asarray(stage, float) - 1) / 7 * 2 * np.pi, 2 * np.pi)


def paper_ellipse(result, metadata, *, period, units, origin=0.0, alignment=None):
    """Summarize raw coordinates; optionally retain an existing phase calibration."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vendor/clrcycle"))
    import clrcycle

    metadata = np.asarray(metadata, float)
    summary = clrcycle.ellipse_summary(result, n_permutations=1000, seed=SEED)
    angles = result.coordinates.clrcycle_angle.to_numpy(float)
    if alignment is None:
        agreement = clrcycle.circular_agreement(angles, metadata-origin, period=period)
        rotation, reflected = agreement["rotation"], agreement["reflection"]
    else:
        rotation, reflected = alignment
    aligned = ((-angles if reflected else angles) + rotation) % (2*np.pi)
    residual = aligned - (metadata-origin)*2*np.pi/period
    summary.update(n_samples=len(angles), metadata_units=units, metadata_period=period,
                   metadata_origin=origin, n_permutations=1000, seed=SEED,
                   circular_agreement_R=float(abs(np.mean(np.exp(1j*residual)))),
                   alignment_rotation=float(rotation), alignment_reflection=bool(reflected))
    measured = np.sort(np.unique((metadata-origin) % period + origin))
    for name in ("major_vertex_1", "major_vertex_2", "minor_vertex_1", "minor_vertex_2"):
        angle = summary[name+"_angle"]
        phase = ((-angle if reflected else angle) + rotation) % (2*np.pi)
        value = phase*period/(2*np.pi)+origin
        summary[name+"_metadata"] = value
        summary[name+"_nearest_measured_metadata"] = measured[np.argmin(abs((measured-value+period/2) % period-period/2))]
    return summary


def ellipse_legend_handles():
    """Consistent keys for the ellipse curve, axes, and vertex markers."""
    from matplotlib.lines import Line2D

    return [
        Line2D([], [], color=".3", ls=":", lw=1.4, label="Fitted ellipse"),
        Line2D([], [], color=".5", lw=.6, label="Major / minor axes"),
        Line2D([], [], color="black", marker="D", ls="", markerfacecolor="white",
               markersize=5, markeredgewidth=.7, label="Vertices (axis endpoints)"),
    ]


def ellipse_stats_box(axis, summary, *, fontsize=9):
    """Show ellipse shape and normalized center displacement inside a panel."""
    displacement = summary["ellipse_center_distance"] / summary["major_semiaxis"]
    axis.text(
        .03, .03,
        rf"Anisotropy $b/a$ = {summary['axis_ratio']:.2f}" + "\n"
        + rf"Displacement $d/a$ = {displacement:.2f}",
        transform=axis.transAxes, ha="left", va="bottom", fontsize=fontsize,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor=".8", alpha=.9),
        zorder=10,
    )


def draw_ellipse(axis, summary, raw_xy, displayed_xy, *, labels=None, colors=None, fontsize=8):
    """Overlay the raw ellipse using only the panel's existing rotation/reflection."""
    transform = np.linalg.lstsq(raw_xy, displayed_xy, rcond=None)[0]
    np.testing.assert_allclose(transform.T @ transform, np.eye(2), atol=1e-10)
    np.testing.assert_allclose(raw_xy @ transform, displayed_xy, atol=1e-10)
    center = np.array([summary["ellipse_center_x"], summary["ellipse_center_y"]])
    names = ("major_vertex_1", "major_vertex_2", "minor_vertex_1", "minor_vertex_2")
    vertices = np.array([[summary[n+"_x"], summary[n+"_y"]] for n in names])
    theta = np.linspace(0, 2*np.pi, 361)
    curve = center + np.outer(np.cos(theta), vertices[0]-center) + np.outer(np.sin(theta), vertices[2]-center)
    axis.plot(*(curve @ transform).T, color=".3", ls=":", lw=1.4, label="Fitted ellipse", zorder=4)
    vertices, center = vertices @ transform, center @ transform
    for pair in (slice(0, 2), slice(2, 4)):
        axis.plot(*vertices[pair].T, color=".5", lw=.6, zorder=3)
    significant = np.isfinite(summary.get("anisotropy_pvalue", np.nan)) and summary["anisotropy_pvalue"] < .05
    shown = np.array([significant, significant, True, True])
    axis.scatter(*vertices[shown].T, c=np.asarray(colors)[shown] if colors is not None else "white", marker="D",
                 s=27, edgecolor="black", linewidth=.7, zorder=6)
    for show, vertex, label in zip(shown, vertices, labels or [f"{summary[n+'_metadata']:.1f}" for n in names]):
        if not show:
            continue
        right, above = vertex[0] >= center[0], vertex[1] >= center[1]
        axis.annotate(label, vertex, xytext=(5 if right else -5, 7 if above else -7),
                      textcoords="offset points", ha="left" if right else "right",
                      va="bottom" if above else "top", fontsize=fontsize, zorder=7,
                      bbox=dict(facecolor="white", edgecolor="none", alpha=.75, pad=.6))
    axis.margins(.22)


def day_to_stage(day: np.ndarray) -> np.ndarray:
    """Convert cycle-day estimates to the seven histological stages."""
    return np.digitize(np.asarray(day, float), [4, 7, 11, 15, 19, 23], right=True) + 1


def ensembl_symbol_map(gtf: Path) -> dict[str, str]:
    """Read an Ensembl-to-gene-symbol mapping from a gzipped GTF."""
    mapping: dict[str, str] = {}
    pattern = re.compile(r'gene_id "([^"]+)";.*gene_name "([^"]+)";')
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if "\tgene\t" not in line:
                continue
            match = pattern.search(line)
            if match:
                mapping[match.group(1).split(".")[0]] = match.group(2)
    return mapping


def frozen_project(
    data: pd.DataFrame, ordered_features: list[str], training_mean: np.ndarray
) -> pd.DataFrame:
    """Project new samples through a fixed clrcycle feature order."""
    selected = data.reindex(columns=ordered_features, fill_value=0.0).astype(float)
    values = selected.to_numpy()
    positive = values[values > 0]
    pseudo = 1e-3 * positive.min() if np.any(values <= 0) and positive.size else 1e-9
    log_values = np.log(values + pseudo)
    clr = log_values - log_values.mean(axis=1, keepdims=True)
    centered = clr - training_mean
    d = len(ordered_features)
    angles = 2 * np.pi * np.arange(d) / d
    cosine = np.sqrt(2 / d) * np.cos(angles)
    sine = np.sqrt(2 / d) * np.sin(angles)
    x, y = centered @ cosine, centered @ sine
    return pd.DataFrame(
        {
            "sample": data.index.astype(str),
            "clrcycle_cosine": x,
            "clrcycle_sine": y,
            "clrcycle_angle": np.mod(np.arctan2(y, x), 2 * np.pi),
            "clrcycle_radius": np.hypot(x, y),
        }
    )
