"""Shared parsing, circular-statistics, and publication-style helpers."""

from __future__ import annotations

import csv
import gzip
import json
import platform
import re
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SEED = 17
STAGE_NAMES = [
    "Menstrual", "Early prolif.", "Mid prolif.", "Late prolif.",
    "Early secr.", "Mid secr.", "Late secr.",
]
STAGE_COLORS = ["#7b5aa6", "#4f81bd", "#3aa39b", "#67a65b", "#e4c64f", "#f39a3c", "#d95b69"]


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9,
        "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "axes.linewidth": 0.7, "lines.linewidth": 1.2,
        "savefig.bbox": "tight", "svg.fonttype": "none", "pdf.fonttype": 42,
    })


def panel_label(axis, label: str, x: float = -0.12, y: float = 1.08) -> None:
    axis.text(x, y, label, transform=axis.transAxes, fontsize=11, fontweight="bold", va="top")


def save_figure(figure, root: Path, stem: str) -> None:
    out = root / "figures"
    out.mkdir(parents=True, exist_ok=True)
    figure.savefig(out / f"{stem}.pdf")
    figure.savefig(out / f"{stem}.svg")
    figure.savefig(out / f"{stem}.png", dpi=350)
    plt.close(figure)


def write_versions(path: Path, packages=("clrcycle", "numpy", "pandas", "matplotlib", "scipy", "scikit-learn", "umap-learn")) -> None:
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
    frames = []
    for path in paths:
        _, meta = parse_geo_matrix(path)
        if len(meta):
            frames.append(meta)
    return pd.concat(frames, axis=0)


def harmonic_score(values: np.ndarray, phase: np.ndarray) -> np.ndarray:
    design = np.c_[np.ones(len(phase)), np.cos(phase), np.sin(phase)]
    fitted = design @ np.linalg.lstsq(design, values, rcond=None)[0]
    residual = ((values - fitted) ** 2).sum(axis=0)
    total = ((values - values.mean(axis=0)) ** 2).sum(axis=0)
    score = np.zeros_like(total, dtype=float)
    valid = total > 0
    score[valid] = 1 - residual[valid] / total[valid]
    return score


def align_angles(predicted: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, float, float]:
    candidates = []
    for direction in (1.0, -1.0):
        oriented = direction * predicted
        offset = np.angle(np.mean(np.exp(1j * (target - oriented))))
        aligned = np.mod(oriented + offset, 2 * np.pi)
        error = np.angle(np.exp(1j * (aligned - target)))
        candidates.append((aligned, float(np.abs(np.mean(np.exp(1j * error)))), float(np.mean(np.abs(error)))))
    return max(candidates, key=lambda item: item[1])


def phase_from_stage(stage: np.ndarray) -> np.ndarray:
    return np.mod((np.asarray(stage, float) - 1) / 7 * 2 * np.pi, 2 * np.pi)


def day_to_stage(day: np.ndarray) -> np.ndarray:
    return np.digitize(np.asarray(day, float), [4, 7, 11, 15, 19, 23], right=True) + 1


def ensembl_symbol_map(gtf: Path) -> dict[str, str]:
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


def frozen_project(data: pd.DataFrame, ordered_features: list[str], training_mean: np.ndarray) -> pd.DataFrame:
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
    return pd.DataFrame({"sample": data.index.astype(str), "clrcycle_cosine": x, "clrcycle_sine": y,
                         "clrcycle_angle": np.mod(np.arctan2(y, x), 2 * np.pi),
                         "clrcycle_radius": np.hypot(x, y)})
