#!/usr/bin/env python3
"""Collate native clrcycle feature weights and render Supplementary Figure 1."""

from __future__ import annotations

import json
import sys
from importlib.metadata import version
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor/clrcycle"))
import clrcycle

OUT = ROOT / "results/supplementary_data"
FIGURES = ROOT / "figures"
EXCEL_MAX_DATA_ROWS = 1_048_575


def metadata_tables(tables):
    """Return compact sample and feature metadata tables."""
    circadian = tables["Circadian"]
    fucci = tables["FUCCI"]
    endometrium = tables["Endometrium"]
    samples = pd.concat(
        [
            circadian[["sample", "tissue", "ct", "ct_mod24"]]
            .drop_duplicates()
            .assign(dataset="Circadian"),
            fucci[["sample", "measured_fucci_theta", "experiment", "chip_id"]]
            .drop_duplicates()
            .assign(dataset="FUCCI"),
            endometrium[["sample", "histology_stage", "model_cycle_time"]]
            .drop_duplicates()
            .assign(dataset="Endometrium"),
        ],
        ignore_index=True,
    )
    features = pd.concat(
        [
            circadian[
                ["feature", "tissue", "clrcycle_position", "feature_angle", "gene_symbol"]
            ]
            .drop_duplicates()
            .assign(dataset="Circadian"),
            fucci[
                [
                    "feature",
                    "clrcycle_position",
                    "feature_angle",
                    "gene",
                    "published_rank",
                    "published_pve",
                ]
            ]
            .drop_duplicates()
            .assign(dataset="FUCCI"),
            endometrium[
                [
                    "feature",
                    "clrcycle_position",
                    "feature_angle",
                    "training_log1p_variance",
                ]
            ]
            .drop_duplicates()
            .assign(dataset="Endometrium"),
        ],
        ignore_index=True,
    )
    samples = samples[["dataset", *[c for c in samples if c != "dataset"]]]
    features = features[["dataset", *[c for c in features if c != "dataset"]]]
    return samples, features


def liver_figure():
    """Render the complete ordered CT10 liver feature-weight decomposition."""
    liver_dir = ROOT / "results/fig1/liver"
    weights = pd.read_csv(
        liver_dir / "liver_feature_weights.csv.gz", dtype={"feature": str}
    )
    order = pd.read_csv(
        liver_dir / "liver_clr_acs_feature_order.csv", dtype={"ID_REF": str}
    )
    samples = pd.read_csv(liver_dir / "liver_clr_acs_sample_coordinates.csv")
    ct10_ids = samples.loc[np.isclose(samples.ct_mod24, 10), "sample"].tolist()
    if len(ct10_ids) != 2:
        raise RuntimeError(f"Expected two CT10 liver samples; found {ct10_ids}")
    full = weights[weights["sample"].isin(ct10_ids)].copy()
    full["sample"] = pd.Categorical(full["sample"], ct10_ids, ordered=True)
    full = full.sort_values(["sample", "clrcycle_position"])
    expected = order.ID_REF.tolist()
    for sample in ct10_ids:
        observed = full.loc[full["sample"].eq(sample), "feature"].tolist()
        if observed != expected:
            raise RuntimeError(f"CT10 rows do not follow clrcycle order for {sample}")
    radius_error = (
        full.groupby("sample", observed=False)["feature_weight"].sum()
        - full.groupby("sample", observed=False)["clrcycle_radius"].first()
    ).abs()
    if not np.allclose(radius_error, 0, rtol=0, atol=1e-12):
        raise RuntimeError("CT10 feature weights do not reconstruct sample radii")

    mean_weight = full.groupby("clrcycle_position")["feature_weight"].mean()
    positive = mean_weight.nlargest(25).index
    negative = mean_weight.nsmallest(25).index
    positions = np.sort(np.r_[positive, negative])
    if len(np.unique(positions)) != 50:
        raise RuntimeError("Positive and negative CT10 marker sets overlap")
    subset = full[full["clrcycle_position"].isin(positions)].copy()
    subset["sample"] = subset["sample"].astype(str)
    subset.to_csv(OUT / "liver_CT10_feature_weights.csv", index=False)

    matrix = np.column_stack(
        [
            subset.loc[subset["sample"].eq(sample), "feature_weight"].to_numpy()
            for sample in ct10_ids
        ]
    )
    selected_order = order.iloc[positions]
    symbols = selected_order.gene_symbol.fillna("").astype(str)
    labels = [
        f"{symbol.split('///')[0].strip() or probe} [{probe}]"
        for symbol, probe in zip(symbols, selected_order.ID_REF)
    ]
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none"}
    )
    fig, ax = plt.subplots(figsize=(5.2, 9), constrained_layout=True)
    image = ax.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=matrix.min(), vcenter=0, vmax=matrix.max()),
    )
    ax.set_xticks(np.arange(2), labels=ct10_ids, fontsize=9)
    ax.set_yticks(np.arange(len(labels)), labels=labels, fontsize=6.5)
    ax.tick_params(axis="y", length=1.5, pad=1)
    ax.set(
        xlabel="CT10 biological replicate",
        ylabel="Selected probe in clrcycle order",
        title="Liver CT10 sample feature weights",
    )
    ax.title.set_fontweight("bold")
    colorbar = fig.colorbar(image, ax=ax, pad=0.025, shrink=0.35)
    colorbar.set_label("Signed feature weight")
    colorbar.set_ticks([matrix.min(), 0, matrix.max() / 2, matrix.max()])
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(
            FIGURES / f"supp_fig1_liver_CT10_feature_weights.{suffix}",
            dpi=300 if suffix == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)
    return ct10_ids, len(positions), float(radius_error.max())


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    sources = {
        "Circadian": ROOT / "results/fig1/all_tissues/all_tissue_feature_weights.csv.gz",
        "FUCCI": ROOT / "results/fig2/clrcycle_feature_weights.csv.gz",
        "Endometrium": ROOT / "results/fig3/endometrium_feature_weights.csv.gz",
    }
    names = {key: f"{key.lower()}_feature_weights.csv.gz" for key in sources}
    tables = {
        key: pd.read_csv(path, dtype={"feature": str})
        for key, path in sources.items()
    }
    reconstruction_errors = {}
    for key, table in tables.items():
        group = ["sample"] if key != "Circadian" else ["tissue", "sample"]
        error = (
            table.groupby(group, sort=False)["feature_weight"].sum()
            - table.groupby(group, sort=False)["clrcycle_radius"].first()
        ).abs()
        reconstruction_errors[key] = float(error.max())
        if not np.allclose(error, 0, rtol=0, atol=1e-12):
            raise RuntimeError(f"{key} feature weights do not reconstruct sample radii")
        table.to_csv(OUT / names[key], index=False)
    sample_metadata, feature_metadata = metadata_tables(tables)
    workbook = OUT / "Supplementary_Data_clrcycle_feature_weights.xlsx"
    sheet_names = {
        "Circadian": "Circadian_feature_weights",
        "FUCCI": "FUCCI_feature_weights",
        "Endometrium": "Endometrium_feature_weights",
    }
    included = {
        key: len(table) <= EXCEL_MAX_DATA_ROWS for key, table in tables.items()
    }
    readme = pd.DataFrame(
        {
            "item": [
                "Definition",
                "Interpretation",
                "Reconstruction",
                *[f"{key} full table" for key in tables],
            ],
            "description": [
                "w_ig = Z_ig (c_g cos(phi_i) + s_g sin(phi_i))",
                "Per-sample feature contribution to the radial first-Fourier projection.",
                "Feature weights sum to the sample clrcycle radius.",
                *[
                    f"Workbook sheet {sheet_names[key]}; complete machine-readable file {names[key]}."
                    if included[key]
                    else f"Excel row limit exceeded; compact index in workbook and complete table in {names[key]}."
                    for key in tables
                ],
            ],
        }
    )
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        for key, table in tables.items():
            if included[key]:
                table.to_excel(writer, sheet_name=sheet_names[key], index=False)
            else:
                pd.DataFrame(
                    {
                        "full_table": [names[key]],
                        "rows": [len(table)],
                        "columns": [len(table.columns)],
                    }
                ).to_excel(writer, sheet_name=sheet_names[key], index=False)
        sample_metadata.to_excel(writer, sheet_name="Sample_metadata", index=False)
        feature_metadata.to_excel(writer, sheet_name="Feature_metadata", index=False)
        for sheet in writer.book.worksheets:
            sheet.freeze_panes = "A2"
            sheet.auto_filter.ref = sheet.dimensions
            sheet.sheet_view.showGridLines = False

    ct10_ids, plotted_features, ct10_error = liver_figure()
    audit = {
        "clrcycle_version": version("clrcycle"),
        "clrcycle_path": str(Path(clrcycle.__file__).resolve()),
        "analyses": {
            key: {
                "samples": int(table["sample"].nunique()),
                "selected_features_per_fit": int(
                    table.groupby("tissue")["feature"].nunique().iloc[0]
                    if "tissue" in table
                    else table["feature"].nunique()
                ),
                "fits": int(table["tissue"].nunique()) if "tissue" in table else 1,
                "feature_weight_rows": len(table),
                "max_radius_reconstruction_error": reconstruction_errors[key],
                "included_in_excel": included[key],
                "csv": str(OUT / names[key]),
            }
            for key, table in tables.items()
        },
        "workbook": str(workbook),
        "liver_CT10_sample_ids": ct10_ids,
        "liver_CT10_features_plotted": plotted_features,
        "liver_CT10_selection": (
            "25 highest and 25 lowest mean feature weights across both replicates"
        ),
        "rows_follow_clrcycle_order": True,
        "liver_CT10_max_radius_error": ct10_error,
        "supplementary_figure_renumbering": (
            "CT10 heatmap is Supplementary Figure 1; "
            "benchmarks are Supplementary Figure 2"
        ),
    }
    (OUT / "feature_weight_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
