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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor/clrcycle"))
import clrcycle

OUT = ROOT / "results/supplementary_data"
FIGURES = ROOT / "figures"
FIG1_OUT = ROOT / "results/supp_fig1"
EXCEL_MAX_DATA_ROWS = 1_048_575
ENDOMETRIUM_STAGES = {
    1: "Menstrual",
    2: "Early proliferative",
    3: "Mid proliferative",
    4: "Late proliferative",
    5: "Early secretory",
    6: "Mid secretory",
    7: "Late secretory",
}


def add_sample_annotations(tables):
    """Add explicit, consistently named sample annotations to each table."""
    circadian = tables["Circadian"]
    circadian["annotation_type"] = "circadian collection time"
    circadian["annotated_stage_or_phase"] = circadian.ct_mod24.map(
        lambda value: f"CT{value:g}"
    )
    circadian["annotated_cycle_position"] = circadian.ct_mod24
    circadian["annotated_cycle_units"] = "CT hours modulo 24"

    fucci = tables["FUCCI"]
    fucci["annotation_type"] = "measured FUCCI phase"
    fucci["annotated_stage_or_phase"] = fucci.measured_fucci_theta.map(
        lambda value: f"{value:.3f} rad"
    )
    fucci["annotated_cycle_position"] = fucci.measured_fucci_theta
    fucci["annotated_cycle_units"] = "radians modulo 2pi"

    endometrium = tables["Endometrium"]
    stage_name = endometrium.histology_stage.map(ENDOMETRIUM_STAGES)
    endometrium["annotation_type"] = "histology stage and model cycle time"
    endometrium["annotated_stage_or_phase"] = (
        endometrium.histology_stage.astype(int).astype(str) + " " + stage_name
    )
    endometrium["annotated_cycle_position"] = endometrium.model_cycle_time
    endometrium["annotated_cycle_units"] = "model cycle time (0-100)"


def metadata_tables(tables):
    """Return compact sample and feature metadata tables."""
    circadian = tables["Circadian"]
    fucci = tables["FUCCI"]
    endometrium = tables["Endometrium"]
    annotation = [
        "annotation_type",
        "annotated_stage_or_phase",
        "annotated_cycle_position",
        "annotated_cycle_units",
    ]
    samples = pd.concat(
        [
            circadian[["sample", "tissue", "ct", "ct_mod24", *annotation]]
            .drop_duplicates()
            .assign(dataset="Circadian"),
            fucci[
                [
                    "sample",
                    "measured_fucci_theta",
                    "experiment",
                    "chip_id",
                    *annotation,
                ]
            ]
            .drop_duplicates()
            .assign(dataset="FUCCI"),
            endometrium[
                ["sample", "histology_stage", "model_cycle_time", *annotation]
            ]
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
    FIG1_OUT.mkdir(parents=True, exist_ok=True)
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
    candidates = order.reset_index(drop=True).copy()
    candidates["display_gene"] = [
        str(symbol).split("///")[0].strip() or str(probe)
        for symbol, probe in zip(candidates.gene_symbol.fillna(""), candidates.ID_REF)
    ]
    candidates["mean_weight"] = mean_weight.reindex(candidates.index)
    candidates["absolute_mean_weight"] = candidates.mean_weight.abs()
    representatives = candidates.sort_values(
        "absolute_mean_weight", ascending=False, kind="stable"
    ).drop_duplicates("display_gene")
    positive = representatives.nlargest(25, "mean_weight").index
    negative = representatives.nsmallest(25, "mean_weight").index
    positions = np.sort(np.r_[positive, negative])
    if len(np.unique(positions)) != 50:
        raise RuntimeError("Expected 50 distinct CT10 genes")
    if candidates.loc[positions, "display_gene"].duplicated().any():
        raise RuntimeError("Repeated gene labels in CT10 heatmap")
    subset = full[full["clrcycle_position"].isin(positions)].copy()
    subset["sample"] = subset["sample"].astype(str)
    subset.to_csv(FIG1_OUT / "liver_CT10_feature_weights.csv", index=False)

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
        vmin=-0.2,
        vmax=0.2,
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
    colorbar.set_ticks(np.linspace(-0.2, 0.2, 5))
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(
            FIGURES / f"supp_fig1_liver_CT10_feature_weights.{suffix}",
            dpi=300 if suffix == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)
    return ct10_ids, len(positions), float(radius_error.max())


def collate_ellipses():
    """Collect all plotted clrcycle cohorts and their vertex interpretations."""
    circadian = pd.read_csv(ROOT / "results/fig1/all_tissues/ellipse_geometry.csv", float_precision="round_trip")
    circadian = circadian.assign(dataset="Circadian", fit_id=circadian.tissue,
                                metadata_reference="CT", figure_panel="Figure1 tissue panels")
    liver = pd.read_csv(ROOT / "results/fig1/liver/ellipse_geometry.csv", float_precision="round_trip").assign(
        dataset="Circadian", fit_id="Liv_liver_panel", metadata_reference="CT", figure_panel="Figure1 A")
    fucci = pd.read_csv(ROOT / "results/fig2/ellipse_geometry.csv", float_precision="round_trip").assign(
        dataset="FUCCI", fit_id="GSE121265", metadata_reference="FUCCI phase", figure_panel="Figure2 B")
    endometrium = pd.read_csv(ROOT / "results/fig3/training_ellipse_geometry.csv", float_precision="round_trip")
    endometrium = endometrium.assign(dataset="Endometrium", fit_id=endometrium.cohort,
                                   figure_panel="Figure3")
    geometry = pd.concat([circadian, liver, fucci, endometrium], ignore_index=True)
    geometry["significant_anisotropy"] = geometry.anisotropy_pvalue.lt(.05)
    phase_columns = [c for c in geometry if c.startswith("major_vertex_")
                     and (c.endswith("_metadata") or c.endswith("_stage"))]
    geometry.loc[~geometry.significant_anisotropy, phase_columns] = np.nan
    ids = ["dataset", "fit_id", "metadata_reference", "figure_panel"]
    geometry = geometry[ids + [c for c in geometry if c not in ids]]
    vertices = []
    for _, row in geometry.iterrows():
        for name in ("major_vertex_1", "major_vertex_2", "minor_vertex_1", "minor_vertex_2"):
            vertex = {key: row[key] for key in ids}
            vertex.update(vertex=name, metadata_units=row.metadata_units)
            vertex.update({key.removeprefix(name+"_"): value for key, value in row.items() if key.startswith(name+"_")})
            vertices.append(vertex)
    vertices = pd.DataFrame(vertices)
    geometry.to_csv(OUT / "ellipse_geometry.csv", index=False)
    vertices.to_csv(OUT / "ellipse_vertices.csv", index=False)
    return geometry, vertices


def write_workbook(path, geometry, table):
    """Write ellipse features and feature weights, separating circadian tissues."""
    from openpyxl.comments import Comment

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        geometry.to_excel(writer, sheet_name="Ellipse_features", index=False)
        if "tissue" in table:
            labels = geometry.dropna(subset=["tissue", "label"]).set_index("tissue")["label"].to_dict()
            for tissue, weights in table.groupby("tissue", sort=False):
                weights.to_excel(writer, sheet_name=labels[tissue], index=False)
        else:
            table.to_excel(writer, sheet_name="Feature_weights", index=False)
        notes = {
            "axis_ratio": "b/a: minor semiaxis divided by major semiaxis; 1 denotes a circle.",
            "displacement_ratio": "d/a: distance of the ellipse center from the origin divided by the major semiaxis.",
            "feature_weight": "Signed contribution to the sample radius. Feature weights sum to the sample clrcycle radius.",
            "metadata_reference": "Reference used to label ellipse vertices. Repeated references for one cohort share the same intrinsic geometry. Endometrium ellipse results describe training biopsies only.",
            "orientation": "Major-axis orientation in radians modulo pi about the fitted center. Geometry retains raw clrcycle units.",
            "anisotropy_pvalue": "Radius-preserving random-angle test of anisotropy: 1000 randomizations, seed 17; no multiple-testing adjustment. Major-vertex phase/stage values are omitted unless p < 0.05. This does not test general cyclicity.",
        }
        for sheet in writer.book.worksheets:
            sheet.freeze_panes = "A2"
            sheet.auto_filter.ref = sheet.dimensions
            sheet.sheet_view.showGridLines = False
            for cell in sheet[1]:
                if cell.value in notes:
                    cell.comment = Comment(notes[cell.value], "clrcycle")
                sheet.column_dimensions[cell.column_letter].width = min(36, max(14, len(cell.value) + 2))


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
    add_sample_annotations(tables)
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
    ellipse_geometry, ellipse_vertices = collate_ellipses()
    workbook_names = {
        "Circadian": "Supplementary_Data_1_GSE54650.xlsx",
        "FUCCI": "Supplementary_Data_2_GSE121265.xlsx",
        "Endometrium": "Supplementary_Data_3_GSE234354_GSE111976.xlsx",
    }
    included = {key: len(table) <= EXCEL_MAX_DATA_ROWS for key, table in tables.items()}
    for key, table in tables.items():
        if not included[key]:
            raise RuntimeError(f"{key} feature weights exceed the Excel row limit")
        geometry = ellipse_geometry.loc[ellipse_geometry.dataset.eq(key)].dropna(axis=1, how="all").copy()
        geometry["displacement_ratio"] = geometry.ellipse_center_distance / geometry.major_semiaxis
        write_workbook(OUT / workbook_names[key], geometry, table)

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
        "workbooks": {key: str(OUT / name) for key, name in workbook_names.items()},
        "ellipse_summary_rows": len(ellipse_geometry),
        "ellipse_vertex_rows": len(ellipse_vertices),
        "liver_CT10_sample_ids": ct10_ids,
        "liver_CT10_features_plotted": plotted_features,
        "liver_CT10_selection": (
            "One probe per displayed gene, chosen by largest absolute mean weight; "
            "25 highest and 25 lowest mean weights across both replicates"
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
