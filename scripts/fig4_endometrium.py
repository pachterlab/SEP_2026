"""Reproduce Figure 4: endometrium discovery and frozen transfer.

The analysis learns a cycle from GSE234354 bulk RNA-seq, then transfers that
fixed geometry to GSE111976 C1 and 10x donor pseudobulks. Inputs are raw counts,
source annotations, the Ensembl 91 GTF, and the platform files under
``data/raw/``. Outputs are written under ``results/fig3/`` and
``figures/Figure4.{pdf,svg,png}``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lib.figure_style import font_rcparams
from matplotlib.colors import BoundaryNorm, ListedColormap, LinearSegmentedColormap
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib import endometrium
from lib.common import write_versions, draw_ellipse, ellipse_legend_handles, ellipse_stats_box


def render_figures() -> None:
    """Render the complete two-tier reference composition from validated artifacts."""
    clean = ROOT / "results/fig3"
    target = ROOT / "figures"
    target.mkdir(exist_ok=True)
    train = pd.read_csv(clean / "rnaseq_train_sample_coordinates.csv")
    held = pd.read_csv(clean / "rnaseq_heldout_sample_coordinates.csv")
    tenx = pd.read_csv(clean / "10x_sample_coordinates.csv")
    c1 = pd.read_csv(clean / "c1_sample_coordinates.csv")
    genes = pd.read_csv(clean / "gene_circle_posthoc_histology_stage.csv")
    metrics = pd.read_csv(clean / "summary_metrics.csv", index_col=0)
    summary = json.loads((clean / "summary.json").read_text())
    ellipses = pd.read_csv(clean / "training_ellipse_geometry.csv").set_index(["cohort", "metadata_reference"])
    cmap = ListedColormap(endometrium.STAGE_COLORS)
    norm = BoundaryNorm(np.arange(0.5, 8.5), 7)

    time_cmap = LinearSegmentedColormap.from_list(
        "cycle_day",
        ["#eeeeee", "#76a5c5", "#63419a", "#8c315a", "#c56b5e", "#d8a98d", "#eeeeee"],
    )
    BASE_SIZE = 18
    TITLE_SIZE = 21
    PANEL_SIZE = 24
    GENE_SIZE = 16
    METRIC_SIZE = 16
    plt.rcParams.update(
        {
            **font_rcparams(),
            "font.size": BASE_SIZE,
            "axes.titlesize": TITLE_SIZE,
            "axes.titleweight": "bold",
            "axes.labelsize": BASE_SIZE,
            "xtick.labelsize": BASE_SIZE,
            "ytick.labelsize": BASE_SIZE,
            "legend.fontsize": BASE_SIZE,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(22, 14), constrained_layout=True)
    fig.get_layout_engine().set(h_pad=0.16, w_pad=0.08)
    outer = fig.add_gridspec(2, 1, height_ratios=(1.05, 0.90), hspace=0.22)
    top = outer[0].subgridspec(1, 5, width_ratios=(0.10, 1, 1, 1, 0.18), wspace=0.06)
    hist_ax = fig.add_subplot(top[0, 1])
    gene_ax = fig.add_subplot(top[0, 2])
    time_ax = fig.add_subplot(top[0, 3])
    hist_sc = hist_ax.scatter(
        train.clrcycle_cosine,
        train.clrcycle_sine,
        c=train.histology_stage,
        cmap=cmap,
        norm=norm,
        s=38,
        alpha=0.9,
        edgecolor=".3",
        lw=0.3,
    )
    time_sc = time_ax.scatter(
        train.clrcycle_cosine,
        train.clrcycle_sine,
        c=train.model_cycle_time,
        cmap=time_cmap,
        vmin=0,
        vmax=100,
        s=38,
        alpha=0.9,
        edgecolor=".3",
        lw=0.3,
    )
    for letter, ax, title in (
        ("A", hist_ax, "Histological stage"),
        ("C", time_ax, "Authors' molecular-model time"),
    ):
        reference = "stage" if letter == "A" else "model_time"
        geometry = ellipses.loc[("rnaseq_train", reference)]
        names = ("major_vertex_1", "major_vertex_2", "minor_vertex_1", "minor_vertex_2")
        values = [geometry[n+"_stage"] if letter == "A" else geometry[n+"_metadata"] for n in names]
        labels = [f"{int(v)} {endometrium.STAGE_NAMES[int(v)-1]}" if letter == "A" else f"Time {v:.1f}" for v in values]
        colors = [cmap(norm(v)) if letter == "A" else time_cmap(v/100) for v in values]
        draw_ellipse(ax, geometry, train[["raw_cosine", "raw_sine"]].to_numpy(),
                     train[["clrcycle_cosine", "clrcycle_sine"]].to_numpy(),
                     labels=labels, colors=colors, fontsize=9)
        ellipse_stats_box(ax, geometry, fontsize=11)
        ax.axhline(0, color=".85", lw=0.5)
        ax.axvline(0, color=".85", lw=0.5)
        ax.set_aspect("equal", adjustable="box")
        ax.set(xlabel="First Fourier cosine coordinate ($u$)", ylabel="First Fourier sine coordinate ($v$)", title=title)
        ax.text(
            -0.10,
            1.34,
            letter,
            transform=ax.transAxes,
            fontsize=PANEL_SIZE,
            fontweight="bold",
        )
        ax.legend(
            handles=ellipse_legend_handles(),
            loc="upper right",
            fontsize=10,
            frameon=True,
            framealpha=0.68,
            facecolor="white",
            edgecolor=".8",
            handlelength=2.0,
        )
    # Each R compares sample phases with its panel-specific metadata reference.
    # Use the same post-hoc alignment as the molecular-time ellipse calibration.
    _training_metrics = (
        rf"$\boldsymbol{{\rho}}_{{\mathbf{{clrcycle}}}}$={summary['training_rho']:.3f}, "
        rf"$R_{{\mathrm{{histology}}}}$={summary['orientation']['training_stage_R']:.3f}"
    )
    hist_ax.set_title("Histological stage", y=1.20, pad=4)
    hist_ax.text(
        0.5,
        1.04,
        _training_metrics,
        transform=hist_ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=METRIC_SIZE,
        fontweight="bold",
    )
    molecular_agreement = endometrium.clrcycle.circular_agreement(
        train.raw_angle.to_numpy(), train.model_cycle_time.to_numpy(), period=100
    )
    _molecular_metrics = (
        rf"$\boldsymbol{{\rho}}_{{\mathbf{{clrcycle}}}}$={summary['training_rho']:.3f}, "
        rf"$R_{{\mathrm{{molecular}}}}$={molecular_agreement['R']:.3f}"
    )
    time_ax.set_title("Authors' molecular-model time", y=1.20, pad=4)
    time_ax.text(
        0.5,
        1.04,
        _molecular_metrics,
        transform=time_ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=METRIC_SIZE,
        fontweight="bold",
    )
    angle = genes.gene_angle.to_numpy()
    stage = genes.posthoc_peak_histology_stage.to_numpy()
    gene_ax.scatter(
        np.cos(angle), np.sin(angle), c=stage, cmap=cmap, norm=norm, s=24, alpha=0.85
    )
    gene_ax.add_patch(plt.Circle((0, 0), 1, fill=False, color=".6", lw=0.8))
    gene_ax.set_aspect("equal")
    gene_ax.axis("off")
    gene_ax.set_title("Learned gene circle", y=1.20, pad=4)
    highlighted = genes.highlighted_top20_training_variance.astype(bool).to_numpy()
    _gene_labels = list(zip(genes.loc[highlighted, "gene"], angle[highlighted]))
    # Place labels in collision-free columns and connect them to their circle
    # positions; this keeps all 20 names readable at the Figure 4 font scale.
    for _side in (-1, 1):
        _side_labels = [
            item for item in _gene_labels if np.sign(np.cos(item[1])) == _side
        ]
        _side_labels.sort(key=lambda item: np.sin(item[1]))
        _label_y = np.linspace(-1.32, 1.32, len(_side_labels))
        for (gene, a), _y in zip(_side_labels, _label_y):
            _x = 1.30 * _side
            gene_ax.plot(
                [1.02 * np.cos(a), (_x - 0.035 * _side)],
                [1.02 * np.sin(a), _y],
                color=".72",
                lw=0.55,
                zorder=1,
            )
            gene_ax.text(
                _x,
                _y,
                gene,
                fontsize=GENE_SIZE,
                fontweight="normal",
                ha="left" if _side > 0 else "right",
                va="center",
                clip_on=False,
            )
    gene_ax.set(xlim=(-1.78, 1.78), ylim=(-1.48, 1.58))
    gene_ax.text(
        -0.12,
        1.34,
        "B",
        transform=gene_ax.transAxes,
        fontsize=PANEL_SIZE,
        fontweight="bold",
    )
    time_bar = fig.colorbar(time_sc, ax=time_ax, orientation="vertical",
                            fraction=0.04, pad=0.04, shrink=0.65)
    time_bar.set_ticks([0, 25, 50, 75, 100])
    time_bar.ax.tick_params(labelsize=14, pad=2)
    time_bar.set_label("Molecular-model time (C)", fontsize=14, labelpad=8)
    bottom = outer[1].subgridspec(
        2,
        7,
        height_ratios=(0.06, 1),
        width_ratios=(0.25, 1, 0.25, 1, 0.25, 1, 0.25),
        wspace=0.05,
        hspace=0.0,
    )
    axes = [fig.add_subplot(bottom[1, i]) for i in (1, 3, 5)]
    datasets = [
        (held, "histology_stage", "GSE234354 held-out", "RNAseq held-out"),
        (tenx, "canonical_stage", "GSE111976 10x pseudobulk", "10x transfer"),
        (c1, "canonical_stage", "GSE111976 C1 pseudobulk", "C1 transfer"),
    ]
    bounds = pd.read_csv(clean / "stage_boundaries.csv")
    limit = 1.08 * max(
        np.abs(
            pd.concat(
                [
                    held[["clrcycle_cosine", "clrcycle_sine"]],
                    tenx[["clrcycle_cosine", "clrcycle_sine"]],
                    c1[["clrcycle_cosine", "clrcycle_sine"]],
                ]
            ).to_numpy()
        ).max(),
        1,
    )
    light = [
        tuple(0.82 + 0.18 * x for x in matplotlib.colors.to_rgb(color))
        for color in endometrium.STAGE_COLORS
    ]
    for letter, ax, (frame, color, title, metric) in zip("DEF", axes, datasets):
        angles = bounds.boundary_radians.to_numpy()
        stages = bounds.stage_before.to_numpy(int)
        for start, end, st in zip(np.roll(angles, 1), angles, stages):
            if end <= start:
                end += 2 * np.pi
            ax.add_patch(
                matplotlib.patches.Wedge(
                    (0, 0),
                    limit * np.sqrt(2),
                    np.degrees(start),
                    np.degrees(end),
                    facecolor=light[st - 1],
                    edgecolor="none",
                    zorder=-3,
                )
            )
        ax.scatter(
            frame.clrcycle_cosine,
            frame.clrcycle_sine,
            c=frame[color],
            cmap=cmap,
            norm=norm,
            s=78,
            alpha=0.92,
            edgecolor=".3",
            lw=0.5,
        )
        row = metrics.loc[metric]
        ax.set_title(title, y=1.20, pad=4)
        _benchmark_metrics = (
            f"n={len(frame)}, R={row.circular_phase_agreement_R:.2f}, "
            f"exact={100*row.stage_exact:.0f}%, "
            rf"within $\pm 1$={100*row.stage_within_one_circular:.0f}%"
        )
        ax.text(
            0.5,
            1.04,
            _benchmark_metrics,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=METRIC_SIZE,
            fontweight="bold",
        )
        ax.axhline(0, c=".86", lw=0.7)
        ax.axvline(0, c=".86", lw=0.7)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_xlabel("First Fourier cosine coordinate ($u$)")
        ax.set_ylabel("First Fourier sine coordinate ($v$)")
        ax.legend(
            handles=[
                Line2D(
                    [],
                    [],
                    marker="o",
                    markersize=7,
                    ls="",
                    color=".45",
                    markerfacecolor=".65",
                    label="Point color: observed stage",
                ),
                matplotlib.patches.Patch(
                    facecolor=".90",
                    edgecolor=".70",
                    label="Background: predicted stage",
                ),
            ],
            loc="upper right",
            fontsize=9,
            frameon=True,
            framealpha=0.68,
            facecolor="white",
            edgecolor=".8",
            borderpad=0.35,
            labelspacing=0.25,
            handlelength=1.5,
        )
        ax.text(
            -0.20,
            1.34,
            letter,
            transform=ax.transAxes,
            fontsize=PANEL_SIZE,
            fontweight="bold",
            ha="left",
            va="bottom",
        )
    full_stage_names = [
        "Menstrual",
        "Early proliferative",
        "Mid proliferative",
        "Late proliferative",
        "Early secretory",
        "Mid secretory",
        "Late secretory",
    ]
    handles = [
        Line2D(
            [0], [0], marker="o", markersize=16, ls="", color=color, label=f"{i} {name}"
        )
        for i, (color, name) in enumerate(
            zip(endometrium.STAGE_COLORS, full_stage_names), 1
        )
    ]
    fig.legend(
        handles=handles,
        loc="outside lower center",
        ncol=7,
        frameon=False,
        fontsize=16,
        title_fontsize=18,
        title="\nUnified histological-stage colors shared by A, B, and D–F",
    )
    fig.suptitle(
        "Menstrual cycle projection, learned feature geometry (GSE234354), "
        "and stage prediction (GSE111976)",
        fontsize=24,
        fontweight="bold",
    )
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(
            target / f"Figure4.{suffix}",
            dpi=300 if suffix == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def main() -> None:
    """Run the shared analysis, copy its audit tables, and render Figure 4."""
    endometrium.run_shared_analysis()
    render_figures()
    write_versions(ROOT / "results/fig3/software_versions.json")


if __name__ == "__main__":
    main()
