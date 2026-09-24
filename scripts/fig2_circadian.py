"""Reproduce Figure 2: circadian structure across GSE54650 tissues.

Inputs
------
``data/raw/GSE54650/GSE54650_series_matrix.txt.gz`` and
``data/raw/GSE54650/GPL6246.annot.gz``.

Outputs
-------
Numerical tables under ``results/fig1/`` and
``figures/Figure2.{pdf,svg,png}``. Use ``--plot-only`` to redraw the figure
from existing numerical results without repeating the analysis.
"""

import argparse
from pathlib import Path
import sys
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lib.figure_style import font_rcparams
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.circadian import analyze, prefixes, read_annotation, read_matrix
from lib.common import draw_ellipse, ellipse_legend_handles

RAW = ROOT / "data/raw/GSE54650"
OUT = ROOT / "results/fig1/all_tissues"


def circadian_cmap():
    """Return the cyclic color map used for circadian time from 0 to 24 h."""
    return LinearSegmentedColormap.from_list(
        "circadian_phase",
        [
            (0, "#D73027"),
            (0.25, "#F9C74F"),
            (0.5, "#2A9D8F"),
            (0.75, "#4361EE"),
            (1, "#D73027"),
        ],
    )


def sample_projection(ax, table, cmap, labels=False, geometry=None):
    """Draw one tissue's aligned sample projection on ``ax``."""
    sc = ax.scatter(
        table.acs_u_aligned,
        table.acs_v_aligned,
        c=table.ct_mod24,
        cmap=cmap,
        vmin=0,
        vmax=24,
        s=70 if labels else 18,
        edgecolor="black" if labels else "none",
        linewidth=0.45,
        zorder=3,
    )
    ordered = table.sort_values("ct")
    ax.plot(
        ordered.acs_u_aligned, ordered.acs_v_aligned, color=".68", lw=0.55, zorder=1
    )
    names = ("major_vertex_1", "major_vertex_2", "minor_vertex_1", "minor_vertex_2")
    phases = [geometry[name+"_metadata"] for name in names]
    draw_ellipse(ax, geometry, table[["acs_u", "acs_v"]].to_numpy(),
                 table[["acs_u_aligned", "acs_v_aligned"]].to_numpy(),
                 labels=[f"CT{value:.1f}" for value in phases],
                 colors=[cmap(value/24) for value in phases], fontsize=9 if labels else 6.5)
    if labels:
        for row in table.itertuples():
            ax.text(
                row.acs_u_aligned,
                row.acs_v_aligned,
                str(int(row.ct)),
                ha="center",
                va="center",
                fontsize=7.5,
                color="white",
                zorder=4,
            )
    ax.axhline(0, color=".88", lw=0.45)
    ax.axvline(0, color=".88", lw=0.45)
    ax.set_aspect("equal", adjustable="box")
    return sc


def gene_circle(ax, features, cmap):
    """Draw the learned liver gene order and label canonical clock genes."""
    scale = 3.0
    theta = 2 * np.pi * features.acs_position / len(features)
    ax.scatter(
        scale * np.cos(theta),
        scale * np.sin(theta),
        c=features.peak_ct,
        cmap=cmap,
        vmin=0,
        vmax=24,
        s=16,
        lw=0,
    )
    ax.add_patch(plt.Circle((0, 0), scale, fill=False, color=".72", lw=0.8))
    labels = {
        "Arntl",
        "Clock",
        "Cry1",
        "Cry2",
        "Per1",
        "Per2",
        "Per3",
        "Nr1d1",
        "Dbp",
        "Rorc",
        "Nampt",
    }
    offsets = {
        "Nr1d1": (-0.08, 0.01),
        "Dbp": (0.08, 0.04),
        "Per2": (-0.02, 0.07),
        "Nampt": (0.03, -0.06),
    }
    for row in features[features.gene_symbol.isin(labels)].itertuples():
        angle = 2 * np.pi * row.acs_position / len(features)
        x, y = 3.18 * np.cos(angle), 3.18 * np.sin(angle)
        dx, dy = offsets.get(row.gene_symbol, (0, 0))
        example = row.gene_symbol.lower() in {"dbp", "rorc", "clock"}
        ax.annotate(
            row.gene_symbol,
            xy=(scale * np.cos(angle), scale * np.sin(angle)),
            xytext=(x + scale * dx, y + scale * dy),
            textcoords="data",
            arrowprops=dict(arrowstyle="-", color=".65", lw=0.5,
                            shrinkA=1, shrinkB=2),
            annotation_clip=False,
            fontsize=14,
            fontweight="bold" if example else "normal",
            ha="left" if x >= 0 else "right",
            va="center",
        )
    ax.set(xlim=(-3.45, 3.45), ylim=(-3.35, 3.55))
    ax.set_aspect("equal")
    ax.axis("off")


def harmonic_panel(ax, samples, expression, gene, cmap, compact=False):
    """Plot one gene's expression and fitted 24-hour first harmonic."""
    ct = samples.ct.to_numpy(float)
    phase = 2 * np.pi * samples.ct_mod24.to_numpy(float) / 24
    elapsed = ct - ct.min()
    design = np.c_[np.ones_like(phase), np.cos(phase), np.sin(phase)]
    beta = np.linalg.lstsq(design, expression, rcond=None)[0]
    grid = np.linspace(0, 48, 800)
    grid_phase = 2 * np.pi * np.mod(grid + ct.min(), 24) / 24
    fitted = beta[0] + beta[1] * np.cos(grid_phase) + beta[2] * np.sin(grid_phase)
    peak_phase = np.mod(np.arctan2(beta[2], beta[1]), 2 * np.pi)
    peak_ct = peak_phase * 24 / (2 * np.pi)
    ax.scatter(
        elapsed,
        expression,
        c=samples.ct_mod24,
        cmap=cmap,
        vmin=0,
        vmax=24,
        s=27 if compact else 30,
        edgecolor="white",
        lw=0.4,
        zorder=3,
    )
    ax.plot(grid, fitted, color="black", lw=1.5, label="Harmonic fit")
    for peak in (
        np.arange(
            peak_ct + 24 * np.ceil((ct.min() - peak_ct) / 24), ct.min() + 48 + 0.01, 24
        )
        - ct.min()
    ):
        ax.axvline(peak, color="black", ls="--", lw=0.55)
    # The C/D labels sit below the fitted curve instead of crossing it.
    peak_y = 0.10 if gene in {"DBP", "RORC"} else 0.91
    ax.text(
        0.97,
        peak_y,
        f"peak CT {peak_ct:.1f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom" if peak_y < 0.5 else "top",
        fontsize=10.5 if compact else 12,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=0.8),
    )
    ax.legend(
        handles=[
            Line2D([], [], color="black", lw=1.1, label="Harmonic fit"),
            Line2D([], [], color="black", ls="--", lw=0.7, label="Fitted peak"),
        ],
        loc="upper left",
        fontsize=9 if compact else 11,
        frameon=True,
        framealpha=0.88,
        facecolor="white",
        edgecolor=".8",
        handlelength=2.1,
        borderpad=0.35,
        labelspacing=0.25,
    )
    ax.set(
        title=gene, xlim=(0, 48), xticks=[0, 12, 24, 36, 48], xlabel="Elapsed time (h)"
    )
    ax.set_title(
        gene,
        fontsize=16 if compact else TITLE_SIZE if "TITLE_SIZE" in globals() else 20,
        fontweight="bold",
        pad=2 if compact else 5,
    )
    ax.spines[["top", "right"]].set_visible(False)


def render_figure2(expression):
    """Assemble Figure 2 from validated result tables and expression data."""
    liver = ROOT / "results/fig1/liver"
    allout = ROOT / "results/fig1/all_tissues"
    samples = pd.read_csv(liver / "liver_clr_acs_sample_coordinates.csv")
    features = pd.read_csv(liver / "liver_clr_acs_feature_order.csv")
    metrics = pd.read_csv(liver / "liver_metrics.csv").iloc[0]
    pooled = pd.read_csv(allout / "all_tissue_sample_coordinates.csv")
    tissue_metrics = pd.read_csv(allout / "all_tissue_metrics.csv")
    cmap = circadian_cmap()
    BASE_SIZE = 15
    TITLE_SIZE = 20
    PANEL_SIZE = 21
    plt.rcParams.update(
        {
            **font_rcparams(),
            "font.size": BASE_SIZE,
            "axes.titlesize": TITLE_SIZE,
            "axes.titleweight": "bold",
            "axes.labelsize": BASE_SIZE,
            "xtick.labelsize": BASE_SIZE - 1,
            "ytick.labelsize": BASE_SIZE - 1,
            "legend.fontsize": BASE_SIZE - 1,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(20, 14.8), facecolor="white")
    # A/B remain dominant; C--E use the otherwise empty right margin as a shared-x stack.
    outer = fig.add_gridspec(
        5,
        1,
        height_ratios=(1.40, 0.16, 1.08, 0.09, 0.06),
        hspace=0.15,
        left=0.055,
        right=0.98,
        top=0.855,
        bottom=0.065,
    )
    fig.suptitle(
        "Liver circadian projection and learnt feature geometry,\n"
        "and analysis of the complete Hogenesch circadian atlas (GSE54650)",
        fontsize=25,
        fontweight="bold",
        x=0.5,
        ha="center",
        y=0.985,
        linespacing=1.08,
    )
    top = outer[0].subgridspec(1, 3, width_ratios=(1.18, 1.02, 0.78), wspace=0.20)
    ax = fig.add_subplot(top[0, 0])
    sc = sample_projection(ax, samples, cmap, True, metrics)
    ax.set(
        title=rf"Liver samples ($\boldsymbol{{\rho}}_{{\mathbf{{clrcycle}}}}$={metrics.rho_clrcycle:.3f}, R={metrics.phase_agreement_r:.3f})",
        xlabel="First Fourier cosine coordinate ($u$)",
        ylabel="First Fourier sine coordinate ($v$)",
    )
    ax.text(
        -0.10,
        1.10,
        "A",
        transform=ax.transAxes,
        fontsize=PANEL_SIZE,
        fontweight="bold",
        va="bottom",
    )
    ax.legend(
        handles=ellipse_legend_handles(),
        loc="upper right",
        frameon=True,
        framealpha=0.75,
        fontsize=12,
    )
    ax2 = fig.add_subplot(top[0, 1])
    gene_circle(ax2, features, cmap)
    ax2.set_title(
        "Learned feature circle (peak CT)",
        fontweight="bold",
        fontsize=TITLE_SIZE,
        y=1.04,
    )
    ax2.text(
        -0.06,
        1.10,
        "B",
        transform=ax2.transAxes,
        fontsize=PANEL_SIZE,
        fontweight="bold",
        va="bottom",
    )
    bpos = ax2.get_position()
    ax2.set_position([bpos.x0 - 0.018, bpos.y0, bpos.width, bpos.height])
    liver_cols = sorted(
        [c for c in expression if c.startswith("Liv_CT")],
        key=lambda x: int(x.split("CT")[1]),
    )
    liver_expr = expression[liver_cols].T
    middle = top[0, 2].subgridspec(3, 1, hspace=0.34)
    cde_axes = []
    for index, (gene, probe) in enumerate(
        [("DBP", "10553092"), ("RORC", "10494023"), ("CLOCK", "10530733")]
    ):
        axis = fig.add_subplot(
            middle[index, 0], sharex=cde_axes[0] if cde_axes else None
        )
        cde_axes.append(axis)
        harmonic_panel(
            axis,
            samples,
            np.log(liver_expr[probe].to_numpy(float)),
            gene,
            cmap,
            compact=True,
        )
        axis.text(
            -0.13,
            1.02,
            chr(ord("C") + index),
            transform=axis.transAxes,
            fontsize=PANEL_SIZE,
            fontweight="bold",
            va="bottom",
        )
        if index == 0:
            axis.set_ylabel("ln expression")
        if index < 2:
            axis.set_xlabel("")
            axis.tick_params(labelbottom=False)
    tissue_grid = outer[2].subgridspec(2, 6, wspace=0.25, hspace=0.88)
    grid_axes = []
    for i, row in tissue_metrics.iterrows():
        axis = fig.add_subplot(tissue_grid[i // 6, i % 6])
        grid_axes.append(axis)
        subset = pooled[pooled.tissue == row.tissue]
        sc2 = sample_projection(axis, subset, cmap, False, row)
        axis.set_xlabel("First Fourier cosine coordinate ($u$)", fontsize=11)
        axis.set_ylabel("First Fourier sine coordinate ($v$)", fontsize=11)
        axis.tick_params(labelsize=10)
        axis.set_title(
            f"{row.label}\n"
            + rf"$\boldsymbol{{\rho}}_{{\mathbf{{clrcycle}}}}$={row.rho_clrcycle:.2f}, R={row.phase_agreement_r:.2f}",
            fontsize=11.5,
            pad=3,
        )
        axis.text(
            -0.55 if row.tissue == "Hyp" else -0.15,
            1.18,
            chr(ord("F") + i),
            transform=axis.transAxes,
            fontsize=PANEL_SIZE,
            fontweight="bold",
            va="bottom",
        )
    # One unified phase scale applies to A, C--E, and every tissue projection.
    cax_bottom = fig.add_subplot(outer[4])
    cb = fig.colorbar(sc2, cax=cax_bottom, orientation="horizontal")
    cb.set_label("Circadian time modulo 24 h", labelpad=4)
    cb.set_ticks([0, 6, 12, 18, 24])

    FIG = ROOT / "figures"
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "Figure2.pdf", bbox_inches="tight")
    fig.savefig(FIG / "Figure2.svg", bbox_inches="tight")
    fig.savefig(FIG / "Figure2.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    """Run the circadian analysis or redraw Figure 2 from saved results."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Render from previously validated numeric tables without recomputing analysis",
    )
    args = parser.parse_args()
    if args.plot_only:
        render_figure2(read_matrix(RAW / "GSE54650_series_matrix.txt.gz"))
        print("Wrote Figure2 PDF/SVG/PNG from validated outputs")
        return
    OUT.mkdir(parents=True, exist_ok=True)
    expression = read_matrix(RAW / "GSE54650_series_matrix.txt.gz")
    annotation = read_annotation(RAW / "GPL6246.annot.gz")
    liver_out = ROOT / "results/fig1/liver"
    liver_out.mkdir(parents=True, exist_ok=True)
    liver_cells, liver_features, liver_scores, liver_metrics, liver_weights = analyze(
        expression, annotation, "Liv", 17
    )
    liver_cells.drop(columns="tissue").to_csv(
        liver_out / "liver_clr_acs_sample_coordinates.csv", index=False
    )
    liver_features.drop(columns="tissue").to_csv(
        liver_out / "liver_clr_acs_feature_order.csv", index=False
    )
    liver_scores.to_csv(liver_out / "liver_harmonic_feature_scores.csv", index=False)
    pd.DataFrame([liver_metrics]).to_csv(liver_out / "liver_metrics.csv", index=False)
    liver_geometry = pd.DataFrame([liver_metrics]).assign(fit_id="Liv_liver_panel")
    liver_geometry.to_csv(liver_out / "ellipse_geometry.csv", index=False)
    liver_weights.to_csv(liver_out / "liver_feature_weights.csv.gz", index=False)
    samples = []
    features = []
    feature_weights = []
    metrics = []
    for index, prefix in enumerate(prefixes(expression.columns)):
        cell, gene, scores, row, weights = analyze(
            expression, annotation, prefix, 17 + index
        )
        cell.to_csv(OUT / f"{prefix.lower()}_sample_coordinates.csv", index=False)
        gene.to_csv(OUT / f"{prefix.lower()}_feature_order.csv", index=False)
        scores.to_csv(
            OUT / f"{prefix.lower()}_harmonic_feature_scores.csv", index=False
        )
        samples.append(cell)
        features.append(gene)
        feature_weights.append(weights)
        metrics.append(row)
    pd.concat(samples, ignore_index=True).to_csv(
        OUT / "all_tissue_sample_coordinates.csv", index=False
    )
    pd.concat(features, ignore_index=True).to_csv(
        OUT / "all_tissue_feature_orders.csv", index=False
    )
    pd.concat(feature_weights, ignore_index=True).to_csv(
        OUT / "all_tissue_feature_weights.csv.gz", index=False
    )
    table = pd.DataFrame(metrics)
    table.to_csv(OUT / "all_tissue_metrics.csv", index=False)
    atlas_geometry = table.assign(fit_id=table.tissue)
    atlas_geometry.to_csv(OUT / "ellipse_geometry.csv", index=False)
    pd.concat([atlas_geometry, liver_geometry], ignore_index=True).to_csv(
        ROOT / "results/fig1/ellipse_geometry.csv", index=False
    )
    render_figure2(expression)
    print(table[["tissue", "rho_clrcycle", "phase_agreement_r", "major_semiaxis",
                 "minor_semiaxis", "axis_ratio", "anisotropy_pvalue"]].to_string(index=False))


if __name__ == "__main__":
    main()
