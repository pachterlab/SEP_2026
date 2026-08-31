"""Methods Figure: standalone clrcycle graphical workflow.

This script contains the complete figure implementation.  It has no runtime
dependency on the original clrcycle figure directory and draws every panel,
equation, and transition arrow natively with Matplotlib.
"""

from __future__ import annotations

from itertools import permutations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "figures"
METHODS_PANEL_DIR = FIGURE_DIR / "methods figs"

BLUE = "#1f77b4"
ORANGE = "#ff7f0e"
GREEN = "#2ca02c"
GRAY = "0.42"
LIGHT_GRAY = "0.88"
GENES = np.array(["Gene A", "Gene B", "Gene C"])


def clr(counts: np.ndarray) -> np.ndarray:
    """Centered log-ratio transform of a positive cell-by-gene matrix."""
    logged = np.log(np.asarray(counts, dtype=float))
    return logged - logged.mean(axis=1, keepdims=True)


def fourier_weights(n_genes: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Angles and normalized first-Fourier cosine/sine weights."""
    theta = 2 * np.pi * np.arange(n_genes) / n_genes
    scale = np.sqrt(2 / n_genes)
    return theta, scale * np.cos(theta), scale * np.sin(theta)


def objective(covariance: np.ndarray, order: np.ndarray) -> float:
    """Variance captured by the first Fourier mode for one gene order."""
    _, cosine, sine = fourier_weights(len(order))
    ordered = covariance[np.ix_(order, order)]
    return float(cosine @ ordered @ cosine + sine @ ordered @ sine)


def composition_to_xy(composition: np.ndarray) -> np.ndarray:
    """Map a three-part composition to an equilateral simplex."""
    vertices = np.array([[0.5, np.sqrt(3) / 2], [0, 0], [1, 0]])
    composition = np.asarray(composition, dtype=float)
    return (composition / composition.sum()) @ vertices


def workflow_values() -> dict[str, np.ndarray]:
    """Compute every numeric quantity displayed in the schematic."""
    counts = np.array(
        [
            [10.0, 20.0, 30.0],
            [100.0, 200.0, 300.0],
            [10.0, 20.0, 300.0],
        ]
    )
    compositions = counts / counts.sum(axis=1, keepdims=True)
    coordinates = clr(counts)
    covariance = np.cov(coordinates, rowvar=False)

    # All three-gene CLR orders tie; retain the authored non-trivial example.
    display_order = np.array([0, 2, 1])
    scored = [
        (np.asarray(order), objective(covariance, np.asarray(order)))
        for order in permutations(range(3))
    ]
    _, cosine, sine = fourier_weights(3)
    centered = coordinates - coordinates.mean(axis=0, keepdims=True)
    ordered = centered[:, display_order]
    projected = np.column_stack((ordered @ cosine, ordered @ sine))
    projected *= 0.78 / np.hypot(projected[:, 0], projected[:, 1]).max()

    return {
        "counts": counts,
        "compositions": compositions,
        "coordinates": coordinates,
        "covariance": covariance,
        "display_order": display_order,
        "projected": projected,
        "objectives": np.array([value for _, value in scored]),
    }


def panel_title(axis: plt.Axes, number: int, title: str) -> None:
    """Add a numbered workflow title above a panel."""
    axis.text(
        0.5,
        1.045,
        f"{number}. {title}",
        transform=axis.transAxes,
        ha="center",
        va="bottom",
        fontsize=24,
        fontweight="bold",
    )


def draw_composition(axis: plt.Axes, values: dict[str, np.ndarray]) -> None:
    """Draw expression as compositions and introduce CLR coordinates."""
    panel_title(axis, 1, "Gene expression in CLR space")
    top = np.array([0.5, np.sqrt(3) / 2])
    left, right = np.array([0.0, 0.0]), np.array([1.0, 0.0])
    triangle = np.vstack((left, right, top, left))
    axis.plot(triangle[:, 0], triangle[:, 1], color=BLUE, lw=3)

    for fraction in (0.2, 0.4, 0.6, 0.8):
        pairs = (
            ((fraction, 1 - fraction, 0), (fraction, 0, 1 - fraction)),
            ((1 - fraction, fraction, 0), (0, fraction, 1 - fraction)),
            ((1 - fraction, 0, fraction), (0, 1 - fraction, fraction)),
        )
        for pair in pairs:
            line = np.array([composition_to_xy(part) for part in pair])
            axis.plot(line[:, 0], line[:, 1], color=GRAY, lw=1, ls="--", alpha=0.48)

    points = np.array([composition_to_xy(row) for row in values["compositions"]])
    styles = ((BLUE, "o", 120), (ORANGE, "x", 155), (GREEN, "^", 135))
    cell_labels = tuple(
        f"Cell {index}:  ({', '.join(map(str, row.astype(int)))})"
        for index, row in enumerate(values["counts"], start=1)
    )
    for point, (color, marker, size), label in zip(
        points, styles, cell_labels
    ):
        axis.scatter(
            *point,
            color=color,
            marker=marker,
            s=size,
            linewidths=3 if marker == "x" else None,
            label=label,
            zorder=5,
        )

    axis.text(
        0.5,
        0.94,
        "Gene A",
        ha="center",
        va="bottom",
        fontsize=19,
        fontweight="bold",
    )
    axis.text(
        -0.03,
        -0.055,
        "Gene B",
        ha="right",
        va="top",
        fontsize=19,
        fontweight="bold",
    )
    axis.text(
        1.03,
        -0.055,
        "Gene C",
        ha="left",
        va="top",
        fontsize=19,
        fontweight="bold",
    )
    axis.annotate(
        "same relative\nexpression",
        xy=points[0],
        xytext=(0.25, 0.50),
        ha="right",
        va="center",
        fontsize=15,
        arrowprops=dict(arrowstyle="->", lw=1.5),
    )
    axis.legend(
        frameon=False,
        fontsize=16,
        loc="center left",
        bbox_to_anchor=(0.51, 0.80),
    )

    axis.text(
        0.73,
        0.50,
        r"$z_{ig}=\log\frac{x_{ig}}{g(\mathbf{x}_i)}$",
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=27,
    )
    axis.text(
        0.73,
        0.28,
        r"$g(\mathbf{x}_i)=\left(\prod_g\, x_{ig}\right)^{1/G}$",
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=24,
    )
    axis.set(xlim=(-0.18, 2.05), ylim=(-0.18, 1.02), aspect="equal")
    axis.axis("off")


def draw_covariance(axis: plt.Axes, covariance: np.ndarray, figure: plt.Figure) -> None:
    """Draw the CLR covariance matrix used by the ordering objective."""
    panel_title(axis, 2, "Gene covariance")
    matrix_axis = axis.inset_axes([0.02, 0.07, 0.72, 0.79])
    limit = np.abs(covariance).max()
    image = matrix_axis.imshow(covariance, cmap="coolwarm", vmin=-limit, vmax=limit)
    for row in range(3):
        for column in range(3):
            label = (
                f"Var({'ABC'[row]})"
                if row == column
                else f"Cov({'ABC'[row]},{'ABC'[column]})"
            )
            color = "white" if abs(covariance[row, column]) > 0.55 * limit else "black"
            matrix_axis.text(
                column, row, label, ha="center", va="center", fontsize=14, color=color
            )
    matrix_axis.set_xticks(range(3), list("ABC"))
    matrix_axis.set_yticks(range(3), list("ABC"))
    matrix_axis.tick_params(labelsize=17)
    for spine in matrix_axis.spines.values():
        spine.set_visible(False)
    axis.text(
        0.38,
        0.90,
        r"$\Sigma_{\mathrm{genes}}$ in CLR space",
        transform=axis.transAxes,
        ha="center",
        fontsize=20,
        fontweight="bold",
    )
    colorbar_axis = axis.inset_axes([0.73, 0.14, 0.038, 0.66])
    colorbar = figure.colorbar(image, cax=colorbar_axis)
    colorbar.set_label("CLR covariance", fontsize=17)
    colorbar.ax.tick_params(labelsize=15)
    axis.axis("off")


def draw_optimize(axis: plt.Axes) -> None:
    """Illustrate evaluation of candidate circular feature orders."""
    panel_title(axis, 3, "Optimize gene order")
    axis.text(
        0.5,
        0.84,
        r"$F(\pi)=\mathbf{c}^{T}\Sigma_{\pi}\mathbf{c}"
        r"+\mathbf{s}^{T}\Sigma_{\pi}\mathbf{s}$",
        transform=axis.transAxes,
        ha="center",
        fontsize=28,
        fontweight="bold",
    )
    axis.text(
        0.5,
        0.68,
        r"$\pi$: candidate gene order",
        transform=axis.transAxes,
        ha="center",
        fontsize=17,
    )
    axis.text(
        0.5,
        0.59,
        r"$\Sigma_\pi$: CLR covariance reordered by $\pi$",
        transform=axis.transAxes,
        ha="center",
        fontsize=17,
    )
    axis.text(
        0.5,
        0.46,
        r"$\theta_g=2\pi(g-1)/G$",
        transform=axis.transAxes,
        ha="center",
        fontsize=22,
    )
    axis.text(
        0.5,
        0.32,
        r"$c_g=\sqrt{\frac{2}{G}}\cos\theta_g,\qquad"
        r"s_g=\sqrt{\frac{2}{G}}\sin\theta_g$",
        transform=axis.transAxes,
        ha="center",
        fontsize=22,
    )
    axis.text(
        0.5,
        0.16,
        r"$\hat{\pi}=\arg\max_{\pi}\;F(\pi)$",
        transform=axis.transAxes,
        ha="center",
        fontsize=28,
        color=BLUE,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=.22", fc="white", ec=BLUE, lw=2),
    )
    axis.text(
        0.5,
        0.025,
        r"$\hat{\pi}$ is the optimal gene order",
        transform=axis.transAxes,
        ha="center",
        fontsize=17,
    )
    axis.axis("off")


def draw_gene_circle(
    axis: plt.Axes, order: np.ndarray, *, projected: np.ndarray | None = None
) -> None:
    """Draw a learned circular feature ordering."""
    theta, _, _ = fourier_weights(3)
    grid = np.linspace(0, 2 * np.pi, 500)
    axis.plot(np.cos(grid), np.sin(grid), color=BLUE, lw=3)
    axis.axhline(0, color=LIGHT_GRAY, lw=1)
    axis.axvline(0, color=LIGHT_GRAY, lw=1)
    axis.annotate(
        "",
        xy=(1.25, 0),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="->", lw=1.5, color="black"),
    )
    axis.annotate(
        "",
        xy=(0, 1.21),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="->", lw=1.5, color="black"),
    )

    gene_x, gene_y = np.cos(theta), np.sin(theta)
    axis.scatter(gene_x, gene_y, color=BLUE, s=100, zorder=5)
    positions = (
        (1.07, 0.08, "left", "bottom"),
        (-0.58, 0.96, "right", "bottom"),
        (-0.58, -0.96, "right", "top"),
    )
    for position, gene_index in enumerate(order):
        x, y, horizontal, vertical = positions[position]
        axis.text(x, y, GENES[gene_index], fontsize=20, ha=horizontal, va=vertical)

    if projected is None:
        angle = theta[1]
        axis.plot([0, gene_x[1]], [0, gene_y[1]], color=BLUE, lw=2.8)
        axis.scatter(
            gene_x[1],
            gene_y[1],
            s=210,
            facecolor="white",
            edgecolor=BLUE,
            lw=2.8,
            zorder=7,
        )
        arc = np.linspace(0, angle, 120)
        axis.plot(0.25 * np.cos(arc), 0.25 * np.sin(arc), color=BLUE, lw=2)
        axis.text(
            0.36 * np.cos(angle / 2),
            0.36 * np.sin(angle / 2),
            r"$\theta_C$",
            fontsize=20,
            color=BLUE,
            ha="center",
            va="center",
        )
        axis.text(1.14, -0.14, r"$\cos\theta$", fontsize=18)
        axis.text(0.05, 1.16, r"$\sin\theta$", fontsize=18)
    else:
        styles = ((BLUE, "o", 145), (ORANGE, "x", 145), (GREEN, "^", 135))
        for point, (color, marker, size), label in zip(
            projected, styles, ("Cell 1", "Cell 2", "Cell 3")
        ):
            axis.plot([0, point[0]], [0, point[1]], color=color, lw=1.5)
            axis.scatter(
                *point,
                color=color,
                marker=marker,
                s=size,
                linewidths=2.5 if marker == "x" else None,
                label=label,
                zorder=7,
            )
        axis.legend(frameon=False, fontsize=15, loc="upper right")
        axis.annotate(
            "Cells 1 & 2\noverlap",
            xy=projected[0],
            xytext=(projected[0, 0] - 0.48, projected[0, 1] - 0.17),
            fontsize=14,
            ha="center",
            va="top",
            arrowprops=dict(arrowstyle="->", lw=1, color=GRAY),
        )
        point = projected[2]
        phi = np.arctan2(point[1], point[0])
        arc = np.linspace(0, phi, 100)
        axis.plot(0.22 * np.cos(arc), 0.22 * np.sin(arc), color=GREEN, lw=1.8)
        axis.text(
            0.31 * np.cos(phi / 2),
            0.31 * np.sin(phi / 2),
            r"$\phi_i$",
            fontsize=17,
            color=GREEN,
            ha="center",
            va="center",
        )
        axis.text(
            0.57 * point[0],
            0.57 * point[1] + 0.04,
            r"$r_i$",
            fontsize=17,
            color=GREEN,
            ha="center",
        )
        axis.text(1.14, -0.14, r"$u_i$", fontsize=17)
        axis.text(0.05, 1.16, r"$v_i$", fontsize=17)

    axis.set(xlim=(-1.30, 1.40), ylim=(-1.22, 1.28), aspect="equal")
    axis.axis("off")


def draw_fourier(axis: plt.Axes, order: np.ndarray) -> None:
    """Show the first circular Fourier weights assigned to ordered genes."""
    panel_title(axis, 4, "First Fourier mode")
    draw_gene_circle(axis, order)
    axis.text(
        0.5,
        -0.07,
        r"$\theta_C=\frac{2\pi(2-1)}{3}=\frac{2\pi}{3}$",
        transform=axis.transAxes,
        ha="center",
        fontsize=20,
    )


def draw_projection(axis: plt.Axes, order: np.ndarray, projected: np.ndarray) -> None:
    """Draw example samples in the learned two-dimensional projection."""
    panel_title(axis, 5, "Project cells")
    draw_gene_circle(axis, order, projected=projected)
    axis.text(
        0.5,
        -0.08,
        r"$u_i=\sum_g Z_{i,\widehat{\pi}(g)}c_g$"
        "   "
        r"$v_i=\sum_g Z_{i,\widehat{\pi}(g)}s_g$",
        transform=axis.transAxes,
        ha="center",
        fontsize=18,
    )


def block_arrow(
    figure: plt.Figure,
    source: plt.Axes,
    target: plt.Axes,
    *,
    downward: bool = False,
    horizontal_offset: float = 0.0,
) -> None:
    """Draw one uniformly sized native block transition arrow."""
    first, second = source.get_position(), target.get_position()
    horizontal_length = 0.032
    if downward:
        x = (first.x0 + first.x1) / 2
        center = (first.y0 + second.y1) / 2 + 0.01
        vertical_length = horizontal_length * figure.get_figwidth() / figure.get_figheight()
        start = (x, center + vertical_length / 2)
        end = (x, center - vertical_length / 2)
    else:
        y = ((first.y0 + first.y1) + (second.y0 + second.y1)) / 4
        if second.x0 > first.x0:
            center = (first.x1 + second.x0) / 2 + horizontal_offset
            start = (center - horizontal_length / 2, y)
            end = (center + horizontal_length / 2, y)
        else:
            center = (first.x0 + second.x1) / 2 + horizontal_offset
            start = (center + horizontal_length / 2, y)
            end = (center - horizontal_length / 2, y)
    figure.add_artist(
        FancyArrowPatch(
            start,
            end,
            transform=figure.transFigure,
            arrowstyle="-|>",
            mutation_scale=22,
            linewidth=3.2,
            color=GRAY,
            zorder=20,
        )
    )


def save_individual_panel(figure: plt.Figure, filename: str) -> None:
    """Save one standalone methods panel in each publication output format."""
    METHODS_PANEL_DIR.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        METHODS_PANEL_DIR / f"{filename}.pdf", bbox_inches="tight", pad_inches=0.15
    )
    figure.savefig(
        METHODS_PANEL_DIR / f"{filename}.svg", bbox_inches="tight", pad_inches=0.15
    )
    figure.savefig(
        METHODS_PANEL_DIR / f"{filename}.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.15,
    )
    plt.close(figure)


def render_individual_panels(values: dict[str, np.ndarray]) -> None:
    """Render standalone copies of workflow panels 1, 4, and 5."""
    composition_figure, composition_axis = plt.subplots(figsize=(10, 5.8))
    composition_figure.subplots_adjust(left=0.04, right=0.96, bottom=0.08, top=0.86)
    draw_composition(composition_axis, values)
    save_individual_panel(composition_figure, "Panel_1_Gene_expression_CLR")

    fourier_figure, fourier_axis = plt.subplots(figsize=(7.0, 6.0))
    fourier_figure.subplots_adjust(left=0.05, right=0.95, bottom=0.16, top=0.86)
    draw_fourier(fourier_axis, values["display_order"])
    save_individual_panel(fourier_figure, "Panel_4_First_Fourier_mode")

    projection_figure, projection_axis = plt.subplots(figsize=(7.0, 6.0))
    projection_figure.subplots_adjust(left=0.05, right=0.95, bottom=0.16, top=0.86)
    draw_projection(
        projection_axis, values["display_order"], values["projected"]
    )
    save_individual_panel(projection_figure, "Panel_5_Project_cells")


def render_figure() -> None:
    """Assemble and save the complete standalone methods figure."""
    values = workflow_values()
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 18,
            "axes.titlesize": 24,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    figure = plt.figure(figsize=(20, 12.8), facecolor="white")
    grid = figure.add_gridspec(
        2,
        6,
        left=0.05,
        right=0.95,
        bottom=0.085,
        top=0.88,
        wspace=0.26,
        hspace=0.36,
    )
    axes = [
        figure.add_subplot(grid[row, start:stop])
        for row, start, stop in (
            (0, 0, 3),
            (0, 3, 6),
            (1, 4, 6),
            (1, 2, 4),
            (1, 0, 2),
        )
    ]

    draw_composition(axes[0], values)
    draw_covariance(axes[1], values["covariance"], figure)
    draw_optimize(axes[2])
    draw_fourier(axes[3], values["display_order"])
    draw_projection(axes[4], values["display_order"], values["projected"])

    block_arrow(figure, axes[0], axes[1])
    block_arrow(figure, axes[1], axes[2], downward=True)
    block_arrow(figure, axes[2], axes[3], horizontal_offset=0.018)
    block_arrow(figure, axes[3], axes[4])

    figure.suptitle("clrcycle workflow", fontsize=30, fontweight="bold", y=0.97)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURE_DIR / "Methods_Figure.pdf")
    figure.savefig(FIGURE_DIR / "Methods_Figure.svg")
    figure.savefig(FIGURE_DIR / "Methods_Figure.png", dpi=300)
    plt.close(figure)
    render_individual_panels(values)


if __name__ == "__main__":
    render_figure()
    print("Wrote Methods Figure PDF/SVG/PNG from one standalone native script")
