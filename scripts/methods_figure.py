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
        fontsize=22,
        fontweight="bold",
    )


def draw_composition(axis: plt.Axes, values: dict[str, np.ndarray]) -> None:
    """Draw the input compositions in simplex and stacked-bar views."""
    panel_title(axis, 1, "Composition")
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
    for point, (color, marker, size), label in zip(
        points, styles, ("Cell 1", "Cell 2", "Cell 3")
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
        "Gene A\n100%",
        ha="center",
        va="bottom",
        fontsize=15,
        fontweight="bold",
    )
    axis.text(
        -0.03,
        -0.055,
        "Gene B\n100%",
        ha="right",
        va="top",
        fontsize=15,
        fontweight="bold",
    )
    axis.text(
        1.03,
        -0.055,
        "Gene C\n100%",
        ha="left",
        va="top",
        fontsize=15,
        fontweight="bold",
    )
    axis.annotate(
        "same relative\nexpression\n"
        r"$\left(\frac{1}{6},\frac{1}{3},\frac{1}{2}\right)$",
        xy=points[0],
        xytext=(0.17, 0.50),
        ha="right",
        va="center",
        fontsize=12,
        arrowprops=dict(arrowstyle="->", lw=1.5),
    )
    axis.legend(frameon=False, fontsize=12, loc="upper right")

    table = axis.table(
        cellText=values["counts"].astype(int),
        rowLabels=("Cell 1", "Cell 2", "Cell 3"),
        colLabels=("A", "B", "C"),
        cellLoc="center",
        rowLoc="center",
        bbox=[0.26, 0.00, 0.58, 0.14],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor(".7")
        cell.set_linewidth(0.6)
        if row == 0:
            cell.set_facecolor(".93")
            cell.set_text_props(fontweight="bold")

    axis.set(xlim=(-0.20, 1.20), ylim=(-0.34, 1.04), aspect="equal")
    axis.axis("off")


def draw_clr(axis: plt.Axes) -> None:
    """Explain the centered log-ratio transform graphically."""
    panel_title(axis, 2, "CLR coordinates")
    axis.text(
        0.5,
        0.70,
        r"$z_{ig}=\log\frac{x_{ig}}{g(\mathbf{x}_i)}$",
        transform=axis.transAxes,
        ha="center",
        fontsize=25,
    )
    axis.text(
        0.5,
        0.49,
        r"gene expression $x_{ig}$ relative to the"
        "\n"
        r"cell geometric mean $g(\mathbf{x}_i)$",
        transform=axis.transAxes,
        ha="center",
        fontsize=17,
    )
    axis.text(
        0.5,
        0.28,
        r"$\mathrm{clr}(c\mathbf{x})=\mathrm{clr}(\mathbf{x})$",
        transform=axis.transAxes,
        ha="center",
        fontsize=22,
        bbox=dict(boxstyle="round,pad=.32", fc="white", ec=BLUE, lw=2),
    )
    axis.text(
        0.5,
        0.10,
        "scale invariant",
        transform=axis.transAxes,
        ha="center",
        fontsize=16,
        fontweight="bold",
    )
    axis.axis("off")


def draw_covariance(axis: plt.Axes, covariance: np.ndarray, figure: plt.Figure) -> None:
    """Draw the CLR covariance matrix used by the ordering objective."""
    panel_title(axis, 3, "Gene covariance")
    matrix_axis = axis.inset_axes([0.04, 0.12, 0.80, 0.70])
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
                column, row, label, ha="center", va="center", fontsize=10.5, color=color
            )
    matrix_axis.set_xticks(range(3), list("ABC"))
    matrix_axis.set_yticks(range(3), list("ABC"))
    matrix_axis.tick_params(labelsize=13)
    for spine in matrix_axis.spines.values():
        spine.set_visible(False)
    axis.text(
        0.44,
        0.86,
        r"$\Sigma_{\mathrm{genes}}$ in CLR space",
        transform=axis.transAxes,
        ha="center",
        fontsize=16,
        fontweight="bold",
    )
    colorbar_axis = axis.inset_axes([0.88, 0.18, 0.035, 0.58])
    colorbar = figure.colorbar(image, cax=colorbar_axis)
    colorbar.set_label("CLR covariance", fontsize=13)
    colorbar.ax.tick_params(labelsize=11)
    axis.axis("off")


def draw_optimize(axis: plt.Axes) -> None:
    """Illustrate evaluation of candidate circular feature orders."""
    panel_title(axis, 4, "Optimize gene order")
    axis.text(
        0.5,
        0.82,
        r"$F(\pi)=\mathbf{c}^{T}\Sigma_{\pi}\mathbf{c}"
        r"+\mathbf{s}^{T}\Sigma_{\pi}\mathbf{s}$",
        transform=axis.transAxes,
        ha="center",
        fontsize=22,
        fontweight="bold",
    )
    axis.text(
        0.5,
        0.66,
        r"$\pi$: candidate gene order"
        r"$\qquad\Sigma_\pi$: CLR covariance reordered by $\pi$",
        transform=axis.transAxes,
        ha="center",
        fontsize=13,
    )
    axis.text(
        0.5,
        0.50,
        r"$\theta_g=2\pi(g-1)/G$",
        transform=axis.transAxes,
        ha="center",
        fontsize=16,
    )
    axis.text(
        0.5,
        0.34,
        r"$c_g=\sqrt{\frac{2}{G}}\cos\theta_g,\qquad"
        r"s_g=\sqrt{\frac{2}{G}}\sin\theta_g$",
        transform=axis.transAxes,
        ha="center",
        fontsize=17,
    )
    axis.text(
        0.5,
        0.16,
        r"$\hat{\pi}=\arg\max_{\pi}\;F(\pi)$",
        transform=axis.transAxes,
        ha="center",
        fontsize=22,
        color=BLUE,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=.22", fc="white", ec=BLUE, lw=2),
    )
    axis.text(
        0.5,
        0.035,
        r"$\hat{\pi}$ is the optimal gene order",
        transform=axis.transAxes,
        ha="center",
        fontsize=13,
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
        xy=(1.28, 0),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="->", lw=1.5, color="black"),
    )
    axis.annotate(
        "",
        xy=(0, 1.24),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="->", lw=1.5, color="black"),
    )

    gene_x, gene_y = np.cos(theta), np.sin(theta)
    axis.scatter(gene_x, gene_y, color=BLUE, s=100, zorder=5)
    positions = (
        (1.10, 0.10, "left", "bottom"),
        (-0.62, 0.98, "right", "bottom"),
        (-0.62, -0.98, "right", "top"),
    )
    for position, gene_index in enumerate(order):
        x, y, horizontal, vertical = positions[position]
        axis.text(x, y, GENES[gene_index], fontsize=15, ha=horizontal, va=vertical)

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
            fontsize=15,
            color=BLUE,
            ha="center",
            va="center",
        )
        axis.text(1.23, -0.15, r"$\cos\theta$", fontsize=14)
        axis.text(0.06, 1.23, r"$\sin\theta$", fontsize=14)
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
        axis.legend(frameon=False, fontsize=11, loc="upper right")
        axis.annotate(
            "Cells 1 & 2\noverlap",
            xy=projected[0],
            xytext=(projected[0, 0] - 0.48, projected[0, 1] - 0.17),
            fontsize=10,
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
            fontsize=13,
            color=GREEN,
            ha="center",
            va="center",
        )
        axis.text(
            0.57 * point[0],
            0.57 * point[1] + 0.04,
            r"$r_i$",
            fontsize=13,
            color=GREEN,
            ha="center",
        )
        axis.text(1.23, -0.15, r"$u_i$", fontsize=13)
        axis.text(0.06, 1.22, r"$v_i$", fontsize=13)

    axis.set(xlim=(-1.38, 1.48), ylim=(-1.30, 1.35), aspect="equal")
    axis.axis("off")


def draw_fourier(axis: plt.Axes, order: np.ndarray) -> None:
    """Show the first circular Fourier weights assigned to ordered genes."""
    panel_title(axis, 5, "First Fourier mode")
    draw_gene_circle(axis, order)
    axis.text(
        0.5,
        -0.07,
        r"$\theta_C=\frac{2\pi(2-1)}{3}=\frac{2\pi}{3}$",
        transform=axis.transAxes,
        ha="center",
        fontsize=16,
    )


def draw_projection(axis: plt.Axes, order: np.ndarray, projected: np.ndarray) -> None:
    """Draw example samples in the learned two-dimensional projection."""
    panel_title(axis, 6, "Project cells")
    draw_gene_circle(axis, order, projected=projected)
    axis.text(
        0.5,
        -0.08,
        r"$u_i=\sum_g Z_{i,\widehat{\pi}(g)}c_g$"
        "   "
        r"$v_i=\sum_g Z_{i,\widehat{\pi}(g)}s_g$",
        transform=axis.transAxes,
        ha="center",
        fontsize=12,
    )


def block_arrow(
    figure: plt.Figure, source: plt.Axes, target: plt.Axes, *, downward: bool = False
) -> None:
    """Draw one native, filled-looking block transition arrow."""
    first, second = source.get_position(), target.get_position()
    if downward:
        x = (first.x0 + first.x1) / 2
        start, end = (x, first.y0 - 0.012), (x, second.y1 + 0.055)
    else:
        y = ((first.y0 + first.y1) + (second.y0 + second.y1)) / 4
        if second.x0 > first.x0:
            start, end = (first.x1 + 0.012, y), (second.x0 - 0.012, y)
        else:
            start, end = (first.x0 - 0.012, y), (second.x1 + 0.012, y)
    figure.add_artist(
        FancyArrowPatch(
            start,
            end,
            transform=figure.transFigure,
            arrowstyle="-|>",
            mutation_scale=28,
            linewidth=4,
            color=GRAY,
            zorder=20,
        )
    )


def render_figure() -> None:
    """Assemble and save the complete standalone methods figure."""
    values = workflow_values()
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 16,
            "axes.titlesize": 22,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    figure = plt.figure(figsize=(20, 12.8), facecolor="white")
    grid = figure.add_gridspec(
        2,
        3,
        left=0.04,
        right=0.96,
        bottom=0.07,
        top=0.87,
        wspace=0.34,
        hspace=0.40,
    )
    axes = [
        figure.add_subplot(grid[row, column])
        for row, column in ((0, 0), (0, 1), (0, 2), (1, 2), (1, 1), (1, 0))
    ]

    draw_composition(axes[0], values)
    draw_clr(axes[1])
    draw_covariance(axes[2], values["covariance"], figure)
    draw_optimize(axes[3])
    draw_fourier(axes[4], values["display_order"])
    draw_projection(axes[5], values["display_order"], values["projected"])

    block_arrow(figure, axes[0], axes[1])
    block_arrow(figure, axes[1], axes[2])
    block_arrow(figure, axes[2], axes[3], downward=True)
    block_arrow(figure, axes[3], axes[4])
    block_arrow(figure, axes[4], axes[5])

    figure.suptitle("clrcycle workflow", fontsize=26, fontweight="bold", y=0.97)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURE_DIR / "Methods_Figure.pdf", bbox_inches="tight")
    figure.savefig(FIGURE_DIR / "Methods_Figure.svg", bbox_inches="tight")
    figure.savefig(FIGURE_DIR / "Methods_Figure.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    render_figure()
    print("Wrote Methods Figure PDF/SVG/PNG from one standalone native script")
