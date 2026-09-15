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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "figures"
METHODS_PANEL_DIR = FIGURE_DIR / "methods figs"

INK = "#303030"
GRAY = "0.42"
LIGHT_GRAY = "0.88"
GENES = np.array(["Gene A", "Gene B", "Gene C"])
GENE_COLORS = np.array(["#c76b16", "#8b59a5", "#008577"])
GRAY_CMAP = LinearSegmentedColormap.from_list(
    "neutral_gray", ["#f3f3f3", "#b0b0b0", "#303030"]
)


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
    # Ten phases on a rotated ellipse. The two major-axis endpoints touch
    # the unit gene circle; the shorter radii illustrate variable amplitude.
    # Three-gene CLR profiles remain entirely first harmonic.
    display_order = np.array([0, 2, 1])
    _, cosine, sine = fourier_weights(3)
    phases = 2 * np.pi * np.arange(10) / 10
    rotation = np.pi / 4
    rotation_matrix = np.array([[np.cos(rotation), -np.sin(rotation)],
                                [np.sin(rotation), np.cos(rotation)]])
    scores = np.column_stack((np.cos(phases), 0.58 * np.sin(phases))) @ rotation_matrix.T
    coordinates = np.empty((len(scores), 3))
    coordinates[:, display_order] = np.outer(scores[:, 0], cosine) + np.outer(scores[:, 1], sine)
    # Positive expression, with geometric mean 100 in every cell. Keep full
    # precision for all computations; only heatmap annotations are rounded.
    counts = 100 * np.exp(coordinates)
    compositions = counts / counts.sum(axis=1, keepdims=True)
    coordinates = clr(counts)
    covariance = np.cov(coordinates, rowvar=False)
    scored = [
        (np.asarray(order), objective(covariance, np.asarray(order)))
        for order in permutations(range(3))
    ]
    centered = coordinates - coordinates.mean(axis=0, keepdims=True)
    ordered = centered[:, display_order]
    projected = np.column_stack((ordered @ cosine, ordered @ sine))
    projected /= np.linalg.norm(projected, axis=1).max()

    return {
        "counts": counts,
        "compositions": compositions,
        "coordinates": coordinates,
        "covariance": covariance,
        "display_order": display_order,
        "projected": projected,
        "objectives": np.array([value for _, value in scored]),
    }


def panel_title(
    axis: plt.Axes,
    number: int,
    title: str,
    *,
    horizontal_position: float = 0.5,
    vertical_position: float = 1.045,
) -> None:
    """Add a numbered workflow title above a panel."""
    axis.text(
        horizontal_position,
        vertical_position,
        f"{number}. {title}",
        transform=axis.transAxes,
        ha="center",
        va="bottom",
        fontsize=24,
        fontweight="bold",
    )


def draw_composition(
    axis: plt.Axes, values: dict[str, np.ndarray], *, show_title: bool = True
) -> None:
    """Draw expression as compositions and introduce CLR coordinates."""
    if show_title:
        panel_title(axis, 1, "Gene expression in CLR space", horizontal_position=0.54)
    top = np.array([0.5, np.sqrt(3) / 2])
    left, right = np.array([0.0, 0.0]), np.array([1.0, 0.0])
    triangle = np.vstack((left, right, top, left))
    axis.plot(triangle[:, 0], triangle[:, 1], color=INK, lw=3)
    vertices = np.array([top, left, right])
    axis.scatter(vertices[:, 0], vertices[:, 1], c=GENE_COLORS,
                 s=110, zorder=6, edgecolors="white", linewidths=1)

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
    for index, point in enumerate(points):
        color = INK
        axis.scatter(*point, color=color, marker="o", s=75, zorder=5)
        axis.annotate(str(index + 1), point, xytext=(5, 5),
                      textcoords="offset points", fontsize=12, color=color)

    axis.text(
        0.5,
        0.94,
        "Gene A",
        color=GENE_COLORS[0],
        ha="center",
        va="bottom",
        fontsize=19,
        fontweight="bold",
    )
    axis.text(
        -0.03,
        -0.055,
        "Gene B",
        color=GENE_COLORS[1],
        ha="right",
        va="top",
        fontsize=19,
        fontweight="bold",
    )
    axis.text(
        1.03,
        -0.055,
        "Gene C",
        color=GENE_COLORS[2],
        ha="left",
        va="top",
        fontsize=19,
        fontweight="bold",
    )
    axis.text(
        0.73,
        0.62,
        r"$z_{ig}=\log\frac{x_{ig}}{g(\mathbf{x}_i)}$",
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=23,
    )
    axis.text(
        0.73,
        0.40,
        r"$g(\mathbf{x}_i)=\left(\prod_g\, x_{ig}\right)^{1/G}$",
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=20,
    )
    axis.set(xlim=(-0.18, 2.05), ylim=(-0.18, 1.02), aspect="equal")
    axis.axis("off")


def draw_expression_heatmap(axis: plt.Axes, values: dict[str, np.ndarray]) -> None:
    """Show normalized expression fractions with the same cell IDs used in panels 1 and 5."""
    # Display per-cell fractions directly, without a logarithm.
    expression_fractions = values["compositions"]
    limit = expression_fractions.max()
    image = axis.imshow(expression_fractions, cmap=GRAY_CMAP, vmin=0, vmax=limit, aspect="auto")
    colorbar_axis = axis.inset_axes([1.05, 0.0, 0.045, 1.0])
    colorbar = axis.figure.colorbar(image, cax=colorbar_axis)
    colorbar.set_label("Relative expression", fontsize=12)
    colorbar.ax.tick_params(labelsize=11)
    for row, expression in enumerate(expression_fractions):
        for column, value in enumerate(expression):
            color = "white" if value > 0.65 * limit else "black"
            axis.text(column, row, f"{value:.2f}", ha="center", va="center",
                      fontsize=12, color=color)
    axis.set_xticks(range(3), GENES, fontsize=13)
    axis.xaxis.tick_top()
    for tick, color in zip(axis.get_xticklabels(), GENE_COLORS):
        tick.set_color(color)
        tick.set_fontweight("bold")
    axis.set_yticks(range(len(expression_fractions)),
                   [f"Cell {i + 1}" for i in range(len(expression_fractions))], fontsize=12)
    for tick in axis.get_yticklabels():
        tick.set_color(INK)
    axis.tick_params(length=0)
    for spine in axis.spines.values():
        spine.set_visible(False)



def draw_covariance(axis: plt.Axes, covariance: np.ndarray, figure: plt.Figure) -> None:
    """Draw the CLR covariance matrix used by the ordering objective."""
    panel_title(axis, 2, "Gene covariance", horizontal_position=0.38)
    matrix_axis = axis.inset_axes([0.02, 0.07, 0.80, 0.79])
    limit = np.abs(covariance).max()
    image = matrix_axis.imshow(covariance, cmap=GRAY_CMAP, vmin=-limit, vmax=limit)
    for row in range(3):
        for column in range(3):
            label = (
                f"Var({'ABC'[row]})"
                if row == column
                else f"Cov({'ABC'[row]},{'ABC'[column]})"
            )
            color = "white" if covariance[row, column] > 0.1 * limit else "black"
            matrix_axis.text(
                column, row, label, ha="center", va="center", fontsize=11, color=color
            )
    matrix_axis.set_xticks(range(3), list("ABC"))
    matrix_axis.set_yticks(range(3), list("ABC"))
    matrix_axis.tick_params(labelsize=17)
    for ticks in (matrix_axis.get_xticklabels(), matrix_axis.get_yticklabels()):
        for tick, color in zip(ticks, GENE_COLORS):
            tick.set_color(color)
            tick.set_fontweight("bold")
    for spine in matrix_axis.spines.values():
        spine.set_visible(False)
    axis.text(
        0.42,
        0.90,
        r"$\Sigma_{\mathrm{genes}}$ in CLR space",
        transform=axis.transAxes,
        ha="center",
        fontsize=20,
        fontweight="bold",
    )
    colorbar_axis = axis.inset_axes([0.86, 0.19, 0.038, 0.55])
    colorbar = figure.colorbar(image, cax=colorbar_axis)
    colorbar.set_label("CLR covariance", fontsize=17)
    colorbar.ax.tick_params(labelsize=15)
    axis.axis("off")


def draw_optimize(axis: plt.Axes) -> None:
    """Illustrate evaluation of candidate circular feature orders."""
    panel_title(axis, 3, "Maximize first Fourier variance")
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
        color=INK,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=.22", fc="white", ec=INK, lw=2),
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
    axis.plot(np.cos(grid), np.sin(grid), color=INK, lw=3)
    axis.plot([-1.30, 1.40], [0, 0], color=LIGHT_GRAY, lw=1)
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
    axis.scatter(gene_x, gene_y, c=GENE_COLORS[order], s=120, zorder=8)
    positions = (
        (1.07, 0.08, "left", "bottom"),
        (-0.58, 0.96, "right", "bottom"),
        (-0.58, -0.96, "right", "top"),
    )
    for position, gene_index in enumerate(order):
        x, y, horizontal, vertical = positions[position]
        axis.text(x, y, GENES[gene_index], fontsize=20, ha=horizontal, va=vertical,
                  color=GENE_COLORS[gene_index])

    if projected is None:
        angle = theta[1]
        gene_color = GENE_COLORS[order[1]]
        axis.plot([0, gene_x[1]], [0, gene_y[1]], color=gene_color, lw=2.8)
        axis.scatter(
            gene_x[1],
            gene_y[1],
            s=210,
            facecolor="white",
            edgecolor=gene_color,
            lw=2.8,
            zorder=9,
        )
        arc = np.linspace(0, angle, 120)
        axis.plot(0.25 * np.cos(arc), 0.25 * np.sin(arc), color=gene_color, lw=2)
        axis.text(
            0.36 * np.cos(angle / 2),
            0.36 * np.sin(angle / 2),
            r"$\theta_C$",
            fontsize=20,
            color=gene_color,
            ha="center",
            va="center",
        )
        axis.text(1.14, -0.14, r"$\cos\theta$", fontsize=18)
        axis.text(0.05, 1.16, r"$\sin\theta$", fontsize=18)
    else:
        # Fit all ten cells. Uniform ellipse phases make this moment fit exact.
        # For uniformly sampled phases, 2 * covariance gives squared axes.
        reference = projected
        center = reference.mean(axis=0)
        eigenvalues, eigenvectors = np.linalg.eigh(np.cov(reference, rowvar=False, bias=True))
        ellipse = center + (np.column_stack((np.cos(grid), np.sin(grid)))
                            * np.sqrt(2 * eigenvalues)) @ eigenvectors.T
        axis.plot(ellipse[:, 0], ellipse[:, 1], color=GRAY, lw=2,
                  ls="--", label="Ellipse fit", zorder=4)
        ellipse_precision = np.linalg.inv(np.cov(reference, rowvar=False, bias=True))
        for index, point in enumerate(projected):
            color = INK
            axis.scatter(*point, color=color, marker="o",
                         s=100, zorder=7, label="Cells" if index == 0 else None)
            # Place labels inside the fitted ellipse along its inward normal.
            # This clears both the ellipse and the enclosing reference circle.
            offset = -ellipse_precision @ (point - center)
            offset /= np.linalg.norm(offset)
            axis.annotate(str(index + 1), point, xytext=tuple(17 * offset),
                          textcoords="offset points", fontsize=12, color=color,
                          ha="center", va="center")
        handles, labels = axis.get_legend_handles_labels()
        axis.legend(handles[::-1], labels[::-1], frameon=False, fontsize=11,
                    loc="center right", bbox_to_anchor=(-1.16, 0),
                    bbox_transform=axis.transData, borderaxespad=0, labelspacing=0.5)
        point = projected[2]
        axis.plot([0, point[0]], [0, point[1]], color=INK, lw=1.5)
        phi = np.arctan2(point[1], point[0])
        arc = np.linspace(0, phi, 100)
        axis.plot(0.22 * np.cos(arc), 0.22 * np.sin(arc), color=INK, lw=1.8)
        axis.text(
            0.31 * np.cos(phi / 2),
            0.31 * np.sin(phi / 2),
            r"$\phi_3$",
            fontsize=20,
            color=INK,
            ha="center",
            va="center",
        )
        axis.text(
            0.57 * point[0] - 0.13,
            0.57 * point[1] + 0.04,
            r"$r_3$",
            fontsize=20,
            color=INK,
            ha="center",
        )
        axis.text(1.14, -0.14, r"$u_i$", fontsize=20)
        axis.text(0.05, 1.16, r"$v_i$", fontsize=20)

    axis.set(xlim=(-1.68, 1.67) if projected is not None else (-1.30, 1.40),
             ylim=(-1.22, 1.28), aspect="equal")
    axis.axis("off")


def draw_fourier(axis: plt.Axes, order: np.ndarray) -> None:
    """Show the first circular Fourier weights assigned to ordered genes."""
    panel_title(axis, 4, "Gene ordered\ncircular coordinates")
    draw_gene_circle(axis, order)
    axis.text(
        0.5,
        -0.07,
        r"$\theta_C=\frac{2\pi(2-1)}{3}=\frac{2\pi}{3}$",
        transform=axis.transAxes,
        ha="center",
        fontsize=20,
        color=GENE_COLORS[order[1]],
    )


def draw_projection(axis: plt.Axes, order: np.ndarray, projected: np.ndarray) -> None:
    """Draw example samples in the learned two-dimensional projection."""
    panel_title(
        axis,
        5,
        "Project cells onto\njoint feature-sample geometry",
        vertical_position=1.075,
    )
    draw_gene_circle(axis, order, projected=projected)
    axis.text(
        0.5,
        -0.065,
        r"$u_i=\sum_g Z_{i,\widehat{\pi}(g)}c_g$"
        "   "
        r"$v_i=\sum_g Z_{i,\widehat{\pi}(g)}s_g$",
        transform=axis.transAxes,
        ha="center",
        fontsize=24,
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
        source_point = ((first.x0 + first.x1) / 2, first.y0)
        target_point = ((second.x0 + second.x1) / 2, second.y1)
        center = (
            (source_point[0] + target_point[0]) / 2,
            (source_point[1] + target_point[1]) / 2 + 0.01,
        )
        width, height = figure.get_figwidth(), figure.get_figheight()
        direction = np.array(
            [
                (target_point[0] - source_point[0]) * width,
                (target_point[1] - source_point[1]) * height,
            ]
        )
        direction /= np.hypot(*direction)
        half_length = horizontal_length * width / 2
        offset = np.array(
            [half_length * direction[0] / width, half_length * direction[1] / height]
        )
        start = tuple(np.asarray(center) - offset)
        end = tuple(np.asarray(center) + offset)
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
        12,
        left=0.05,
        right=0.95,
        bottom=0.055,
        top=0.88,
        wspace=0.12,
        hspace=0.44,
    )
    axes = [
        figure.add_subplot(grid[row, start:stop])
        for row, start, stop in (
            (0, 4, 9),
            (0, 9, 12),
            (1, 8, 12),
            (1, 4, 8),
            (1, 0, 4),
        )
    ]
    for axis, horizontal_shift in zip(axes[3:5], (-0.018, -0.030)):
        position = axis.get_position()
        axis.set_position(
            [position.x0 + horizontal_shift, position.y0, position.width, position.height]
        )

    # Center the heatmap and its colorbar within the top-left grid cell.
    expression_panel = figure.add_subplot(grid[0, 0:4])
    expression_panel.axis("off")
    for axis in (expression_panel, axes[0]):
        position = axis.get_position()
        axis.set_position([position.x0 - 0.03, position.y0,
                           position.width, position.height])
    heatmap_axis = expression_panel.inset_axes([0.16, 0.05, 0.57, 0.71])
    draw_expression_heatmap(heatmap_axis, values)
    draw_composition(axes[0], values, show_title=False)
    draw_covariance(axes[1], values["covariance"], figure)
    draw_optimize(axes[2])
    draw_fourier(axes[3], values["display_order"])
    draw_projection(axes[4], values["display_order"], values["projected"])

    # Align the number 1 with the number 5 on the centered first title line.
    figure.canvas.draw()
    projection_title = axes[4].texts[0]
    renderer = figure.canvas.get_renderer()
    first_line = projection_title.get_text().split("\n")[0]
    first_line_width, _, _ = renderer.get_text_width_height_descent(
        first_line, projection_title.get_fontproperties(), ismath=False
    )
    title_center_x = projection_title.get_transform().transform(
        projection_title.get_position()
    )[0]
    title_x = figure.transFigure.inverted().transform(
        (title_center_x - first_line_width / 2, 0)
    )[0]
    figure.text(title_x, 0.885, "1. Gene expression in CLR space",
                ha="left", va="bottom", fontsize=24, fontweight="bold")

    block_arrow(figure, axes[0], axes[1], horizontal_offset=-0.030)
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
