"""Reproduce Supplementary Figure 1: circular-ordering benchmarks.

This downstream analysis reads validated Figure 3 artifacts, compares the
clrcycle order with an ElPiGraph cycle, and measures harmonic distortion across
linear, Fourier, and cyclic-spline models. Tables are written to
``results/supp_fig1/``; ``--plot`` additionally renders the figure.
"""

from pathlib import Path
import argparse, json, platform, sys
from importlib.metadata import version
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd
import elpigraph
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler
import clrcycle

ROOT = Path(__file__).resolve().parents[1]
SHARED = ROOT / "results/fig3"
OUT = ROOT / "results/supp_fig1"
OUT.mkdir(parents=True, exist_ok=True)
TWO_PI = 2 * np.pi
ALPHAS = np.logspace(-5, 5, 31)
KNOTS = np.linspace(1.0, 29.0, 8)[:, None]


def align_phase(inferred, reference):
    """Rotate and optionally reflect inferred phases to a reference circle."""
    candidates = []
    for sign in (1, -1):
        oriented = np.mod(sign * inferred, TWO_PI)
        shift = np.angle(np.mean(np.exp(1j * (reference - oriented))))
        aligned = np.mod(oriented + shift, TWO_PI)
        candidates.append(
            (float(abs(np.mean(np.exp(1j * (aligned - reference))))), aligned)
        )
    return max(candidates, key=lambda x: x[0])[::-1]


def cycle_node_order(edges, n):
    """Traverse an undirected ElPiGraph cycle and return its node order."""
    neighbors = {i: [] for i in range(n)}
    for a, b in edges:
        neighbors[int(a)].append(int(b))
        neighbors[int(b)].append(int(a))
    if any(len(v) != 2 for v in neighbors.values()):
        raise ValueError("ElPiGraph output is not one cycle")
    order = [0]
    previous = None
    current = 0
    while True:
        following = [x for x in neighbors[current] if x != previous][0]
        if following == 0:
            break
        order.append(following)
        previous, current = current, following
        if len(order) > n:
            raise ValueError("Failed cycle traversal")
    if len(order) != n:
        raise ValueError("Cycle traversal omitted nodes")
    return np.asarray(order)


def project_to_cycle(points, nodes, edges):
    """Project points to the closest edge of a polygonal cycle."""
    order = cycle_node_order(edges, len(nodes))
    starts, ends = nodes[order], nodes[np.roll(order, -1)]
    vectors = ends - starts
    lengths = np.linalg.norm(vectors, axis=1)
    cumulative = np.r_[0.0, np.cumsum(lengths[:-1])]
    best = np.full(len(points), np.inf)
    phase = np.zeros(len(points))
    projection = np.zeros_like(points)
    for edge, (start, vector, length) in enumerate(zip(starts, vectors, lengths)):
        fraction = np.clip((points - start) @ vector / (length * length), 0, 1)
        candidate = start + fraction[:, None] * vector
        distance = np.square(points - candidate).sum(1)
        update = distance < best
        phase[update] = (
            TWO_PI * (cumulative[edge] + fraction[update] * length) / lengths.sum()
        )
        projection[update] = candidate[update]
        best[update] = distance[update]
    return phase, projection


def circular_ranks(phase):
    """Convert phase values to uniformly spaced circular rank angles."""
    order = np.argsort(phase, kind="stable")
    ranks = np.empty(len(phase), int)
    ranks[order] = np.arange(len(phase))
    return TWO_PI * ranks / len(phase)


def load_shared():
    """Load and validate the Figure 3 artifacts required by this benchmark."""
    paths = [
        SHARED / "training_selected_counts.csv.gz",
        SHARED / "rnaseq_train_sample_coordinates.csv",
        SHARED / "ranked_top250_features.csv",
        SHARED / "frozen_gene_order.csv",
    ]
    missing = [str(x) for x in paths if not x.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing validated Figure 3 outputs:\n" + "\n".join(missing)
        )
    counts = pd.read_csv(paths[0], index_col=0)
    samples = pd.read_csv(paths[1])
    panel = pd.read_csv(paths[2])
    learned = pd.read_csv(paths[3])
    sample_col = "sample_id" if "sample_id" in samples else "sample"
    gene_col = "gene_symbol" if "gene_symbol" in panel else panel.columns[0]
    genes = panel[gene_col].astype(str).tolist()
    counts.index = counts.index.astype(str)
    if set(genes) == set(counts.columns):
        expression = counts.loc[:, genes].to_numpy(float)
    elif set(genes) == set(counts.index):
        expression = counts.loc[genes].T.to_numpy(float)
    else:
        raise ValueError("Shared counts do not contain the fixed top-250 panel")
    if sample_col in samples and set(samples[sample_col].astype(str)) == set(
        counts.index
    ):
        counts = counts.loc[samples[sample_col].astype(str)]
        expression = counts.loc[:, genes].to_numpy(float)
    clr = clrcycle.clr_transform(expression)
    clr -= clr.mean(axis=0, keepdims=True)
    if clr.shape != (188, 250):
        raise ValueError(f"Expected 188x250 shared training matrix; got {clr.shape}")
    return clr, samples, pd.Series(genes, name="gene"), learned


def ordering_benchmark(clr, samples, genes, learned):
    """Compare clrcycle and ElPiGraph circular gene orderings."""
    # Histological peak stage supplies an external reference for both methods.
    stage = samples.histology_stage.to_numpy(int)
    learned = learned.set_index("feature")
    clrcycle_phase = learned.loc[genes, "clrcycle_angle"].to_numpy(float)
    clrcycle_rank = circular_ranks(clrcycle_phase)
    stage_means = np.vstack([clr[stage == x].mean(0) for x in range(1, 8)])
    peak_stage = 1 + np.argmax(stage_means, axis=0)
    reference = TWO_PI * (peak_stage - 1) / 7
    np.random.seed(17)
    # Fit ElPiGraph once to genes represented by their across-sample CLR profiles.
    graph = elpigraph.computeElasticPrincipalCircle(
        clr.T,
        NumNodes=30,
        Lambda=0.01,
        Mu=0.1,
        MaxNumberOfIterations=30,
        Do_PCA=True,
        ReduceDimension=np.arange(20),
        CenterData=True,
        n_cores=1,
        verbose=False,
    )[0]
    elpi_phase, projected = project_to_cycle(
        clr.T, graph["NodePositions"], graph["Edges"][0]
    )
    elpi_rank = circular_ranks(elpi_phase)
    between, agreement = align_phase(elpi_rank, clrcycle_rank)
    residual = np.angle(np.exp(1j * (between - clrcycle_rank)))
    rows = []
    errors = {}
    for method, rank in (("clrcycle", clrcycle_rank), ("ElPiGraph", elpi_rank)):
        aligned, R = align_phase(rank, reference)
        error = 7 / TWO_PI * np.abs(np.angle(np.exp(1j * (aligned - reference))))
        errors[method] = error
        rows.append(
            {
                "method": method,
                "n_genes": 250,
                "R_vs_histology_peak_stage": R,
                "mean_absolute_circular_error_stages": error.mean(),
                "median_absolute_circular_error_stages": np.median(error),
            }
        )
    benchmark = pd.DataFrame(rows)
    supplementary = pd.DataFrame(
        [
            {
                "comparison": "menstrual_top250_between_method_gene_ordering",
                "n_genes": 250,
                "circular_rank_agreement_R": agreement,
                "circular_rank_MAE_radians": np.mean(np.abs(residual)),
            }
        ]
    )
    benchmark.to_csv(OUT / "ordering_metrics.csv", index=False)
    supplementary.to_csv(OUT / "ordering_rank_agreement.csv", index=False)
    detail = pd.DataFrame(
        {
            "gene": genes,
            "histology_peak_stage": peak_stage,
            "clrcycle_rank_angle": clrcycle_rank,
            "elpigraph_rank_angle": elpi_rank,
        }
    )
    detail.to_csv(OUT / "gene_order_benchmark.csv", index=False)
    return benchmark, supplementary, graph, projected, peak_stage, clrcycle_rank, errors


def fourier_design(day, k):
    """Construct an intercept plus ``k`` Fourier harmonic pairs."""
    angle = TWO_PI * (day - 1) / 28
    return np.column_stack(
        [
            np.ones(len(angle)),
            *[f(h * angle) for h in range(1, k + 1) for f in (np.cos, np.sin)],
        ]
    )


def linear_predict(day, y):
    """Return in-sample predictions from a linear trend."""
    return (
        np.column_stack([np.ones(len(day)), day])
        @ np.linalg.lstsq(np.column_stack([np.ones(len(day)), day]), y, rcond=None)[0]
    )


def fourier_predict(day, y, k):
    """Return in-sample predictions from a ``k``-harmonic Fourier model."""
    x = fourier_design(day, k)
    return x @ np.linalg.lstsq(x, y, rcond=None)[0]


def spline_predict(day, y):
    """Return cyclic spline predictions using ridge-selected regularization."""
    pred = np.empty_like(y)
    for gene in range(y.shape[1]):
        model = make_pipeline(
            SplineTransformer(
                degree=3, knots=KNOTS, extrapolation="periodic", include_bias=False
            ),
            StandardScaler(),
            RidgeCV(alphas=ALPHAS),
        )
        model.fit(day[:, None], y[:, gene])
        pred[:, gene] = model.predict(day[:, None])
    return pred


def r2(y, pred):
    """Compute the coefficient of determination for fitted values."""
    total = np.square(y)
    residual = np.square(y - pred)
    return float(1 - residual.sum() / total.sum()), 1 - residual.sum(0) / total.sum(0)


def distortion_benchmark(clr, samples, genes):
    """Measure how model complexity changes fitted expression trajectories."""
    # Compare increasingly flexible trajectory models on the same selected genes.
    day_map = {1: 2.5, 2: 6.0, 3: 9.5, 4: 13.5, 5: 17.5, 6: 21.5, 7: 26.0}
    day = samples.histology_stage.astype(int).map(day_map).to_numpy(float)
    recon = {
        "linear": linear_predict(day, clr),
        "k=1": fourier_predict(day, clr, 1),
        "k=1-2": fourier_predict(day, clr, 2),
        "k=1-3": fourier_predict(day, clr, 3),
        "cyclic spline": spline_predict(day, clr),
    }
    evaluation = {k: r2(clr, v) for k, v in recon.items()}
    k1 = evaluation["k=1"][0]
    summary = pd.DataFrame(
        {"model": list(recon), "in_sample_r2": [evaluation[x][0] for x in recon]}
    )
    summary["gain_vs_k1"] = summary.in_sample_r2 - k1
    metrics = pd.DataFrame(
        {
            "gene": genes,
            "r2_linear": evaluation["linear"][1],
            "r2_k1": evaluation["k=1"][1],
            "r2_k1_2": evaluation["k=1-2"][1],
            "r2_k1_3": evaluation["k=1-3"][1],
            "r2_cyclic_spline": evaluation["cyclic spline"][1],
        }
    )
    for label in ("linear", "k1_2", "k1_3", "cyclic_spline"):
        metrics[f"delta_r2_{label}"] = metrics[f"r2_{label}"] - metrics.r2_k1
    summary.to_csv(OUT / "harmonic_distortion_summary.csv", index=False)
    metrics.to_csv(OUT / "harmonic_distortion_cv.csv", index=False)
    return summary, metrics


def plot(clr, ordering, summary, metrics, peak, clrrank, graph, projected, errors):
    """Render all panels of Supplementary Figure 1 from saved metrics."""
    BASE_SIZE = 18
    TITLE_SIZE = 21
    PANEL_SIZE = 24
    LEGEND_SIZE = 15
    ANNOTATION_SIZE = 15
    METRIC_SIZE = 19
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": BASE_SIZE,
            "axes.titlesize": TITLE_SIZE,
            "axes.titleweight": "bold",
            "axes.labelsize": BASE_SIZE,
            "xtick.labelsize": BASE_SIZE,
            "ytick.labelsize": BASE_SIZE,
            "legend.fontsize": LEGEND_SIZE,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    colors = [
        "#7A4EAB",
        "#3B73B9",
        "#2A9D8F",
        "#59A14F",
        "#E1C542",
        "#F28E2B",
        "#D94F5C",
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(0.5, 8.5), 7)
    fig = plt.figure(figsize=(14, 14), constrained_layout=True)
    gs = fig.add_gridspec(
        2,
        3,
        width_ratios=(1, 1, 0.12),
        height_ratios=(1, 1.02),
        hspace=0.05,
        wspace=0.08,
    )
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
    ]
    clrcycle_r = float(
        ordering.set_index("method").loc["clrcycle", "R_vs_histology_peak_stage"]
    )
    elpi_r = float(
        ordering.set_index("method").loc["ElPiGraph", "R_vs_histology_peak_stage"]
    )
    axes[0].scatter(
        np.cos(clrrank),
        np.sin(clrrank),
        c=peak,
        cmap=cmap,
        norm=norm,
        s=24,
        linewidths=0,
    )
    axes[0].add_patch(plt.Circle((0, 0), 1, fill=False, color=".65", lw=0.9))
    axes[0].set_aspect("equal")
    axes[0].set(
        xlim=(-1.18, 1.18), ylim=(-1.18, 1.18), title="clrcycle learnt gene order"
    )
    axes[0].set_title("clrcycle learnt gene order", y=1.18, pad=0)
    axes[0].text(
        0.5,
        1.065,
        f"R={clrcycle_r:.3f}; mean phase error={errors['clrcycle'].mean():.2f} stages",
        transform=axes[0].transAxes,
        ha="center",
        va="bottom",
        fontsize=METRIC_SIZE,
    )
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    for spine in axes[0].spines.values():
        spine.set_visible(False)
    view = PCA(n_components=2, svd_solver="full").fit(clr.T)
    gene_xy = view.transform(clr.T)
    node_xy = view.transform(graph["NodePositions"])
    projected_xy = view.transform(projected)
    order = cycle_node_order(graph["Edges"][0], len(node_xy))
    closed = np.r_[order, order[0]]
    axes[1].scatter(
        gene_xy[:, 0],
        gene_xy[:, 1],
        color=".82",
        s=10,
        linewidths=0,
        label="genes in 2D view",
        zorder=1,
    )
    axes[1].add_collection(
        LineCollection(
            np.stack([gene_xy, projected_xy], axis=1),
            colors=".75",
            linewidths=0.35,
            alpha=0.45,
            zorder=1,
        )
    )
    axes[1].plot(
        node_xy[closed, 0],
        node_xy[closed, 1],
        color="black",
        lw=1.6,
        label="fitted principal circle",
        zorder=2,
    )
    axes[1].scatter(
        projected_xy[:, 0],
        projected_xy[:, 1],
        c=peak,
        cmap=cmap,
        norm=norm,
        s=25,
        linewidths=0,
        label="projected genes",
        zorder=3,
    )
    axes[1].set_aspect("equal", adjustable="datalim")
    axes[1].set(xlabel="PC1", ylabel="PC2")
    axes[1].set_title("ElPiGraph principal circle", y=1.18, pad=0)
    axes[1].text(
        0.5,
        1.065,
        f"R={elpi_r:.3f}; mean phase error={errors['ElPiGraph'].mean():.2f} stages",
        transform=axes[1].transAxes,
        ha="center",
        va="bottom",
        fontsize=METRIC_SIZE,
    )
    axes[1].legend(frameon=False, fontsize=LEGEND_SIZE, loc="best")
    positions = np.arange(len(summary))
    values = summary.in_sample_r2.to_numpy(float)
    gains = summary.gain_vs_k1.to_numpy(float)
    axes[2].plot(positions, values, marker="o", color="#3569B0", lw=2.2, markersize=7)
    for position, value, gain in zip(positions, values, gains):
        axes[2].annotate(
            f"{value:.3f}\n({gain:+.3f})",
            (position, value),
            xytext=(0, -16),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=ANNOTATION_SIZE,
        )
    axes[2].set_xticks(positions, ["linear", "k=1", "k=1–2", "k=1–3", "cyclic\nspline"])
    axes[2].set_ylim(0, 0.5)
    axes[2].set(
        xlabel="Model fitted to histological\nmenstrual stage",
        ylabel=r"panel-wide $R^2$",
        title="Panel-wide variance explained\nby temporal models",
    )
    columns = [
        "delta_r2_linear",
        "delta_r2_k1_2",
        "delta_r2_k1_3",
        "delta_r2_cyclic_spline",
    ]
    labels = ["linear", "k=1–2", "k=1–3", "cyclic spline"]
    distortion_colors = ["#6B6B6B", "#56B4E9", "#D55E00", "#7A5195"]
    combined = np.concatenate([metrics[column].to_numpy(float) for column in columns])
    finite = combined[np.isfinite(combined)]
    bins = np.linspace(finite.min(), finite.max(), 28)
    for column, label, color in zip(columns, labels, distortion_colors):
        delta = metrics[column].to_numpy(float)
        axes[3].hist(
            delta,
            bins=bins,
            alpha=0.42,
            color=color,
            label=f"{label} (median={np.nanmedian(delta):+.3f})",
        )
    axes[3].axvline(0, color="black", ls="--", lw=1.1)
    axes[3].set(
        xlabel=r"gene-wise $\Delta R^2$ relative to k=1",
        ylabel="genes",
        title="Gene-wise gain relative to\nthe first-harmonic model",
    )
    axes[3].legend(frameon=False, fontsize=LEGEND_SIZE)
    for letter, ax in zip("ABCD", axes):
        ax.text(
            -0.32 if letter == "A" else -0.16,
            1.20,
            letter,
            transform=ax.transAxes,
            fontsize=PANEL_SIZE,
            fontweight="bold",
            ha="left",
            va="bottom",
        )
    shared = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    colorbar = fig.colorbar(
        shared,
        ax=axes[:2],
        orientation="horizontal",
        fraction=0.065,
        pad=0.07,
        aspect=50,
    )
    colorbar.set_ticks(np.arange(1, 8))
    colorbar.set_ticklabels(
        [
            "1 Menstrual",
            "2 Early proliferative",
            "3 Mid proliferative",
            "4 Late proliferative",
            "5 Early secretory",
            "6 Mid secretory",
            "7 Late secretory",
        ],
        rotation=32,
        ha="right",
        rotation_mode="anchor",
    )
    colorbar.ax.tick_params(labelsize=14, pad=4)
    colorbar.set_label("Histology stage with maximum mean gene expression", labelpad=24)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(
            ROOT / "figures" / f"Supplementary_Figure1.{ext}",
            dpi=300 if ext == "png" else None,
        )
    plt.close(fig)


def main():
    """Run the benchmarks and optionally render the supplementary figure."""
    # Numerical tables are always generated; rendering remains opt-in.
    p = argparse.ArgumentParser()
    p.add_argument("--plot", action="store_true")
    args = p.parse_args()
    clr, samples, genes, learned = load_shared()
    ordering, rank, graph, projected, peak, clrrank, errors = ordering_benchmark(
        clr, samples, genes, learned
    )
    summary, metrics = distortion_benchmark(clr, samples, genes)
    packages = {
        name: version(name)
        for name in (
            "clrcycle",
            "numpy",
            "pandas",
            "matplotlib",
            "scikit-learn",
            "elpigraph-python",
        )
    }
    (OUT / "software_versions.json").write_text(
        json.dumps({"python": platform.python_version(), **packages}, indent=2) + "\n"
    )
    if args.plot:
        plot(clr, ordering, summary, metrics, peak, clrrank, graph, projected, errors)
    print(ordering.to_string(index=False))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
