"""Exact numerical analysis shared by the GSE121265 FUCCI panels."""

from __future__ import annotations
from pathlib import Path
import subprocess
import gzip
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.manifold import SpectralEmbedding
from umap import UMAP
import clrcycle

TWO_PI = 2 * np.pi
SEED = 17
N_SWAPS = 4000


def extract_source_tables(eset: Path, fit_quant: Path, out: Path):
    """Use R to extract FUCCI metadata, feature annotations, and ranking tables."""
    out.mkdir(parents=True, exist_ok=True)
    meta = out / "GSE121265_metadata.tsv"
    features = out / "GSE121265_features.tsv"
    ranking = out / "GSE121265_published_top100.tsv"
    inner = out / "GSE121265_eset-final.rds"
    with gzip.open(eset, "rb") as source, inner.open("wb") as target:
        target.write(source.read())
    code = r"""suppressPackageStartupMessages(library(Biobase)); a<-commandArgs(TRUE)
x<-readRDS(gzfile(a[1])); pd<-pData(x); pc<-prcomp(pd[,c("rfp.median.log10sum.adjust","gfp.median.log10sum.adjust")],center=TRUE,scale.=TRUE)
ang<-atan2(pc$x[,2],pc$x[,1]); ang[ang<0]<-ang[ang<0]+2*pi
m<-data.frame(cell_id=rownames(pd),experiment=pd$experiment,chip_id=pd$chip_id,filter_all=pd$filter_all,theta=pd$theta,theta_reconstructed=((2*pi)-ang)%%(2*pi),check.names=FALSE)
write.table(m,a[3],sep="\t",row.names=FALSE,quote=FALSE); fd<-fData(x)
f<-data.frame(feature_id=rownames(fd),name=fd$name,source=fd$source,check.names=FALSE); write.table(f,a[4],sep="\t",row.names=FALSE,quote=FALSE)
fit<-readRDS(a[2]); pve<-vapply(fit,function(g) g$trend.pve,numeric(1)); o<-order(pve,decreasing=TRUE,na.last=NA)[seq_len(100)]
r<-data.frame(rank=seq_len(100),ENSG=names(pve)[o],pve=pve[o]); write.table(r,a[5],sep="\t",row.names=FALSE,quote=FALSE)"""
    subprocess.run(
        [
            "Rscript",
            "-e",
            code,
            str(inner),
            str(fit_quant),
            str(meta),
            str(features),
            str(ranking),
        ],
        check=True,
    )
    return meta, features, ranking


def load_inputs(counts_path, meta_path, features_path, ranking_path):
    """Load, validate, and align FUCCI counts, cells, and the top-100 panel."""
    meta = pd.read_csv(meta_path, sep="\t")
    meta = meta.loc[meta.filter_all.astype(bool)].set_index("cell_id")
    if len(meta) != 888:
        raise ValueError(f"Expected 888 cells; found {len(meta)}")
    err = np.angle(np.exp(1j * (meta.theta - meta.theta_reconstructed)))
    if np.max(np.abs(err)) > 1e-12:
        raise ValueError("FUCCI theta reconstruction failed")
    ranking = pd.read_csv(ranking_path, sep="\t")
    features = pd.read_csv(features_path, sep="\t").set_index("feature_id")
    ranking["gene"] = ranking.ENSG.map(features.name)
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    expression = counts.loc[ranking.ENSG].T
    expression.index = ["-".join(str(x).split(".")[-2:]) for x in expression.index]
    return expression.loc[meta.index], ranking, meta


def align_phase(raw, truth):
    """Align inferred phase to measured FUCCI phase by reflection and rotation."""
    candidates = []
    for reflected in (False, True):
        base = (-raw if reflected else raw) % TWO_PI
        rotation = np.angle(np.mean(np.exp(1j * (truth - base))))
        aligned = (base + rotation) % TWO_PI
        error = np.angle(np.exp(1j * (aligned - truth)))
        candidates.append(
            (
                float(abs(np.mean(np.exp(1j * error)))),
                float(np.mean(abs(error))),
                reflected,
                float(rotation),
                aligned,
            )
        )
    return max(candidates, key=lambda x: (x[0], -x[1]))


def run_clrcycle(expression, ranking, meta, out):
    """Fit clrcycle, align phases, and write FUCCI numerical artifacts."""
    out.mkdir(parents=True, exist_ok=True)
    result = clrcycle.fit(expression, max_features=None, seed=SEED, n_swaps=N_SWAPS)
    coords = result.coordinates.set_index("sample").loc[expression.index.astype(str)]
    theta = meta.loc[expression.index, "theta"].to_numpy(float)
    R, mae, reflected, rotation, _ = align_phase(
        coords.clrcycle_angle.to_numpy(), theta
    )
    u = coords.clrcycle_cosine.to_numpy(float)
    v = coords.clrcycle_sine.to_numpy(float)
    v_oriented = -v if reflected else v
    au = u * np.cos(rotation) - v_oriented * np.sin(rotation)
    av = u * np.sin(rotation) + v_oriented * np.cos(rotation)
    aa = np.mod(np.arctan2(av, au), TWO_PI)
    cells = (
        coords.rename(
            columns={
                "clrcycle_cosine": "clrcycle_cosine_raw",
                "clrcycle_sine": "clrcycle_sine_raw",
                "clrcycle_angle": "clrcycle_angle_raw",
            }
        )
        .rename_axis("sample")
        .reset_index()
    )
    cells["clrcycle_cosine_aligned"] = au
    cells["clrcycle_sine_aligned"] = av
    cells["clrcycle_angle_aligned"] = aa
    cells["measured_fucci_theta"] = theta
    cells["phase_residual_radians"] = np.angle(np.exp(1j * (aa - theta)))
    extra = (
        meta.loc[expression.index, ["experiment", "chip_id"]]
        .reset_index()
        .rename(columns={"cell_id": "sample"})
    )
    cells = cells.merge(extra, on="sample", validate="one_to_one")
    order = result.feature_order.merge(
        ranking, left_on="feature", right_on="ENSG", validate="one_to_one"
    )
    ga = order.clrcycle_angle.to_numpy()
    ga = -ga if reflected else ga
    order["clrcycle_angle_aligned"] = (ga + rotation) % TWO_PI
    cells.to_csv(out / "clrcycle_sample_coordinates.csv", index=False)
    order.to_csv(out / "clrcycle_gene_order.csv", index=False)
    ranking.to_csv(out / "published_gene_panel.csv", index=False)
    clr = clrcycle.clr_transform(expression.to_numpy(float))
    centered = clr - clr.mean(axis=0, keepdims=True)
    cov = centered.T @ centered / (len(centered) - 1)
    eig = np.linalg.eigvalsh(cov)
    pd.DataFrame(cov, index=expression.columns, columns=expression.columns).to_csv(
        out / "clr_gene_covariance.csv"
    )
    metrics = {
        "R": R,
        "mean_absolute_circular_error_radians": mae,
        "rho_clrcycle": result.rho,
        "objective_F": result.objective,
        "lambda_1": eig[-1],
        "lambda_2": eig[-2],
        "alignment_reflected": reflected,
        "alignment_rotation_radians": rotation,
        "n_cells": len(cells),
        "n_genes": len(order),
        "seed": SEED,
        "n_swaps": N_SWAPS,
    }
    pd.DataFrame([metrics]).to_csv(out / "clrcycle_metrics.csv", index=False)
    return {
        "cells": cells,
        "order": order,
        "ranking": ranking,
        "centered_clr": centered,
        "theta": theta,
        "metrics": metrics,
    }


def fitted_circle_phase(xy):
    """Map a two-dimensional embedding to phase around its fitted ellipse."""
    d = np.column_stack((2 * xy[:, 0], 2 * xy[:, 1], np.ones(len(xy))))
    cx, cy, off = np.linalg.lstsq(d, np.square(xy).sum(1), rcond=None)[0]
    radius = np.sqrt(max(off + cx * cx + cy * cy, np.finfo(float).eps))
    fit = least_squares(
        lambda p: np.linalg.norm(xy - p[:2], axis=1) - np.exp(p[2]),
        [cx, cy, np.log(radius)],
        method="lm",
    )
    z = xy - fit.x[:2]
    return np.mod(np.arctan2(z[:, 1], z[:, 0]), TWO_PI)


def distance_stress(reference, xy):
    """Measure normalized distortion of pairwise distances in an embedding."""
    represented = pdist(xy)
    scale = np.dot(reference, represented) / np.dot(represented, represented)
    residual = reference - scale * represented
    return float(np.sqrt(np.dot(residual, residual) / np.dot(reference, reference)))


def run_embedding_benchmark(analysis, out):
    """Fit comparison embeddings and write their phase/distance metrics."""
    clr = analysis["centered_clr"]
    cells = analysis["cells"]
    theta = analysis["theta"]
    embeddings = {
        "CLRcycle first Fourier": cells[
            ["clrcycle_cosine_raw", "clrcycle_sine_raw"]
        ].to_numpy(),
        "PCA (2 PCs)": PCA(n_components=2, svd_solver="full").fit_transform(clr),
        "Spectral embedding (2D)": SpectralEmbedding(
            n_components=2,
            affinity="nearest_neighbors",
            n_neighbors=15,
            random_state=SEED,
            eigen_solver="arpack",
        ).fit_transform(clr),
        "UMAP (2D)": UMAP(
            n_components=2,
            n_neighbors=15,
            min_dist=0.1,
            metric="euclidean",
            random_state=SEED,
            n_jobs=1,
        ).fit_transform(clr),
    }
    reference = pdist(clr)
    rows = []
    coordinates = []
    for name, xy in embeddings.items():
        phase = (
            np.mod(np.arctan2(xy[:, 1], xy[:, 0]), TWO_PI)
            if name.startswith("CLRcycle")
            else fitted_circle_phase(xy)
        )
        R, *_ = align_phase(phase, theta)
        rows.append(
            {
                "method": name,
                "circular_agreement": R,
                "distance_stress": distance_stress(reference, xy),
            }
        )
        coordinates.append(
            pd.DataFrame(
                {
                    "sample": cells["sample"],
                    "method": name,
                    "x": xy[:, 0],
                    "y": xy[:, 1],
                    "measured_fucci_theta": theta,
                }
            )
        )
    metrics = pd.DataFrame(rows)
    table = pd.concat(coordinates, ignore_index=True)
    metrics.to_csv(out / "embedding_metrics.csv", index=False)
    table.to_csv(out / "embedding_coordinates.csv", index=False)
    return embeddings, metrics
