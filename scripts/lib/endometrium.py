#!/usr/bin/env python3
"""RNA-seq-first clrcycle discovery and frozen Quake donor transfer.

Expression alone determines feature selection and clrcycle geometry. Training
histology is loaded only after fitting, for the permitted reflection/rotation
and stage calibration. External annotations are loaded only after projection.
"""

from __future__ import annotations

import gzip, json, re, shutil, sys, zipfile
from xml.etree import ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import train_test_split

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
Q10 = RAW / "GSE111976"
DATA = RAW / "GSE234354"
C1 = RAW / "GSE111976"
GTF = RAW / "reference" / "Homo_sapiens.GRCh38.91.gtf.gz"
OUT = ROOT / "results" / "fig3" / "shared"
DIRS = {
    x: OUT / x
    for x in (
        "discovery_rnaseq",
        "transfer_10x",
        "transfer_c1",
        "figures",
        "tables",
        "metadata",
    )
}
sys.path.insert(0, str(ROOT / "vendor/clrcycle"))
import clrcycle
from lib.common import paper_ellipse, run_r

SEED, TOP_N, TEST_SIZE, N_SWAPS = 17, 250, 0.20, 4000
TWO_PI = 2 * np.pi
STAGES = np.arange(1, 8)
STAGE_DAY = {1: 2.5, 2: 6.0, 3: 9.5, 4: 13.5, 5: 17.5, 6: 21.5, 7: 26.0}
STAGE_COLORS = [
    "#7A4EAB",
    "#3B73B9",
    "#2A9D8F",
    "#59A14F",
    "#E1C542",
    "#F28E2B",
    "#D94F5C",
]
STAGE_NAMES = [
    "Menstrual",
    "Early prolif.",
    "Mid prolif.",
    "Late prolif.",
    "Early secr.",
    "Mid secr.",
    "Late secr.",
]
PHASE_TO_STAGE = {
    "menstrual": 1,
    "proliferative_early": 2,
    "proliferative": 3,
    "proliferative_late": 4,
    "secretory_early": 5,
    "secretory_early-mid": 6,
    "secretory_mid": 6,
    "secretory_late": 7,
}


def sheet_table(name: str, columns: list[str]) -> pd.DataFrame:
    """Extract a named, header-delimited table from the source workbook."""
    raw = _read_xlsx_sheet(DATA / "teh_source_data.xlsx", name)
    row = raw.index[raw.iloc[:, 0].eq(columns[0])][0]
    out = raw.iloc[row + 1 :, : len(columns)].copy()
    out.columns = columns
    return out.dropna(how="all")


def _read_xlsx_sheet(path: Path, sheet_name: str) -> pd.DataFrame:
    """Read the small source-data workbook without an optional Excel package."""
    ns = {
        "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
        "p": "http://schemas.openxmlformats.org/package/2006/relationships",
    }
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = [
                "".join(t.text or "" for t in si.findall(".//m:t", ns))
                for si in root.findall("m:si", ns)
            ]
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relation_id = next(
            s.attrib[f"{{{ns['r']}}}id"]
            for s in workbook.findall(".//m:sheet", ns)
            if s.attrib["name"] == sheet_name
        )
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target = next(
            r.attrib["Target"]
            for r in rels.findall("p:Relationship", ns)
            if r.attrib["Id"] == relation_id
        )
        sheet_path = (
            "xl/" + target.lstrip("/") if not target.startswith("xl/") else target
        )
        sheet = ET.fromstring(archive.read(sheet_path))
        rows = []
        for row in sheet.findall(".//m:sheetData/m:row", ns):
            values = {}
            for cell in row.findall("m:c", ns):
                ref = cell.attrib["r"]
                col = 0
                for char in re.match(r"[A-Z]+", ref).group(0):
                    col = col * 26 + ord(char) - 64
                value = cell.find("m:v", ns)
                if value is None:
                    parsed = None
                elif cell.attrib.get("t") == "s":
                    parsed = shared[int(value.text)]
                else:
                    try:
                        parsed = float(value.text)
                    except ValueError:
                        parsed = value.text
                values[col - 1] = parsed
            if values:
                rows.append([values.get(i) for i in range(max(values) + 1)])
        width = max(map(len, rows))
        return pd.DataFrame([r + [None] * (width - len(r)) for r in rows])


def gtf_map() -> dict[str, str]:
    """Build the Ensembl-to-symbol map used to aggregate RNA-seq counts."""
    ans = {}
    pattern = re.compile(r'gene_id "([^"]+)".*gene_name "([^"]+)"')
    with gzip.open(GTF, "rt") as handle:
        for line in handle:
            if "\tgene\t" not in line:
                continue
            hit = pattern.search(line.split("\t", 8)[8])
            if hit:
                ans[hit.group(1).split(".")[0]] = hit.group(2)
    assert len(ans) > 50000
    return ans


def rnaseq_expression() -> pd.DataFrame:
    """Load GSE234354 counts and aggregate Ensembl rows to gene symbols."""
    x = pd.read_csv(DATA / "GSE234354_gene_count_matrix.txt.gz", sep="\t", index_col=0)
    x.index = x.index.astype(str).str.split(".").str[0]
    names = gtf_map()
    x = x.loc[x.index.isin(names)].copy()
    x["gene_symbol"] = [names[g] for g in x.index]
    x = x.groupby("gene_symbol", sort=True).sum().T.astype(float)
    x.index.name = "sample_id"
    return x.loc[:, x.sum(axis=0) > 0]


def rnaseq_annotations() -> pd.DataFrame:
    """Load the authors' histology stages and independent molecular time."""
    # Figure 2e contains the authors' reclassified seven-stage pathology;
    # Figure 2g supplies the independent 0--100 molecular-model coordinate.
    stage = sheet_table("Figure 2e", ["sample_id", "histology_stage", "model_stage"])
    time = sheet_table(
        "Figure 2g", ["sample_id", "model_cycle_time", "coarse_pathology_stage"]
    )
    x = stage[["sample_id", "histology_stage"]].merge(
        time[["sample_id", "model_cycle_time"]], on="sample_id", validate="one_to_one"
    )
    x["model_cycle_time"] = pd.to_numeric(x.model_cycle_time, errors="coerce")
    x["histology_stage"] = pd.to_numeric(x.histology_stage, errors="coerce")
    return x.dropna().drop_duplicates("sample_id").set_index("sample_id")


def r_pseudobulk(platform: str, genes: list[str]) -> pd.DataFrame:
    """Aggregate single-cell expression to donor pseudobulks for one platform."""
    temp = DIRS[f"transfer_{platform}"]
    req, out = temp / "requested_genes.txt", temp / "pseudobulk.tsv"
    pd.Series(genes).to_csv(req, index=False, header=False)
    if platform == "10x":
        matrix, md = (
            PROCESSED / "GSE111976_ct_endo_10x.rds",
            Q10 / "GSE111976_summary_10x_day_donor_ctype.csv.gz",
        )
        code = r"""suppressPackageStartupMessages(library(Matrix));a<-commandArgs(TRUE);x<-readRDS(a[1]);
        rownames(x)<-x@Dimnames[[1]];colnames(x)<-x@Dimnames[[2]];md<-read.csv(gzfile(a[2]),check.names=FALSE);
        req<-readLines(a[3]);ok<-req%in%rownames(x);stopifnot(sum(ok)>=3);
        md<-md[match(colnames(x),md$cell_name),];stopifnot(all(md$cell_name==colnames(x)),!anyNA(md$donor));
        z<-sparse.model.matrix(~0+factor(md$donor));pb<-x[req[ok],,drop=FALSE]%*%z;
        colnames(pb)<-sub("factor(md$donor)","",colnames(pb),fixed=TRUE);write.table(as.matrix(t(pb)),a[4],sep="\t",quote=FALSE,col.names=NA)"""
        run_r(code, matrix, md, req, out)
        pb = pd.read_csv(out, sep="\t", index_col=0)
        pb.index = pb.index.astype(str)
        return pb
    c1 = C1
    md = pd.read_csv(c1 / "GSE111976_summary_C1_day_donor_ctype.csv.gz")
    donor = md.assign(donor=md.donor.astype(str)).set_index("cell_name").donor
    pieces = []
    for chunk in pd.read_csv(c1 / "GSE111976_ct.csv.gz", index_col=0, chunksize=2000):
        hit = chunk.loc[chunk.index.intersection(genes)]
        if len(hit):
            pieces.append(hit.T.groupby(donor.loc[hit.columns], sort=False).sum().T)
    pb = pd.concat(pieces).T
    pb.to_csv(out, sep="\t", index_label="donor")
    return pb


def tenx_gene_names() -> set[str]:
    """Return gene names measurable in the 10x source matrix."""
    out = DIRS["metadata"] / "tenx_gene_names.txt"
    code = r"""suppressPackageStartupMessages(library(Matrix));a<-commandArgs(TRUE);
    x<-readRDS(a[1]);write.table(data.frame(gene=x@Dimnames[[1]]),a[2],
    row.names=FALSE,col.names=FALSE,quote=FALSE)"""
    run_r(code, PROCESSED / "GSE111976_ct_endo_10x.rds", out)
    genes = set(out.read_text().splitlines())
    assert len(genes) > 10000
    return genes


def c1_gene_names() -> set[str]:
    """Return gene names measurable in the C1 source matrix."""
    path = C1 / "GSE111976_ct.csv.gz"
    return set(pd.read_csv(path, usecols=[0]).iloc[:, 0].dropna().astype(str))


def external_meta(platform: str) -> pd.DataFrame:
    """Load donor cycle-day, cell-count, and phase annotations."""
    root = Q10 if platform == "10x" else C1
    cells = pd.read_csv(root / f"GSE111976_summary_{platform}_day_donor_ctype.csv.gz")
    phase = pd.read_csv(root / f"GSE111976_summary_{platform}_donor_phase.csv.gz")
    cells["donor"] = cells.donor.astype(str)
    phase["donor"] = phase.donor.astype(str)
    d = cells.groupby("donor").agg(
        reported_cycle_day=("day", "first"), n_cells=("cell_name", "size")
    )
    return d.join(phase.set_index("donor")[["phase_canonical"]])


def orient(raw: np.ndarray, stage: np.ndarray) -> dict:
    """Choose the post-fit reflection and rotation that align training stages."""
    known = TWO_PI * (np.array([STAGE_DAY[int(s)] for s in stage]) - 1) / 28
    choices = []
    for reflected in (False, True):
        base = (-raw if reflected else raw) % TWO_PI
        rotation = float(np.angle(np.mean(np.exp(1j * (known - base)))))
        aligned = (base + rotation) % TWO_PI
        r = float(abs(np.mean(np.exp(1j * (aligned - known)))))
        choices.append((r, reflected, rotation, aligned))
    r, reflected, rotation, aligned = max(choices, key=lambda q: q[0])
    return {
        "training_stage_R": r,
        "reflected": reflected,
        "rotation_radians": rotation,
        "aligned": aligned,
    }


def apply_orientation(angle: np.ndarray, p: dict) -> np.ndarray:
    """Apply a previously learned reflection and rotation to phase angles."""
    return ((-angle if p["reflected"] else angle) + p["rotation_radians"]) % TWO_PI


def centers(angle, stage):
    """Compute the circular center of each histological stage."""
    return {
        int(s): float(np.mod(np.angle(np.mean(np.exp(1j * angle[stage == s]))), TWO_PI))
        for s in STAGES
    }


def predict(angle, ctr):
    """Assign each angle to its nearest circular stage center."""
    c = np.array([ctr[int(s)] for s in STAGES])
    d = abs(np.angle(np.exp(1j * (angle[:, None] - c))))
    return STAGES[np.argmin(d, axis=1)]


def project(x: pd.DataFrame, train: pd.DataFrame, order: pd.DataFrame) -> pd.DataFrame:
    """Project new samples through the frozen training feature geometry."""
    assert (
        list(x.columns) == list(train.columns) == sorted(train.columns)
    ), "frozen feature mismatch"
    a = clrcycle.clr_transform(train.to_numpy())
    b = clrcycle.clr_transform(x.to_numpy())
    b -= a.mean(axis=0, keepdims=True)
    pos = {g: i for i, g in enumerate(train.columns)}
    idx = np.array([pos[g] for g in order.feature])
    n = len(idx)
    theta = TWO_PI * np.arange(n) / n
    u = b[:, idx] @ (np.sqrt(2 / n) * np.cos(theta))
    v = b[:, idx] @ (np.sqrt(2 / n) * np.sin(theta))
    return pd.DataFrame(
        {
            "raw_cosine": u,
            "raw_sine": v,
            "raw_angle": np.mod(np.arctan2(v, u), TWO_PI),
            "radius": np.hypot(u, v),
        },
        index=x.index,
    )


def add_oriented(z, p):
    """Add aligned angles, coordinates, and 0--100 cycle time to a table."""
    z = z.copy()
    z["clrcycle_angle"] = apply_orientation(z.raw_angle.to_numpy(), p)
    z["clrcycle_cosine"] = z.radius * np.cos(z.clrcycle_angle)
    z["clrcycle_sine"] = z.radius * np.sin(z.clrcycle_angle)
    z["clrcycle_time_0_100"] = z.clrcycle_angle * 100 / TWO_PI
    return z


def training_ellipse_table(training, orientation, stage_centers):
    """Describe only the training ellipse displayed in Figure 4 A and C."""
    model = clrcycle.circular_agreement(training.raw_angle, training.model_cycle_time, period=100)
    coordinates = training[["raw_cosine", "raw_sine", "raw_angle"]].rename(columns={
        "raw_cosine": "clrcycle_cosine", "raw_sine": "clrcycle_sine", "raw_angle": "clrcycle_angle"})
    result = clrcycle.ClrCycleResult(coordinates, pd.DataFrame(), np.nan, np.nan)
    references = [
        ("stage", training.histology_stage.map(STAGE_DAY), 28, 1,
         (orientation["rotation_radians"], orientation["reflected"]), "cycle day calibrated to training histology"),
        ("model_time", training.model_cycle_time, 100, 0,
         (model["rotation"], model["reflection"]), "molecular-model time (0-100)"),
    ]
    rows = []
    for reference, metadata, period, origin, alignment, units in references:
        row = paper_ellipse(result, metadata, period=period, units=units, origin=origin, alignment=alignment)
        row.update(cohort="rnaseq_train", metadata_reference=reference, alignment_source="RNA-seq training calibration")
        for name in ("major_vertex_1", "major_vertex_2", "minor_vertex_1", "minor_vertex_2"):
            phase = apply_orientation(np.array([row[name+"_angle"]]), orientation)
            stage = int(predict(phase, stage_centers)[0])
            row[name+"_stage"] = stage
            row[name+"_stage_label"] = STAGE_NAMES[stage-1]
        rows.append(row)
    return pd.DataFrame(rows)


def score(z, stage_col=None, time_col=None, ctr=None, period=100.0) -> dict:
    """Calculate circular-time and stage-agreement validation metrics."""
    ans = {"n": int(len(z))}
    if time_col and z[time_col].notna().any():
        d = z.dropna(subset=[time_col])
        known = TWO_PI * d[time_col].to_numpy(float) / period
        residual = np.angle(np.exp(1j * (d.clrcycle_angle.to_numpy() - known)))
        rho, p = spearmanr(d[time_col], d.clrcycle_angle * period / TWO_PI)
        ans.update(
            reference_n=len(d),
            circular_phase_agreement_R=float(abs(np.mean(np.exp(1j * residual)))),
            mean_absolute_circular_error=float(
                np.mean(abs(residual)) * period / TWO_PI
            ),
            spearman_rho=float(rho),
            spearman_p=float(p),
        )
    if stage_col and ctr is not None and z[stage_col].notna().any():
        d = z.dropna(subset=[stage_col])
        true = d[stage_col].to_numpy(int)
        pred = predict(d.clrcycle_angle.to_numpy(), ctr)
        delta = np.abs(true - pred)
        ans.update(
            stage_n=len(d),
            stage_exact=float(np.mean(delta == 0)),
            stage_within_one_circular=float(np.mean(np.minimum(delta, 7 - delta) <= 1)),
        )
    return ans


def phase_compatible(phase, stage) -> bool:
    """Return whether a reported phase is compatible with a predicted stage."""
    allowed = {
        "menstrual": {1},
        "proliferative_early": {2},
        "proliferative_late": {3, 4},
        "proliferative": {2, 3, 4},
        "secretory_early": {5},
        "secretory_early-mid": {5, 6},
        "secretory_mid": {6},
        "secretory_late": {7},
    }
    return int(stage) in allowed.get(str(phase), set())


def stage_boundaries(ctr):
    """Compute circular midpoints between adjacent learned stage centers."""
    order = np.argsort(list(ctr.values()))
    stages = np.array(list(ctr))[order]
    c = np.array(list(ctr.values()))[order]
    b = np.mod(c + np.mod(np.roll(c, -1) - c, TWO_PI) / 2, TWO_PI)
    return pd.DataFrame({"stage_before": stages, "boundary_radians": b})


def gene_circle_annotations(train, train_meta, feature_order, orientation, selected):
    """Annotate the frozen gene order after fitting, using training histology."""
    hist_phase = TWO_PI * (
        np.array([STAGE_DAY[int(stage)] for stage in train_meta.histology_stage]) - 1
    ) / 28
    design = np.column_stack([np.ones(len(hist_phase)), np.cos(hist_phase), np.sin(hist_phase)])
    beta = np.linalg.lstsq(
        design, np.log1p(train[feature_order.feature]).to_numpy(), rcond=None
    )[0]
    peak = np.mod(np.arctan2(beta[2], beta[1]), TWO_PI)
    stage_angle = TWO_PI * (np.array([STAGE_DAY[int(stage)] for stage in STAGES]) - 1) / 28
    peak_stage = STAGES[np.argmin(abs(np.angle(np.exp(1j * (peak[:, None] - stage_angle)))), axis=1)]
    amplitude = np.hypot(beta[1], beta[2])
    table = pd.DataFrame({
        "gene": feature_order.feature,
        "gene_angle": apply_orientation(feature_order.clrcycle_angle.to_numpy(), orientation),
        "posthoc_peak_histology_stage": peak_stage,
        "posthoc_harmonic_amplitude": amplitude,
        "posthoc_amplitude_rank": pd.Series(amplitude).rank(ascending=False, method="first").astype(int),
        "highlighted_top20_training_variance": feature_order.feature.isin(selected.index[:20]),
    })
    return table


def run_shared_analysis():
    """Run discovery, frozen transfers, validation, and artifact export."""
    # Stage 1: load RNA-seq discovery data and make a reproducible held-out split.
    for d in DIRS.values():
        d.mkdir(parents=True, exist_ok=True)
    expr = rnaseq_expression()
    ids = expr.index.intersection(rnaseq_annotations().index)
    # Split uses histology only for sampling balance; labels do not enter selection or fit.
    split_meta = rnaseq_annotations().loc[ids]
    train_ids, test_ids = train_test_split(
        ids.to_numpy(),
        test_size=TEST_SIZE,
        random_state=SEED,
        stratify=split_meta.histology_stage.astype(int).to_numpy(),
    )
    train_all = expr.loc[train_ids]
    test_all = expr.loc[test_ids]

    # Stage 2: define a cross-platform measurable universe, then select genes
    # using only the RNA-seq training samples.
    # Establish measurability first, then rank the shared universe using only
    # RNA-seq training expression. External sample structure and labels never
    # enter either operation.
    tenx_genes = tenx_gene_names()
    c1_genes = c1_gene_names()
    shared_universe = sorted(set(expr.columns) & tenx_genes & c1_genes)
    assert len(shared_universe) >= TOP_N
    variance = (
        np.log1p(train_all[shared_universe]).var(ddof=1).sort_values(ascending=False)
    )
    selected = variance.head(TOP_N)
    common = sorted(selected.index)
    selected.rename("training_log1p_variance").to_csv(
        DIRS["discovery_rnaseq"] / "ranked_top250_features.csv"
    )
    pd.DataFrame(
        {
            "gene": selected.index,
            "training_rank": np.arange(1, TOP_N + 1),
            "training_log1p_variance": selected.to_numpy(),
        }
    ).to_csv(DIRS["tables"] / "final_common_genes.csv", index=False)
    pd.DataFrame(
        {
            "platform": ["RNA-seq", "10x", "C1", "three-way shared"],
            "n_genes": [
                len(expr.columns),
                len(tenx_genes),
                len(c1_genes),
                len(shared_universe),
            ],
        }
    ).to_csv(DIRS["metadata"] / "feature_universe_counts.csv", index=False)
    tenx = r_pseudobulk("10x", common)
    c1 = r_pseudobulk("c1", common)
    train = train_all[common]
    test = test_all[common]

    # Stage 3: learn the clrcycle geometry once from RNA-seq training data.
    fit = clrcycle.fit(train, max_features=None, seed=SEED, n_swaps=N_SWAPS)
    assert set(fit.feature_order.feature) == set(common)

    # Only now attach training labels for orientation and calibration.
    train_meta = rnaseq_annotations().loc[train.index]
    test_meta = rnaseq_annotations().loc[test.index]
    raw_train = (
        fit.coordinates.set_index("sample")
        .loc[train.index]
        .rename(
            columns={
                "clrcycle_cosine": "raw_cosine",
                "clrcycle_sine": "raw_sine",
                "clrcycle_angle": "raw_angle",
                "clrcycle_radius": "radius",
            }
        )
    )
    p = orient(raw_train.raw_angle.to_numpy(), train_meta.histology_stage.to_numpy(int))
    aligned = p.pop("aligned")
    ctr = centers(aligned, train_meta.histology_stage.to_numpy(int))
    bounds = stage_boundaries(ctr)
    train_z = add_oriented(raw_train, p).join(train_meta)
    test_z = add_oriented(project(test, train, fit.feature_order), p).join(test_meta)
    train_z["predicted_stage"] = predict(train_z.clrcycle_angle.to_numpy(), ctr)
    test_z["predicted_stage"] = predict(test_z.clrcycle_angle.to_numpy(), ctr)

    # Stage 4: transfer the frozen feature order and orientation to both
    # single-cell platforms; neither transfer dataset is used for refitting.
    # External expression projection precedes the corresponding annotation joins.
    tenx_raw = project(tenx[common], train, fit.feature_order)
    c1_raw = project(c1[common], train, fit.feature_order)
    tenx_z = add_oriented(tenx_raw, p)
    c1_z = add_oriented(c1_raw, p)
    tenx_z = tenx_z.join(external_meta("10x"))
    c1_z = c1_z.join(external_meta("C1"))
    for z in (tenx_z, c1_z):
        z["predicted_stage"] = predict(z.clrcycle_angle.to_numpy(), ctr)
    for z in (tenx_z, c1_z):
        z["reported_day_angle"] = TWO_PI * (z.reported_cycle_day - 1) / 28
        z["reported_day_stage"] = (
            np.searchsorted(
                [4, 7, 11, 15, 19, 23, 28], z.reported_cycle_day, side="left"
            )
            + 1
        )
        z["canonical_stage"] = z.phase_canonical.map(PHASE_TO_STAGE)
        z["canonical_phase_compatible"] = [
            phase_compatible(a, b) for a, b in zip(z.phase_canonical, z.predicted_stage)
        ]

    # Stage 5: quantify held-out and cross-platform agreement.
    metrics = {
        "RNAseq held-out": score(
            test_z, "histology_stage", "model_cycle_time", ctr, 100
        ),
        "10x transfer": score(tenx_z, "canonical_stage", "reported_cycle_day", ctr, 28),
        "C1 transfer": score(c1_z, "canonical_stage", "reported_cycle_day", ctr, 28),
    }
    metrics["10x transfer"]["canonical_phase_compatible_fraction"] = float(
        tenx_z.canonical_phase_compatible.mean()
    )
    metrics["C1 transfer"]["canonical_phase_compatible_fraction"] = float(
        c1_z.canonical_phase_compatible.mean()
    )

    # Stage 6: write the complete numerical audit trail used by Figure 4.
    training_ellipse_table(train_z, p, ctr).to_csv(
        DIRS["tables"] / "training_ellipse_geometry.csv", index=False)
    for platform, z in (
        ("rnaseq_train", train_z),
        ("rnaseq_heldout", test_z),
        ("10x", tenx_z),
        ("c1", c1_z),
    ):
        z.assign(dataset=platform).to_csv(
            DIRS["tables"] / f"{platform}_sample_coordinates.csv",
            index_label="sample_id",
        )
    train_z.to_csv(
        DIRS["discovery_rnaseq"] / "training_sample_coordinates.csv",
        index_label="sample_id",
    )
    test_z.to_csv(
        DIRS["discovery_rnaseq"] / "heldout_sample_coordinates.csv",
        index_label="sample_id",
    )
    tenx_z.to_csv(DIRS["transfer_10x"] / "sample_coordinates.csv", index_label="donor")
    c1_z.to_csv(DIRS["transfer_c1"] / "sample_coordinates.csv", index_label="donor")
    fit.feature_order.to_csv(DIRS["tables"] / "frozen_gene_order.csv", index=False)
    weights = fit.feature_weights.copy()
    weights["histology_stage"] = weights["sample"].map(train_meta.histology_stage)
    weights["model_cycle_time"] = weights["sample"].map(train_meta.model_cycle_time)
    weights["training_log1p_variance"] = weights["feature"].map(selected)
    weights.to_csv(DIRS["tables"] / "feature_weights.csv.gz", index=False)
    fit.feature_order[["feature", "clrcycle_angle"]].rename(
        columns={"feature": "gene", "clrcycle_angle": "gene_angle"}
    ).to_csv(DIRS["tables"] / "gene_angles.csv", index=False)
    pd.DataFrame({"stage": list(ctr), "center_radians": list(ctr.values())}).to_csv(
        DIRS["tables"] / "stage_centers.csv", index=False
    )
    bounds.to_csv(DIRS["tables"] / "stage_boundaries.csv", index=False)
    pd.DataFrame.from_dict(metrics, orient="index").rename_axis("dataset").to_csv(
        DIRS["tables"] / "summary_metrics.csv"
    )
    orientation = {**p, "training_stage_R": p["training_stage_R"]}
    gene_circle_annotations(train, train_meta, fit.feature_order, orientation, selected).to_csv(
        DIRS["tables"] / "gene_circle_posthoc_histology_stage.csv", index=False
    )
    (DIRS["tables"] / "training_orientation.json").write_text(
        json.dumps(orientation, indent=2) + "\n"
    )
    train.to_csv(OUT / "training_selected_counts.csv.gz", compression="gzip")
    test.to_csv(OUT / "heldout_selected_counts.csv.gz", compression="gzip")

    # Record every fixed setting and headline metric in one machine-readable file.
    summary = {
        "settings": {
            "seed": SEED,
            "top_n": TOP_N,
            "test_size": TEST_SIZE,
            "n_swaps": N_SWAPS,
            "preprocessing": {
                "rnaseq": "raw gene counts; duplicate symbols summed; CLR pseudocount internal to clrcycle",
                "single_cell": "donor sums across all cells; CLR",
            },
        },
        "n_rnaseq_train": len(train),
        "n_rnaseq_heldout": len(test),
        "n_rnaseq_10x_c1_shared_genes": len(shared_universe),
        "n_final_common_genes": len(common),
        "training_rho": fit.rho,
        "orientation": orientation,
        "metrics": metrics,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"RNA-seq train/held-out: {len(train)}/{len(test)}\nretained genes: {len(common)}\ntraining rho: {fit.rho:.6f}"
    )
    for k, v in metrics.items():
        print(k, json.dumps(v, sort_keys=True))
    public = ROOT / "results" / "fig3"
    public.mkdir(parents=True, exist_ok=True)
    exports = {
        DIRS["tables"] / "training_ellipse_geometry.csv": "training_ellipse_geometry.csv",
        DIRS["discovery_rnaseq"]
        / "ranked_top250_features.csv": "ranked_top250_features.csv",
        DIRS["tables"] / "frozen_gene_order.csv": "frozen_gene_order.csv",
        DIRS["tables"] / "feature_weights.csv.gz": "endometrium_feature_weights.csv.gz",
        DIRS["tables"]
        / "gene_circle_posthoc_histology_stage.csv": "gene_circle_posthoc_histology_stage.csv",
        DIRS["tables"]
        / "rnaseq_train_sample_coordinates.csv": "rnaseq_train_sample_coordinates.csv",
        DIRS["tables"]
        / "rnaseq_heldout_sample_coordinates.csv": "rnaseq_heldout_sample_coordinates.csv",
        DIRS["tables"] / "10x_sample_coordinates.csv": "10x_sample_coordinates.csv",
        DIRS["tables"] / "c1_sample_coordinates.csv": "c1_sample_coordinates.csv",
        DIRS["tables"] / "stage_centers.csv": "stage_centers.csv",
        DIRS["tables"] / "stage_boundaries.csv": "stage_boundaries.csv",
        DIRS["tables"] / "training_orientation.json": "training_orientation.json",
        DIRS["tables"] / "summary_metrics.csv": "summary_metrics.csv",
        OUT / "training_selected_counts.csv.gz": "training_selected_counts.csv.gz",
        OUT / "heldout_selected_counts.csv.gz": "heldout_selected_counts.csv.gz",
        OUT / "summary.json": "summary.json",
    }
    for source, name in exports.items():
        shutil.copyfile(source, public / name)
    return summary


if __name__ == "__main__":
    run_shared_analysis()
