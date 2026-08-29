#!/usr/bin/env python3
"""Download every public source file required to reproduce the paper figures."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import sys
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"

FILES = {
    "reference/Homo_sapiens.GRCh38.91.gtf.gz":
        "https://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz",
    "GSE54650/GSE54650_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54650/matrix/GSE54650_series_matrix.txt.gz",
    "GSE54650/GPL6246.annot.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6246/annot/GPL6246.annot.gz",
    "GSE121265/GSE121265_fucci-counts.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121265/suppl/GSE121265_fucci-counts.txt.gz",
    "GSE121265/GSE121265_fucci-annotation.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121265/suppl/GSE121265_fucci-annotation.txt.gz",
    "GSE121265/GSE121265_fucci-annotation-description.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121265/suppl/GSE121265_fucci-annotation-description.txt.gz",
    "GSE121265/GSE121265_eset-final.rds.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121265/suppl/GSE121265_eset-final.rds.gz",
    "GSE121265/fit.quant.rds":
        "https://media.githubusercontent.com/media/jhsiao999/peco-paper/master/data/fit.quant.rds",
    "GSE121265/phenotypes_singletonsamples.txt":
        "https://raw.githubusercontent.com/jdblischak/fucci-seq/master/output/data-overview.Rmd/phenotypes_singletonsamples.txt",
    "GSE234354/GSE234354_gene_count_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234354/suppl/GSE234354_gene_count_matrix.txt.gz",
    "GSE234354/GSE234354-GPL11154_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234354/matrix/GSE234354-GPL11154_series_matrix.txt.gz",
    "GSE234354/GSE234354-GPL20301_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234354/matrix/GSE234354-GPL20301_series_matrix.txt.gz",
    "GSE234354/teh_source_data.xlsx":
        "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-41979-z/MediaObjects/41467_2023_41979_MOESM8_ESM.xlsx",
    "GSE111976/GSE111976_README.txt":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111976/suppl/GSE111976_README.txt",
    "GSE111976/GSE111976_ct.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111976/suppl/GSE111976_ct.csv.gz",
    "GSE111976/GSE111976_ct_endo_10x.rds.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111976/suppl/GSE111976_ct_endo_10x.rds.gz",
    "GSE111976/GSE111976_summary_10x_day_donor_ctype.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111976/suppl/GSE111976_summary_10x_day_donor_ctype.csv.gz",
    "GSE111976/GSE111976_summary_10x_donor_phase.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111976/suppl/GSE111976_summary_10x_donor_phase.csv.gz",
    "GSE111976/GSE111976_summary_C1_day_donor_ctype.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111976/suppl/GSE111976_summary_C1_day_donor_ctype.csv.gz",
    "GSE111976/GSE111976_summary_C1_donor_phase.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111976/suppl/GSE111976_summary_C1_donor_phase.csv.gz",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path, force: bool = False) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not force:
        print(f"exists  {destination.relative_to(ROOT)}")
        return
    partial = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "SEP-2026-reproducibility/1.0"})
    print(f"fetch   {url}")
    try:
        with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as output:
            shutil.copyfileobj(response, output, length=1024 * 1024)
        partial.replace(destination)
    finally:
        if partial.exists():
            partial.unlink()


def prepare_inputs(force: bool = False) -> None:
    """Create the small number of local files needed from downloaded archives."""
    source = RAW / "GSE111976" / "GSE111976_ct_endo_10x.rds.gz"
    destination = PROCESSED / "GSE111976_ct_endo_10x.rds"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not force:
        print(f"exists  {destination.relative_to(ROOT)}")
        return
    print(f"prepare {destination.relative_to(ROOT)}")
    with gzip.open(source, "rb") as compressed, destination.open("wb") as output:
        shutil.copyfileobj(compressed, output, length=1024 * 1024)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="replace existing downloads")
    args = parser.parse_args()
    for relative, url in FILES.items():
        download(url, RAW / relative, force=args.force)
    prepare_inputs(force=args.force)

    manifest = {
        relative: {
            "url": url,
            "bytes": (RAW / relative).stat().st_size,
            "sha256": sha256(RAW / relative),
        }
        for relative, url in FILES.items()
    }
    manifest_path = RAW / "download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote   {manifest_path.relative_to(ROOT)}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, urllib.error.URLError) as error:
        print(f"download failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
