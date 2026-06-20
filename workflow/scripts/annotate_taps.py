#!/usr/bin/env python3
"""
annotate_taps.py — Annotate TAPS call TSVs with gene names from Ensembl GTF.

Adds three columns to each TAPS TSV:
    gene_name     e.g. "hsa-miR-21-5p", "SNHG14", "RNA45S5"
    gene_id       Ensembl gene ID, e.g. "ENSG00000284190"
    gene_biotype  GTF biotype, e.g. "miRNA", "rRNA", "lncRNA"

Chromosome positions are retained unchanged for IGV compatibility.
When no GTF feature overlaps a site, all three columns are set to ".".

Priority for multi-overlap: annotations whose gene_biotype matches the
biotype wildcard (i.e. the BAM the site came from) are preferred over
annotations from other biotypes at the same position.

Called by Snakemake via `script:` directive.
"""

import os
import re
import subprocess
import sys

import pandas as pd


# ── helpers ───────────────────────────────────────────────────────────────────

def _parse_attr(attr_str: str, key: str) -> str:
    """Extract a quoted value from a GTF attributes string."""
    m = re.search(rf'{key} "([^"]+)"', attr_str)
    return m.group(1) if m else "."


BIOTYPE_GTF_MAP = {
    "miRNA":  "miRNA",
    "tRNA":   "tRNA",
    "rRNA":   "rRNA",
    "snoRNA": "snoRNA",
    "snRNA":  "snRNA",
    "lncRNA": "lncRNA",
    "piRNA":  "piRNA",
    "other":  "",          # no preferred biotype for "other"
}


# ── main ──────────────────────────────────────────────────────────────────────

def annotate_taps(tsv_in: str, gtf: str, tsv_out: str,
                  biotype: str, log_fh=None):
    """
    Annotate one TAPS TSV with gene names from the GTF.

    Parameters
    ----------
    tsv_in   : path to input TAPS TSV (output of caller.py)
    gtf      : path to Ensembl GFF/GTF
    tsv_out  : path for annotated output TSV
    biotype  : Snakemake wildcard, e.g. "miRNA" — used to prioritise matching GTF entries
    log_fh   : file handle for log output (optional)
    """

    def log(msg):
        if log_fh:
            print(msg, file=log_fh, flush=True)

    # ── read input ────────────────────────────────────────────────────────────
    try:
        df = pd.read_csv(tsv_in, sep="\t", dtype={"chrom": str})
    except pd.errors.EmptyDataError:
        log(f"[SKIP] Empty file: {tsv_in}")
        # write empty annotated header
        pd.DataFrame(columns=["chrom", "start", "end", "strand", "context",
                               "mod_count", "unmod_count", "coverage",
                               "mod_rate", "pvalue", "padj", "snp_flag",
                               "gene_name", "gene_id", "gene_biotype"]
                     ).to_csv(tsv_out, sep="\t", index=False)
        return
    except Exception as exc:
        log(f"[ERROR] Could not read {tsv_in}: {exc}")
        raise

    if df.empty:
        df["gene_name"]    = pd.Series(dtype=str)
        df["gene_id"]      = pd.Series(dtype=str)
        df["gene_biotype"] = pd.Series(dtype=str)
        df.to_csv(tsv_out, sep="\t", index=False)
        return

    log(f"Annotating {len(df)} sites  biotype={biotype}")

    # ── temp file paths ───────────────────────────────────────────────────────
    bed_tmp      = tsv_out + ".tmp.bed"
    gtf_gene_tmp = tsv_out + ".tmp_genes.gtf"

    try:
        # ── write BED (row index as name for round-trip merge) ────────────────
        with open(bed_tmp, "w") as bf:
            for i, row in df.iterrows():
                bf.write(f"{row['chrom']}\t{int(row['start'])}\t{int(row['end'])}\t{i}\n")

        # ── filter GTF to gene features (much smaller, much faster) ───────────
        log("Filtering GTF to gene features...")
        with open(gtf_gene_tmp, "w") as gf:
            awk = subprocess.run(
                ["awk", '$3 == "gene"', gtf],
                stdout=gf, stderr=subprocess.PIPE, check=True
            )
        gene_count = int(subprocess.run(
            ["wc", "-l", gtf_gene_tmp],
            capture_output=True, text=True
        ).stdout.split()[0])
        log(f"  {gene_count} gene features extracted")

        # ── bedtools intersect: left outer join ───────────────────────────────
        log("Running bedtools intersect...")
        result = subprocess.run(
            ["bedtools", "intersect",
             "-a", bed_tmp,
             "-b", gtf_gene_tmp,
             "-wa", "-wb", "-loj"],
            capture_output=True, text=True, check=True
        )

        # ── parse intersect output ────────────────────────────────────────────
        #  BED: 4 cols (chrom start end idx)
        #  GTF: 9 cols (chrom src feature start end score strand frame attrs)
        #  Total: 13 cols → attrs at index 12

        target_biotype = BIOTYPE_GTF_MAP.get(biotype, "")

        # dict: idx → (gene_name, gene_id, gene_biotype, biotype_matched)
        gene_info: dict = {}

        for line in result.stdout.splitlines():
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue

            try:
                idx = int(parts[3])
            except ValueError:
                continue

            # No intersection (-loj fills with ".")
            if len(parts) < 13 or parts[4] == ".":
                gene_info.setdefault(idx, (".", ".", ".", False))
                continue

            attrs = parts[12] if len(parts) > 12 else ""
            gname = _parse_attr(attrs, "gene_name")
            gid   = _parse_attr(attrs, "gene_id")
            gbio  = _parse_attr(attrs, "gene_biotype")
            is_match = bool(target_biotype) and (gbio == target_biotype)

            if idx not in gene_info:
                gene_info[idx] = (gname, gid, gbio, is_match)
            elif is_match and not gene_info[idx][3]:
                # Upgrade to biotype-matching annotation
                gene_info[idx] = (gname, gid, gbio, True)

        # fill any index not seen (e.g. edge case)
        for i in df.index:
            gene_info.setdefault(i, (".", ".", ".", False))

        df["gene_name"]    = [gene_info[i][0] for i in df.index]
        df["gene_id"]      = [gene_info[i][1] for i in df.index]
        df["gene_biotype"] = [gene_info[i][2] for i in df.index]

        n_annotated = (df["gene_name"] != ".").sum()
        log(f"  {n_annotated}/{len(df)} sites annotated with a gene name")

    finally:
        for tmp in (bed_tmp, gtf_gene_tmp):
            try:
                os.unlink(tmp)
            except OSError:
                pass

    df.to_csv(tsv_out, sep="\t", index=False)
    log(f"Written: {tsv_out}")


# ── Snakemake entry point ─────────────────────────────────────────────────────

if __name__ == "__main__":
    with open(str(snakemake.log[0]), "w") as _log:
        annotate_taps(
            tsv_in  = str(snakemake.input.tsv),
            gtf     = str(snakemake.input.gtf),
            tsv_out = str(snakemake.output.tsv),
            biotype = str(snakemake.wildcards.biotype),
            log_fh  = _log,
        )
