# -*- coding: utf-8 -*-
"""
srnataps.utils
Shared utilities used across all sRNA-TAPS modules.
"""

import logging
import os
import sys
from pathlib import Path
from typing import Optional

import pysam


def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Return a consistently formatted logger."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(
            "[%(asctime)s] %(levelname)s %(name)s: %(message)s",
            datefmt="%H:%M:%S",
        ))
        logger.addHandler(handler)
    logger.setLevel(level)
    return logger


def ensure_bam_indexed(bam_path: str) -> None:
    """Index BAM if .bai does not exist."""
    bai = bam_path + ".bai"
    if not os.path.exists(bai):
        pysam.index(bam_path)


def ensure_fasta_indexed(fasta_path: str) -> None:
    """Index FASTA with samtools faidx if .fai does not exist."""
    fai = fasta_path + ".fai"
    if not os.path.exists(fai):
        pysam.faidx(fasta_path)


def count_mapped_reads(bam_path: str) -> int:
    """Return number of primary mapped reads in BAM."""
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        n   = bam.mapped
        bam.close()
        return max(n, 1)
    except Exception:
        return 1


def get_chroms(bam_path: str) -> list[str]:
    """Return list of chromosome names from BAM header."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return [sq["SN"] for sq in bam.header.to_dict()["SQ"]]


def detect_cell_line(sample_name: str) -> str:
    """Detect cell line from sample name string."""
    if "HEK"   in sample_name: return "HEK"
    if "Caco2" in sample_name: return "Caco2"
    return "unknown"


def detect_condition(sample_name: str) -> str:
    """Detect TAPS condition from sample name string."""
    if "pb_Ctrl"  in sample_name: return "pb_ctrl"
    if "no-treat" in sample_name: return "no_treat"
    if "treat"    in sample_name: return "treat"
    return "unknown"
