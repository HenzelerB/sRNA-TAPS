# -*- coding: utf-8 -*-
"""
srnataps.utils
Shared utilities used across all sRNA-TAPS modules.
"""

import logging
import os
import re
import sys
from pathlib import Path
from typing import Optional


CANONICAL_CONDITIONS = ("no_treat", "pb_ctrl", "treat")

CONDITION_ALIASES = {
    "no_treat": {
        "no_treat", "no-treat", "notreat", "no treat", "untreated",
        "untr", "ctrl_untreated", "control_untreated",
    },
    "pb_ctrl": {
        "pb_ctrl", "pb-ctrl", "pbctrl", "pb control", "pb_control",
        "pb only", "pb_only", "pyridine_borane", "pyridine borane",
        "pb", "ctrl_pb",
    },
    "treat": {
        "treat", "treated", "tet_pb", "tet+pb", "tet-pb", "tet pb",
        "taps", "full_taps", "tet", "tet_plus_pb",
    },
}


def _condition_key(value: str) -> str:
    """Return a lowercase comparison key for condition aliases."""
    text = str(value or "").strip().lower()
    text = re.sub(r"[\s\-.]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text


def normalize_condition(value: str, *, strict: bool = False) -> str:
    """
    Normalize user/sample condition labels to sRNA-TAPS canonical names.

    Canonical internal names are:
        no_treat  - no chemistry / untreated baseline
        pb_ctrl   - pyridine borane only
        treat     - TET oxidation plus pyridine borane

    Unknown labels are returned unchanged unless strict=True.
    """
    original = str(value or "").strip()
    key = _condition_key(original)
    alias_keys = {
        canonical: {_condition_key(alias) for alias in aliases}
        for canonical, aliases in CONDITION_ALIASES.items()
    }
    for canonical, aliases in alias_keys.items():
        if key in aliases:
            return canonical

    # Conservative fallback for common lab-specific labels.
    if re.search(r"^(no|un).*(treat|tr|chem|tx)", key):
        return "no_treat"
    if key.startswith("pb") and "treat" not in key and "tet" not in key:
        return "pb_ctrl"
    if "tet" in key or (("treat" in key or "treated" in key) and not key.startswith("no")):
        return "treat"

    if strict:
        allowed = ", ".join(CANONICAL_CONDITIONS)
        raise ValueError(f"Unknown condition label '{original}'. Expected one of: {allowed}")
    return original


def is_condition(value: str, canonical: str) -> bool:
    """Return True if value normalizes to the requested canonical condition."""
    return normalize_condition(value) == canonical


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
    import pysam

    bai = bam_path + ".bai"
    if not os.path.exists(bai):
        pysam.index(bam_path)


def ensure_fasta_indexed(fasta_path: str) -> None:
    """Index FASTA with samtools faidx if .fai does not exist."""
    import pysam

    fai = fasta_path + ".fai"
    if not os.path.exists(fai):
        pysam.faidx(fasta_path)


def count_mapped_reads(bam_path: str) -> int:
    """Return number of primary mapped reads in BAM."""
    import pysam

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        n   = bam.mapped
        bam.close()
        return max(n, 1)
    except Exception:
        return 1


def get_chroms(bam_path: str) -> list[str]:
    """Return list of chromosome names from BAM header."""
    import pysam

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return [sq["SN"] for sq in bam.header.to_dict()["SQ"]]


def detect_cell_line(sample_name: str) -> str:
    """Detect cell line from sample name string."""
    if "HEK"   in sample_name: return "HEK"
    if "Caco2" in sample_name: return "Caco2"
    return "unknown"


def detect_condition(sample_name: str) -> str:
    """Detect TAPS condition from sample name string."""
    sample = str(sample_name)
    if re.search(r"pb[-_]?ctrl|pb[-_]?control|pb[-_]?only", sample, re.IGNORECASE):
        return "pb_ctrl"
    if re.search(r"no[-_]?treat|untreated|untr", sample, re.IGNORECASE):
        return "no_treat"
    if re.search(r"(^|[_-])treat(ed)?([_-]|$)|tet[-_+]?pb|full[-_]?taps", sample, re.IGNORECASE):
        return "treat"
    return "unknown"
