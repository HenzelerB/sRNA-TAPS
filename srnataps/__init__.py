"""
sRNA-TAPS: TAPS-based m5C detection pipeline for small RNA sequencing.

Detects 5-methylcytosine (m5C) and 5-hydroxymethylcytosine (5hmC) in small RNA
using TET-assisted pyridine borane sequencing (TAPS). Supports miRNA, tRNA,
rRNA, snoRNA, snRNA, piRNA, and lncRNA biotypes from human samples.

Chemistry:
    TAPS converts m5C and 5hmC → T via TET oxidation + pyridine borane reduction.
    Unmodified C stays as C. C→T in reads = modification signal (opposite of bisulfite).

Usage:
    srnataps init --outdir my_project
    srnataps run  --configfile my_project/config.yaml
    srnataps module fastqc --configfile my_project/config.yaml
"""

__version__ = "0.2.4"
__author__  = "Bennett Henzeler"
__email__   = "bennett.henzeler@cup.lmu.de"
