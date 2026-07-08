"""Supported Ensembl species and reference-resolution helpers."""

from copy import deepcopy
from pathlib import Path


SPECIES = {
    "human": {
        "scientific_name": "Homo sapiens", "ensembl_name": "homo_sapiens",
        "file_prefix": "Homo_sapiens", "assembly": "GRCh38",
        "index_stem": "hg38", "aliases": ("homo_sapiens", "hsapiens", "hg38"),
    },
    "mouse": {
        "scientific_name": "Mus musculus", "ensembl_name": "mus_musculus",
        "file_prefix": "Mus_musculus", "assembly": "GRCm39",
        "index_stem": "mm39", "aliases": ("mus_musculus", "mmusculus", "mm39"),
    },
    "rat": {
        "scientific_name": "Rattus norvegicus", "ensembl_name": "rattus_norvegicus",
        "file_prefix": "Rattus_norvegicus", "assembly": "mRatBN7.2",
        "index_stem": "mratbn7.2", "aliases": ("rattus_norvegicus", "rnorvegicus"),
    },
    "zebrafish": {
        "scientific_name": "Danio rerio", "ensembl_name": "danio_rerio",
        "file_prefix": "Danio_rerio", "assembly": "GRCz11",
        "index_stem": "grcz11", "aliases": ("danio_rerio", "drerio"),
    },
    "fruit_fly": {
        "scientific_name": "Drosophila melanogaster",
        "ensembl_name": "drosophila_melanogaster",
        "file_prefix": "Drosophila_melanogaster", "assembly": "BDGP6.46",
        "index_stem": "bdgp6.46",
        "aliases": ("fruit-fly", "fly", "drosophila", "drosophila_melanogaster"),
    },
    "c_elegans": {
        "scientific_name": "Caenorhabditis elegans",
        "ensembl_name": "caenorhabditis_elegans",
        "file_prefix": "Caenorhabditis_elegans", "assembly": "WBcel235",
        "index_stem": "wbcel235",
        "aliases": ("c-elegans", "worm", "caenorhabditis_elegans"),
    },
    "chicken": {
        "scientific_name": "Gallus gallus", "ensembl_name": "gallus_gallus",
        "file_prefix": "Gallus_gallus",
        "assembly": "bGalGal1.mat.broiler.GRCg7b", "index_stem": "grcg7b",
        "aliases": ("gallus_gallus", "ggallus"),
    },
}


def supported_species():
    return deepcopy(SPECIES)


def canonical_species(value):
    normalized = str(value or "human").strip().lower().replace(" ", "_")
    for key, record in SPECIES.items():
        if normalized == key or normalized in record["aliases"]:
            return key
    raise ValueError(
        f"Unsupported species '{value}'. Supported species: {', '.join(SPECIES)}"
    )


def reference_details(reference, genome_dir):
    """Resolve local reference paths and Ensembl URLs without mutating config."""
    reference = dict(reference or {})
    key = canonical_species(reference.get("species", "human"))
    record = SPECIES[key]
    release = int(reference.get("ensembl_release", 112))
    genome_dir = Path(genome_dir)
    base = f"{record['file_prefix']}.{record['assembly']}"

    def unset(value):
        return not value or str(value).startswith("<REQUIRED:")

    genome_fa = reference.get("genome_fa")
    gtf = reference.get("gtf")
    index = reference.get("bowtie1_index")
    if unset(genome_fa):
        genome_fa = str(genome_dir / f"{base}.dna.toplevel.fa")
    if unset(gtf):
        gtf = str(genome_dir / f"{base}.{release}.gtf")
    if unset(index):
        index = str(genome_dir / base)

    ftp = f"https://ftp.ensembl.org/pub/release-{release}"
    return {
        **record,
        "species": key,
        "ensembl_release": release,
        "genome_fa": genome_fa,
        "gtf": gtf,
        "bowtie1_index": index,
        "genome_url": (
            f"{ftp}/fasta/{record['ensembl_name']}/dna/"
            f"{base}.dna.toplevel.fa.gz"
        ),
        "gtf_url": f"{ftp}/gtf/{record['ensembl_name']}/{base}.{release}.gtf.gz",
    }


def resolve_reference_config(reference, genome_dir):
    """Populate derived reference paths in a workflow configuration."""
    details = reference_details(reference, genome_dir)
    for key in ("species", "ensembl_release", "genome_fa", "gtf", "bowtie1_index"):
        reference[key] = details[key]
    return details
