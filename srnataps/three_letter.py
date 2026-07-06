# -*- coding: utf-8 -*-
"""Three-letter preprocessing and alignment restoration for RNA TAPS reads."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def convert_reference(input_path, c2t_path, g2a_path):
    """Write C-to-T and G-to-A reference FASTAs without changing coordinates."""
    Path(c2t_path).parent.mkdir(parents=True, exist_ok=True)
    Path(g2a_path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(input_path) as source, open(c2t_path, "w") as c2t, open(
        g2a_path, "w"
    ) as g2a:
        for line in source:
            if line.startswith(">"):
                c2t.write(line)
                g2a.write(line)
                continue
            sequence = line.rstrip("\r\n").upper()
            c2t.write(sequence.translate(str.maketrans({"C": "T"})) + "\n")
            g2a.write(sequence.translate(str.maketrans({"G": "A"})) + "\n")


def convert_fastq(input_path, output_path):
    """Convert read C bases to T while preserving names and qualities."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open_text(input_path) as source, open_text(output_path, "wt") as output:
        while True:
            name = source.readline()
            if not name:
                break
            sequence = source.readline()
            separator = source.readline()
            quality = source.readline()
            if not sequence or not separator or not quality:
                raise ValueError(f"incomplete FASTQ record after read {count}")
            output.write(name)
            output.write(sequence.rstrip("\r\n").upper().replace("C", "T") + "\n")
            output.write(separator)
            output.write(quality)
            count += 1
    return count


def load_original_reads(fastq_path):
    """Return read name to original sequence and quality."""
    originals = {}
    with open_text(fastq_path) as source:
        while True:
            name = source.readline()
            if not name:
                break
            sequence = source.readline().rstrip("\r\n")
            separator = source.readline()
            quality = source.readline().rstrip("\r\n")
            if not separator:
                raise ValueError(f"incomplete FASTQ record for {name.rstrip()}")
            key = name[1:].strip().split()[0]
            if key in originals:
                raise ValueError(f"duplicate FASTQ read name: {key}")
            originals[key] = (sequence, quality)
    return originals


def alignment_score(read):
    """Lower scores are preferred when selecting between conversion branches."""
    try:
        mismatches = int(read.get_tag("NM"))
    except KeyError:
        try:
            mismatches = int(read.get_tag("XM"))
        except KeyError:
            mismatches = 10**9
    return (
        mismatches,
        -int(read.mapping_quality),
        int(read.reference_id),
        int(read.reference_start),
        int(read.is_reverse),
    )


def restore_best_alignments(original_fastq, plus_bam, minus_bam, output_bam):
    """Select the best transformed alignment and restore original read bases."""
    import pysam

    originals = load_original_reads(original_fastq)
    best = {}

    plus = pysam.AlignmentFile(plus_bam, "rb")
    minus = pysam.AlignmentFile(minus_bam, "rb")
    try:
        for branch, bam in (("C2T", plus), ("G2A", minus)):
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped:
                    continue
                score = alignment_score(read)
                current = best.get(read.query_name)
                if current is None or score < current[0]:
                    best[read.query_name] = (score, branch, read)

        Path(output_bam).parent.mkdir(parents=True, exist_ok=True)
        with pysam.AlignmentFile(output_bam, "wb", header=plus.header) as output:
            for name, (_, branch, read) in best.items():
                if name not in originals:
                    raise ValueError(f"alignment read absent from FASTQ: {name}")
                sequence, quality = originals[name]
                read.query_sequence = sequence
                read.query_qualities = pysam.qualitystring_to_array(quality)
                for tag in ("MD", "NM", "XM"):
                    if read.has_tag(tag):
                        read.set_tag(tag, None)
                read.set_tag("XC", branch, value_type="Z")
                output.write(read)
    finally:
        plus.close()
        minus.close()
    return len(best), len(originals)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    reference = subparsers.add_parser("convert-reference")
    reference.add_argument("--input", required=True)
    reference.add_argument("--c2t", required=True)
    reference.add_argument("--g2a", required=True)

    fastq = subparsers.add_parser("convert-fastq")
    fastq.add_argument("--input", required=True)
    fastq.add_argument("--output", required=True)

    restore = subparsers.add_parser("restore")
    restore.add_argument("--original-fastq", required=True)
    restore.add_argument("--plus-bam", required=True)
    restore.add_argument("--minus-bam", required=True)
    restore.add_argument("--output-bam", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.command == "convert-reference":
        convert_reference(args.input, args.c2t, args.g2a)
    elif args.command == "convert-fastq":
        count = convert_fastq(args.input, args.output)
        print(f"Converted reads: {count:,}")
    elif args.command == "restore":
        mapped, total = restore_best_alignments(
            args.original_fastq,
            args.plus_bam,
            args.minus_bam,
            args.output_bam,
        )
        print(f"Restored alignments: {mapped:,}/{total:,}")


if __name__ == "__main__":
    main()
