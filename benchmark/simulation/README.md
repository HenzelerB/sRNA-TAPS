# Rebuilding the benchmark truth set

These scripts regenerate `benchmark/fig2c_truth_set.tsv` and the simulated reads used in the manuscript benchmark. `probes/` holds the vetted loci. Only an indexed Ensembl GRCh38 release 112 `dna.toplevel` FASTA is needed in addition.

```bash
cd benchmark/simulation
FASTA=/path/to/Homo_sapiens.GRCh38.dna.toplevel.fa

python build_truth.py --fasta "$FASTA" --out truth_shared.tsv
python build_truth_v3.py --shared truth_shared.tsv --fasta "$FASTA" --out fig2c_truth_set.tsv
cmp fig2c_truth_set.tsv ../fig2c_truth_set.tsv

python simulate_taps_srna_v3.py --probe-dir probes --fasta "$FASTA" \\
    --truth fig2c_truth_set.tsv --outdir fastq --reads 1000000
```

`cmp` prints nothing when the rebuilt table is identical to the published one. Requires `pysam`.
