# Figure 2C benchmark ground truth

`fig2c_truth_set.tsv` is the synthetic ground-truth table used to benchmark
sRNA-TAPS against Bismark, asTAIR, and raSTAIR (manuscript Figure 2C). Each
row is one cytosine position planted at a known modification rate in the
simulated small-RNA validation libraries.

**738 planted sites** in two batches, distinguished by the `truth_class`
column:

- **`discovery`** (n=100) — the planted signal exists only in the treat
  (TET+PB) condition. At the same position in the PB-only control, the
  simulated rate is flat background noise (1–3%), indistinguishable from an
  unmodified site.
- **`shared`** (n=638) — the planted signal exists in *both* conditions: the
  PB-only rate is deliberately scaled to 5–15% of the site's own treat rate,
  rather than flat background. This is the harder, more realistic case,
  since a caller must correctly subtract a genuine chemistry-associated
  background rather than simply distinguish signal from noise.

## Columns

| Column | Description |
|---|---|
| `chrom` | Genomic chromosome of the planted site |
| `genomic_pos` | 1-based genomic coordinate of the planted cytosine |
| `strand` | Genomic strand (+/-) |
| `biotype` | Small-RNA biotype of the locus (miRNA, tRNA, rRNA, snRNA, snoRNA) |
| `gene_id` | Ensembl gene identifier of the locus |
| `planted_rate` | Ground-truth modification rate simulated at this site |
| `truth_class` | `discovery` or `shared` — see above |

## Usage

This file follows the `truth.tsv` coordinate schema expected by
`srnataps evaluate` (see the main README's Test Dataset section):

```bash
srnataps evaluate \
    --truth benchmark/fig2c_truth_set.tsv \
    --calls-dir 07d.replicate_calls \
    --samples-tsv samples.tsv \
    --outdir 10.truth_evaluation \
    --condition treat
```

Truth-based evaluation is for synthetic validation only. It is not used by
the caller and must not be used to tune thresholds on biological
experiments.
