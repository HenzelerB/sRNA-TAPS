# TAPSy knowledge layer

TAPSy retrieves reference material at question time; this does not fine-tune model weights.

## Coverage

- Text from all nine public HTML pages, with scripts, styles, navigation and SVG markup excluded.
- All exported atlas records, retaining biotype, identifiers, genomic/mature positions, modification percentage, coverage, status and structured citations.
- Thirteen curated research notes in `literature.json`: miRNA evidence, tRNA writers/fragmentation, mitochondrial oxidation pathways, rRNA and assay chemistry. This is a starting collection, not a systematic or exhaustive review. The notes distinguish directly checked publication summaries from interpretations of the project's curated master.
- Atlas totals computed from exported records rather than hand-written prose.

## Retrieval

The browser loads `index.json` only when a visitor sends a question. It selects up to four relevant current-page chunks, up to eight additional matching chunks, and at most eight exact-identifier site matches. A selected atlas detail record is prioritized. Site matching preserves miRNA arms and prioritizes requested mature positions. Generic biotype questions do not generate arbitrary gene matches. Limited site results explicitly report the matching and returned counts.

The system context instructs TAPSy to cite supplied sources, distinguish assay and evidence resolutions, acknowledge missing evidence, and avoid inferring global statistics from record subsets. These instructions guide model behavior but do not guarantee accuracy. Each answer includes an expandable list of reference material supplied to the model; this is not a claim that every listed source was cited or used in the answer.

## Refresh and validation

From the repository root:

```
python3 scripts/build_tapsy_knowledge.py
gjs tests/test_tapsy_knowledge.js
```

Rebuild after changing site content, atlas exports, or literature notes. Commit the generated `docs/knowledge/index.json` with the source edits. The index contains source-page hashes and only public site data plus the curated summaries. No local unpublished analysis directory is automatically ingested.

## Chat backend discovery and verification

The existing frontend calls `https://srna-taps-chat.bennett-henzeler.workers.dev` using its existing model identifier and API request format. Repository/workspace searches found the frontend endpoint reference but no Worker implementation or Wrangler configuration. The connected GitHub account exposes Huhu, sRNA-TAPS and bioconda-recipes; GitHub code search did not find a Wrangler configuration or another occurrence of the endpoint. These searches cannot establish whether the Worker was created in the Cloudflare dashboard or lives in an inaccessible repository.

No Cloudflare management connection is currently available in this session. A minimal test request on 2026-09-21 returned HTTP 403 from this environment. It did not verify model responses or whether the Worker forwards the `system` field. The knowledge payload and local retrieval tests are complete, but live behavior requires testing through the website and, if necessary, inspection of Cloudflare Workers & Pages → srna-taps-chat. Keep provider credentials server-side; none are needed in these static knowledge files.
