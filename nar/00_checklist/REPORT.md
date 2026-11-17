# 00_checklist - report

**What this folder is:** the full list of everything Kong et al. 2020 (Nat Commun 11:5485, DOI 10.1038/s41467-020-19313-8)
reports: every main figure panel, supplementary figure panel, supplementary table, source-data sheet and text-only number.

- `PAPER_CHECKLIST.md`: for people. It has the key claims with their statistics, the target cancers for generalising,
  the checklist table, places where the paper does not match itself, and which items the original `python/` code covers.
- `checklist.tsv`: for machines. It has the same 77 rows, plus a `search_token` and a `status` column (`todo`).
- `tables/source_data_facts.tsv`: numbers I worked out from the paper's source data (tagged [derived] in the checklist).
- `code/build_checklist.py`: rebuilds both checklist files. `code/derive_source_data_facts.py` rebuilds the facts table.
  Both need no randomness and no network. They read from `nar/_cache/`.

**Counts:** 31 main-figure items, 37 supplementary-figure items, 4 supplementary tables, 3 source-data items and 2 text items, so 77 in total.

**How I got the data:** nature.com and media.springernature.com are blocked here. The main article PDF and figures came from
the PMC open-access S3 bucket (`pmc-oa-opendata.s3.amazonaws.com/PMC7599252.1/`). The figure legends and the text with
figure numbers came from NCBI efetch (JATS XML). The PubMed MCP full text strips figure numbers, so I used it only as a cross-check.
The supplementary PDF and the source-data xlsx came from static-content.springer.com.

**Data file versions (sha256):**
- 41467_2020_19313_MOESM1_ESM.pdf  baf1c7f4d603eb2c87aa6c461fef30626e55c8e4b0e8108b8c7fe97f282c7cf2 (2,645,264 bytes)
- 41467_2020_19313_MOESM2_ESM.pdf  bfc626963bae8f300c1d2bde3fadfa5829e873e3528be6286dc71023f646f8ce (2,199,157 bytes)
- 41467_2020_19313_MOESM3_ESM.xlsx c1e065aa9490694abb9183c51388c29881190180691bbc2d71506032e3cbc657 (23,499,402 bytes)
- 41467_2020_Article_19313.pdf     1d4d34cad1b88adc88c0c61928160a1b963e77deb69840f7c61695b229a247d3 (1,535,353 bytes)
- PMC7599252_jats.xml              f031c8ab4074db96f6096d4acaabf518d1c792fe39b01f072bf60790471e2394

**Limits:** Some values were read by eye from bar heights or line points. These are the bar heights in Fig 2b/2f and Supp Fig 2,
the curve values in Supp Fig 18, and the weights in Supp Fig 11. They are tagged "[read from plot, approx]". I did not copy
the explained-variance numbers printed on the Supp Fig 14 axes (16 PCA panels).

What I think: the most important thing to know before reproducing the paper is in section 5 of the checklist. The bladder
biomarker is rank 2, not rank 1. The colon biomarker only just passes the proximity cut-off (z = -1.2857). The text calls
"5FU-treated" patients n = 114, which almost certainly includes FOLFOX.
