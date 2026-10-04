# Genetic-interaction networks (network priors)

Built from GitHub-hosted copies (NCBI/EBI/zenodo/figshare blocked). Each `{name}.tsv.gz` has columns
`gene_a, gene_b, type, score, source` (+ extras). Load with `obd.trial_cohorts.load_gi_network(name)`.

| name | pairs | genes | types | description |
|---|---|---|---|---|
| synlethdb_human_sl | 35943 | 9856 | SL:35943 | SynLethDB 2.0 human SL pairs (all evidence types: CRISPR/RNAi screens, computational prediction, text mining, ...). score = SynLethDB statistic_score (0-1). |
| slkb_crispr_pairs | 16059 | 2568 | SL:16059 | SLKB combinatorial-CRISPR SL calls (as mirrored in bhklab/synthLethal SLKB_pairs.csv). score = SL_score (more negative = stronger SL). |
| isle_clinical_sl_fdr0.1 | 2326 | 2153 | SL:2326 | ISLE genome-wide clinically relevant SL network, FDR 0.1 (Lee et al. Nat Commun 2018; 2,326 pairs). Edges parsed from the Cytoscape session (.cys) XGMML; no per-edge scores in the session. |
| isle_clinical_sl_fdr0.2 | 21534 | 8510 | SL:21534 | ISLE genome-wide clinically relevant SL network, FDR 0.2 (Lee et al. 2018). Edges parsed from the Cytoscape session (.cys) XGMML; no per-edge scores in the session. |
| isle_drug_target_csl | 232 | 237 | SL:232 | ISLE drug-target clinical SL (cSL) network used for TCGA drug-response prediction (gene_b = drug target). Edges parsed from the Cytoscape session (.cys) XGMML; no per-edge scores in the session. |
| isle_gold_standard_sl | 6033 | 2944 | SL:6033 | ISLE experimentally reported SL gold-standard set (literature screens: shRNA/sgRNA/drug/mutation screens), sl.golden.set.RData. source = screen(s). |
| bhklab_pancancer_sl | 24243 | 11358 | SL:24243 | bhklab/synthLethal ISLE-style pan-cancer SL inference (TCGA co-inactivation depletion + CRISPR + survival + phylogeny); score = q_value. Columns SR_DD_* are the synthetic-rescue (down-down) test. |

## Sources (sha256 of the downloaded source file)

- **synlethdb_human_sl**: https://raw.githubusercontent.com/bhklab/synthLethal/main/data/Human_SL.csv (sha256 43f7f66bfad325639ab25284419c8ed792ec00de36f6596e633db1c9baa774d0)
- **slkb_crispr_pairs**: https://raw.githubusercontent.com/bhklab/synthLethal/main/data/SLKB_pairs.csv (sha256 9cd60624f66661dcd297e36fe99b50879f2db31b87866352e3facdb4d5fc5c01)
- **isle_clinical_sl_fdr0.1**: https://github.com/jooslee/ISLE/blob/master/networks/ISLE_clinical_SL_network_FDR_0.1.cys (sha256 321afb73a2f7b47d89f48c3236e72e4c1aa3fd7ace82a97d414c540f63fdc39b)
- **isle_clinical_sl_fdr0.2**: https://github.com/jooslee/ISLE/blob/master/networks/ISLE_clinical_SL_network_FDR_0.2.cys (sha256 fac3df58aecd937ffbd050a71a042578f9f1b8df93bf916d74cf904a4f5b7a7d)
- **isle_drug_target_csl**: https://github.com/jooslee/ISLE/blob/master/networks/ISLE_drug_cSL_network.cys (sha256 176633025ddc3c888fdc5287d4e2737fa74a732cdfa3f69888df589a59cc104d)
- **isle_gold_standard_sl**: https://github.com/jooslee/ISLE/blob/master/data/sl.golden.set.RData (sha256 afbee9a2b876c6fadf18b75fc8e215dd7049bbb58c7b58abe716c3abff28f0e2)
- **bhklab_pancancer_sl**: https://raw.githubusercontent.com/bhklab/synthLethal/main/data/SL_pairs_pan-cancer_binarize_expression_cox.csv (sha256 0c41e9c45f4c5d10731634aa0912f9bd8df28049a1181b22d7e4663682da58a9)

## sha256 of the deliverables

- bhklab_pancancer_sl.tsv.gz: `d52623385d38f0a17e3660a6217aadfb95acfa88ea141a183153042344e779a3`
- isle_clinical_sl_fdr0.1.tsv.gz: `dba622fafd713da4781b2f00b1dfb1f7e091f38e8634c70e5a6b238cf364926b`
- isle_clinical_sl_fdr0.2.tsv.gz: `95cf065347c67a3f00b10b6b69553348e7a6a02f4db01f2dcda8e55e471c58ac`
- isle_drug_target_csl.tsv.gz: `80bae2e23f97bf35adf1ab149288f332226571c065dc987ed08c27bb38a50b44`
- isle_gold_standard_sl.tsv.gz: `3ea9d5bd629489105102038afa6483490e657531b26eea91265aa0a6a1b486b6`
- slkb_crispr_pairs.tsv.gz: `d35b44b87064f9535c1dc4573d19d03461fb5690b333bc05427a7e5b666895af`
- synlethdb_human_sl.tsv.gz: `b79e2cb72edbedfe915ae6b428cd4251195282a147e373a8d3a989ee56f135ad`

## Notes
- ISLE = Lee JS, et al. Harnessing synthetic lethality to predict the response to cancer treatment. Nat Commun 2018;9:2546 (github jooslee/ISLE). FDR 0.1 / 0.2 clinical SL networks and the drug-target cSL network were parsed from the Cytoscape `.cys` sessions (XGMML edge labels `A (interacts with) B`); the sessions carry no per-edge statistics, so `score` is NA. Direction is not meaningful (SL is symmetric), except in `isle_drug_target_csl` where gene_b is the drug target.
- `isle_gold_standard_sl`: the experimental SL gold standard used by ISLE for benchmarking (sl.golden.set.RData, all flag = 1); `source` lists the screen(s)/dataset(s) supporting the pair.
- `synlethdb_human_sl`: SynLethDB 2.0 human SL table as mirrored in bhklab/synthLethal (`Human_SL.csv`, Neo4j export). Mixed evidence: GenomeRNAi, CRISPR/CRISPRi, computational predictions, text mining. Filter on `source` for experimental-only priors.
- `slkb_crispr_pairs`: SLKB (combinatorial CRISPR screens, Gu et al. NAR 2023) SL calls as mirrored in bhklab/synthLethal (`SLKB_pairs.csv`; only SL-positive rows were present).
- `bhklab_pancancer_sl`: bhklab/synthLethal ISLE-style re-implementation (pan-cancer, binarised expression, Cox survival); `score` = q_value; `SR_DD_*` columns are the synthetic-rescue (down-down) screen p/q values.
- bhklab/synthLethal `CSL_networks/CSL_pairs.csv` (21,534 pairs) and `TCGA_Drug_response_CSL_pairs.csv` (232) are identical in size to the ISLE FDR-0.2 and drug cSL networks and were not duplicated.

## Not found (no GitHub-hosted copy located)
- SELECT (Lee et al. Nat Commun 2021) SL/SR partner lists for targeted and immune-checkpoint drugs (supplementary tables only; journal/zenodo blocked).
- ENLIGHT (Dinstag et al. Med 2023) genetic-interaction (SL/SR) partner lists per drug: the PangeaResearch/enlight-data repo holds only expression + response labels.
- SynLethDB native download (synlethdb.sist.shanghaitech.edu.cn) and SR/SDL-specific tables (e.g. Sinha 2017 SR, DAISY SDL) — not mirrored on GitHub as far as searched.
