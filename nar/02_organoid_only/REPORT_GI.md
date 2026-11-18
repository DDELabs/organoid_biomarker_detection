# GI organoid-only report (colon, rectum, pancreas, liver)

## What this is
We re-ran the Kong et al. 2020 method using **only lab models** (organoids, and one PDX set) to train, and TCGA
patients only to test. Steps: (1) pick Reactome pathways close to the drug's targets in the STRING protein network
(proximity z <= -1.2816, the paper's cut-off); (2) score these pathways with ssGSEA (a per-sample pathway activity
score); (3) train Ridge, SVR and LinearRegression on lab-model drug response; (4) score TCGA patients who got the
drug, split them at the median prediction, and compare overall survival (Kaplan-Meier + log-rank test). The same
split is applied to patients with **no drug record** ("untreated") as a check for general prognosis.

Extra checks for every set x drug:
- **Bootstrap** (paper Fig. 4e/f): 10,000 random pathway sets of the same size. "Share" = fraction of random sets
  whose rank-1 pathway gives a log-rank p as small or smaller. Small share = the real proximal pathways are special.
- **Response shuffle**: 1,000 runs where the lab-model drug responses are shuffled before training. "rank-1 share"
  compares the main result; "best-of-30 share" compares the best p over all 30 settings (3 algorithms x top 1-10
  pathways), which shows what "picking the best setting" can give by chance.
- **Leave-one-out rho**: how well the model predicts drug response inside the lab models (Spearman correlation).
- **Cut-off sensitivity**: proximity cut-offs -1.1824, -1.2816 (paper), -1.6449, -2.1925.

HR = hazard ratio of predicted non-responders vs predicted responders. HR > 1 means predicted non-responders die
sooner (the expected direction). "Not estimable" = no deaths in one group.

## Colon reproduction of the paper (checked from the pickles, not copied)
Using the authors' own committed files (organoid and TCGA NES from MOESM3, 2019 TCGA drug file):

| Paper item | Paper number | Our number |
|---|---|---|
| Fig 2d left: 5-FU-treated COAD | n = 114 (57/57), p = 0.014, no deaths in responders | n = 114 (57/57), 8 deaths (0 vs 8), p = 0.014, rank-1 = Activation of BH3-only proteins |
| Fig 2d right: untreated COAD | n = 298, p = 0.16 | n = 298, 41 deaths, p = 0.161, HR 0.64 (0.34-1.20) |
| Fig 4b: HCT-8 sensitive vs resistant (BH3-only) | t-test p = 0.0085, sensitive higher | p = 0.0085 (t = 3.27), sensitive higher |
| Fig 4e: top-rank bootstrap (paper statistic) | 12 / 10,000 random sets reach rank 1 (p = 0.0012); 8,927 none | 74 / 10,000 reach rank 1 (p = 0.0074); 7,209 none |
| Fig 5a: cetuximab model vs BRAF V600E | one-sided MWU p = 0.037 (WT 328, V600E 47) | p = 0.031 (WT 369, V600E 49), Gastrin-CREB via PKC and MAPK |
| Supp Table 1: clinical features | age p = 0.26, sex p = 0.19, T p = 0.44, N p = 0.44 | age 0.260, sex 0.188, T 0.435, N 0.440 (same counts) |
| Supp Fig 10a: pathway size vs p | R2 = 0.0026 | R2 = 0.012 |

What I think: the paper's colon 5-FU result reproduces exactly with the authors' inputs. Our bootstrap gives a
larger share than the paper (0.0074 vs 0.0012); it is still small. The response-shuffle test adds a view the paper
did not have: 3.4% of shuffled runs give the rank-1 p as small, and 16.4% when the best of 30 settings is picked.
With **current** GDC data (new expression pipeline, more patients, more follow-up: COAD 5-FU n = 129, 25 deaths) the
rank-1 pathway changes (Purine ribonucleoside monophosphate biosynthesis) and p = 0.866. With current data the
cetuximab-BRAF test gets stronger (one-sided p = 0.00072).

## Colon extra analyses
- **5-FU regimen split** (`COAD_vandewetering/figures/Fig2d_regimen_split_5FU_vs_FOLFOX`, table
  `COAD_vandewetering/tables/regimen_split_5FU_FOLFOX.tsv`). Authors' files: 94 of 114 "5-FU" patients also had
  oxaliplatin (FOLFOX-type): p = 0.036, 6 deaths, all in predicted non-responders. Only 11 had no oxaliplatin or
  irinotecan record (0 deaths) and 9 had irinotecan without oxaliplatin (2 deaths). Current data, COAD: FOLFOX-type
  n = 106, 18 deaths, p = 0.583, HR 1.30 (0.51-3.29). What I think: the paper's "5-FU" signal sits almost entirely in
  FOLFOX patients; a true 5-FU-alone group is too small in TCGA to test.
- **Fig 4b isogenic cells, current data**: recomputed ssGSEA on GSE81005 gives BH3-only p = 0.039 (sensitive higher,
  same direction); purine ribonucleoside pathway p = 0.016. Samples (audit I24, from `code/iso.py`, as the paper's
  Fig. 4b): sensitive = parental HCT-8 at 0 h (3 samples, titles "HCT8-0h*"); resistant = HCT-8/5-FU at 0 h, 24 h and
  48 h pooled (9 samples, "HCT8/5Fu-0h/24h/48h*"); two-sided t-test, 3 vs 9. nar/08 instead compares 3 vs 3 at each
  time point (0 h p = 0.18, 24 h p = 0.007, 48 h p = 0.019), so the two p-values answer different questions.
- **Supp Fig 15 (GDSC cell lines -> TCGA COAD, 5-FU)**: all 903 GDSC1 lines, rank-1 CDC6/ORC: p = 0.708 (current), 0.583
  (authors' patients); 43 colorectal lines, rank-1 E2F-enabled inhibition: p = 0.403 / 0.506. Paper: 0.58 and 0.51.
  What I think: as in the paper, cell-line models do not give the organoid result; this supports the organoid angle.
- **Fig 5a, PDX (Isella)**: Signaling by NOTCH3, one-sided p = 0.266.

## Sets x drugs: all numbers
Below, one block per set x drug (generated by `code/report_sections.py` from `tables/cross_set_summary.tsv`).
Rank-1 Ridge pathway, median split, overall survival.

### Colorectal organoids (van de Wetering 2015)  (`COAD_vandewetering/`)

#### TCGA colon, authors' 2019 files (exact paper inputs) | fluorouracil
- Lab models: 19 (response = ln IC50); proximal pathways: 37; rank-1 Ridge pathway: Activation of bh3 only proteins.
- Drug-treated patients: n = 114, deaths = 8 (0 in predicted responders, 8 in predicted non-responders); log-rank p = 0.014; HR not estimable (no deaths in one group).
- Untreated patients (same split): n = 298, deaths = 41; p = 0.161; HR 0.64 (95% CI 0.34-1.20).
- Other algorithms (rank-1): SVR p = 0.014, LinearRegression p = 0.014. Best of all 30 settings: p = 0.014 (Ridge k=1); settings with p < 0.05: 3 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.007. Response shuffle (1,000 runs): rank-1 share = 0.034; best-of-30 share = 0.164 (best setting HR not estimable, expected direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.44 (p = 0.063), SVR 0.47, LinearRegression 0.49.
- What I think: The pathway splits drug-treated patients into groups with different survival, and the predicted non-responders die more often, as expected. In the response-shuffle test 3.4% of shuffled runs did as well, and 0.7% of random pathway sets did as well. Untreated patients are not separated, which fits a drug-specific signal if one exists. Only 8 deaths among treated patients, so any estimate here is fragile.

#### TCGA colon, authors' 2019 files (exact paper inputs) | irinotecan
- Lab models: 19 (response = ln IC50); proximal pathways: 76; rank-1 Ridge pathway: Hyaluronan uptake and degradation.
- Drug-treated patients: n = 22, deaths = 3 (3 in predicted responders, 0 in predicted non-responders); log-rank p = 0.182; HR not estimable (no deaths in one group).
- Untreated patients (same split): n = 298, deaths = 41; p = 0.718; HR 0.89 (95% CI 0.48-1.67).
- Other algorithms (rank-1): SVR p = 0.182, LinearRegression p = 0.182. Best of all 30 settings: p = 0.083 (Ridge k=3); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.251. Response shuffle (1,000 runs): rank-1 share = 0.317; best-of-30 share = 0.495 (best setting HR not estimable, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.21 (p = 0.379), SVR -0.22, LinearRegression -0.21.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Only 3 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model does not predict drug response (leave-one-out rho = -0.21, negative).

#### TCGA colon, authors' 2019 files (exact paper inputs) | oxaliplatin
- Lab models: 19 (response = ln IC50); proximal pathways: 43; rank-1 Ridge pathway: Bile acid and bile salt metabolism.
- Drug-treated patients: n = 100, deaths = 7 (2 in predicted responders, 5 in predicted non-responders); log-rank p = 0.844; HR 1.18 (95% CI 0.23-6.11).
- Untreated patients (same split): n = 298, deaths = 41; p = 0.660; HR 0.87 (95% CI 0.47-1.61).
- Other algorithms (rank-1): SVR p = 0.327, LinearRegression p = 0.844. Best of all 30 settings: p = 0.127 (LinearRegression k=5); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.844. Response shuffle (1,000 runs): rank-1 share = 0.946; best-of-30 share = 0.569 (best setting HR 4.49, expected direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.15 (p = 0.533), SVR 0.16, LinearRegression 0.15.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Only 7 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.15).

#### TCGA colon (current GDC data) | cetuximab
- Lab models: 19 (response = ln IC50); proximal pathways: 32; rank-1 Ridge pathway: Gastrin creb signalling pathway via pkc and mapk.
- Drug-treated patients: n = 15, deaths = 7 (4 in predicted responders, 3 in predicted non-responders); log-rank p = 0.329; HR 0.44 (95% CI 0.08-2.42).
- Untreated patients (same split): n = 302, deaths = 68; p = 0.193; HR 1.37 (95% CI 0.85-2.21).
- Other algorithms (rank-1): SVR p = 0.329, LinearRegression p = 0.329. Best of all 30 settings: p = 0.018 (Ridge k=10); settings with p < 0.05: 3 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.403. Response shuffle (1,000 runs): rank-1 share = 0.539; best-of-30 share = 0.311 (best setting HR not estimable, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.40 (p = 0.090), SVR -0.33, LinearRegression -0.32.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Some of the other 29 model settings reach p < 0.05 (3 of 30), but 31.1% of shuffled runs find a best p at least this small, so picking the best setting is not enough on its own. Only 7 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model does not predict drug response (leave-one-out rho = -0.40, negative).

#### TCGA colon (current GDC data) | fluorouracil
- Lab models: 19 (response = ln IC50); proximal pathways: 30; rank-1 Ridge pathway: Purine ribonucleoside monophosphate biosynthesis.
- Drug-treated patients: n = 129, deaths = 25 (12 in predicted responders, 13 in predicted non-responders); log-rank p = 0.866; HR 1.07 (95% CI 0.49-2.35).
- Untreated patients (same split): n = 302, deaths = 68; p = 0.419; HR 1.22 (95% CI 0.75-1.96).
- Other algorithms (rank-1): SVR p = 0.866, LinearRegression p = 0.866. Best of all 30 settings: p = 0.179 (SVR k=4); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.885. Response shuffle (1,000 runs): rank-1 share = 0.941; best-of-30 share = 0.790 (best setting HR 0.58, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.58 (p = 0.009), SVR 0.61, LinearRegression 0.53.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either.

#### TCGA colon (current GDC data) | irinotecan
- Lab models: 19 (response = ln IC50); proximal pathways: 66; rank-1 Ridge pathway: Hyaluronan uptake and degradation.
- Drug-treated patients: n = 32, deaths = 15 (9 in predicted responders, 6 in predicted non-responders); log-rank p = 0.257; HR 0.55 (95% CI 0.19-1.56).
- Untreated patients (same split): n = 302, deaths = 68; p = 0.043; HR 0.60 (95% CI 0.36-0.99).
- Other algorithms (rank-1): SVR p = 0.257, LinearRegression p = 0.257. Best of all 30 settings: p = 0.257 (Ridge k=1); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.199. Response shuffle (1,000 runs): rank-1 share = 0.238; best-of-30 share = 0.787 (best setting HR 0.55, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.22 (p = 0.355), SVR -0.21, LinearRegression -0.19.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. The same split also separates untreated patients, so part of the signal may be general prognosis rather than drug response. Inside the lab models the pathway model does not predict drug response (leave-one-out rho = -0.22, negative).

#### TCGA colon (current GDC data) | oxaliplatin
- Lab models: 19 (response = ln IC50); proximal pathways: 36; rank-1 Ridge pathway: Metabolism of polyamines.
- Drug-treated patients: n = 112, deaths = 21 (12 in predicted responders, 9 in predicted non-responders); log-rank p = 0.440; HR 0.71 (95% CI 0.30-1.70).
- Untreated patients (same split): n = 302, deaths = 68; p = 0.470; HR 1.19 (95% CI 0.74-1.93).
- Other algorithms (rank-1): SVR p = 0.440, LinearRegression p = 0.880. Best of all 30 settings: p = 0.001 (Ridge k=4); settings with p < 0.05: 14 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.511. Response shuffle (1,000 runs): rank-1 share = 0.448; best-of-30 share = 0.008 (best setting HR 0.17, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.09 (p = 0.700), SVR 0.06, LinearRegression 0.06.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Some of the other 29 model settings reach p < 0.05 (14 of 30). The best setting (Ridge k=4, p = 0.001) beats 99.2% of shuffled runs, but its HR is 0.17 (opposite direction): predicted responders do worse, the opposite of the drug hypothesis. No setting points in the expected direction. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.09).

#### TCGA rectum (current GDC data) | fluorouracil
- Lab models: 19 (response = ln IC50); proximal pathways: 30; rank-1 Ridge pathway: Purine ribonucleoside monophosphate biosynthesis.
- Drug-treated patients: n = 59, deaths = 5 (2 in predicted responders, 3 in predicted non-responders); log-rank p = 0.409; HR 2.12 (95% CI 0.34-13.17).
- Untreated patients (same split): n = 89, deaths = 18; p = 0.024; HR 3.34 (95% CI 1.09-10.20).
- Other algorithms (rank-1): SVR p = 0.409, LinearRegression p = 0.409. Best of all 30 settings: p = 0.062 (SVR k=6); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.398. Response shuffle (1,000 runs): rank-1 share = 0.219; best-of-30 share = 0.536 (best setting HR 7.28, expected direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.58 (p = 0.009), SVR 0.61, LinearRegression 0.53.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. The same split also separates untreated patients, so part of the signal may be general prognosis rather than drug response. Only 5 deaths among treated patients, so any estimate here is fragile.

#### TCGA rectum (current GDC data) | oxaliplatin
- Lab models: 19 (response = ln IC50); proximal pathways: 36; rank-1 Ridge pathway: Metabolism of polyamines.
- Drug-treated patients: n = 42, deaths = 5 (1 in predicted responders, 4 in predicted non-responders); log-rank p = 0.093; HR 5.37 (95% CI 0.60-48.34).
- Untreated patients (same split): n = 89, deaths = 18; p = 0.326; HR 0.62 (95% CI 0.23-1.65).
- Other algorithms (rank-1): SVR p = 0.093, LinearRegression p = 0.133. Best of all 30 settings: p = 0.093 (Ridge k=1); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.084. Response shuffle (1,000 runs): rank-1 share = 0.126; best-of-30 share = 0.391 (best setting HR 5.37, expected direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.09 (p = 0.700), SVR 0.06, LinearRegression 0.06.
- What I think: There is a weak trend in drug-treated patients. Untreated patients are not separated, which fits a drug-specific signal if one exists. Only 5 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.09).

#### TCGA colon + rectum pooled (current GDC data) | fluorouracil
- Lab models: 19 (response = ln IC50); proximal pathways: 30; rank-1 Ridge pathway: Purine ribonucleoside monophosphate biosynthesis.
- Drug-treated patients: n = 188, deaths = 30 (14 in predicted responders, 16 in predicted non-responders); log-rank p = 0.603; HR 1.21 (95% CI 0.59-2.48).
- Untreated patients (same split): n = 391, deaths = 86; p = 0.104; HR 1.42 (95% CI 0.93-2.19).
- Other algorithms (rank-1): SVR p = 0.603, LinearRegression p = 0.603. Best of all 30 settings: p = 0.343 (SVR k=4); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.564. Response shuffle (1,000 runs): rank-1 share = 0.775; best-of-30 share = 0.978 (best setting HR 0.71, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.58 (p = 0.009), SVR 0.61, LinearRegression 0.53.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either.

#### TCGA colon + rectum pooled (current GDC data) | irinotecan
- Lab models: 19 (response = ln IC50); proximal pathways: 66; rank-1 Ridge pathway: Hyaluronan uptake and degradation.
- Drug-treated patients: n = 39, deaths = 17 (10 in predicted responders, 7 in predicted non-responders); log-rank p = 0.392; HR 0.65 (95% CI 0.25-1.74).
- Untreated patients (same split): n = 391, deaths = 86; p = 0.028; HR 0.61 (95% CI 0.39-0.95).
- Other algorithms (rank-1): SVR p = 0.392, LinearRegression p = 0.392. Best of all 30 settings: p = 0.392 (Ridge k=1); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.439. Response shuffle (1,000 runs): rank-1 share = 0.703; best-of-30 share = 0.983 (best setting HR 0.65, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.22 (p = 0.355), SVR -0.21, LinearRegression -0.19.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. The same split also separates untreated patients, so part of the signal may be general prognosis rather than drug response. Inside the lab models the pathway model does not predict drug response (leave-one-out rho = -0.22, negative).

#### TCGA colon + rectum pooled (current GDC data) | oxaliplatin
- Lab models: 19 (response = ln IC50); proximal pathways: 36; rank-1 Ridge pathway: Metabolism of polyamines.
- Drug-treated patients: n = 154, deaths = 26 (14 in predicted responders, 12 in predicted non-responders); log-rank p = 0.695; HR 0.86 (95% CI 0.40-1.86).
- Untreated patients (same split): n = 391, deaths = 86; p = 0.975; HR 0.99 (95% CI 0.65-1.52).
- Other algorithms (rank-1): SVR p = 0.695, LinearRegression p = 0.911. Best of all 30 settings: p = 0.00046 (LinearRegression k=6); settings with p < 0.05: 21 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.802. Response shuffle (1,000 runs): rank-1 share = 0.822; best-of-30 share = 0.007 (best setting HR 0.18, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.09 (p = 0.700), SVR 0.06, LinearRegression 0.06.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Some of the other 29 model settings reach p < 0.05 (21 of 30). The best setting (LinearRegression k=6, p = 0.00046) beats 99.3% of shuffled runs, but its HR is 0.18 (opposite direction): predicted responders do worse, the opposite of the drug hypothesis. The best setting in the expected direction has p = 0.911; at most 100.0% of shuffles do as well (upper bound). Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.09).

### Pancreatic organoids (Tiriac 2018)  (`PAAD_tiriac/`)

#### TCGA pancreas | fluorouracil
- Lab models: 36 (response = AUC); proximal pathways: 38; rank-1 Ridge pathway: Tetrahydrobiopterin bh4 synthesis recycling salvage and regulation.
- Drug-treated patients: n = 37, deaths = 18 (10 in predicted responders, 8 in predicted non-responders); log-rank p = 0.932; HR 0.96 (95% CI 0.38-2.44).
- Untreated patients (same split): n = 61, deaths = 35; p = 0.281; HR 0.69 (95% CI 0.34-1.37).
- Other algorithms (rank-1): SVR p = 0.420, LinearRegression p = 0.793. Best of all 30 settings: p = 0.277 (Ridge k=4); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.926. Response shuffle (1,000 runs): rank-1 share = 0.873; best-of-30 share = 0.913 (best setting HR 0.58, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.01 (p = 0.948), SVR -0.04, LinearRegression -0.19.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.01).

#### TCGA pancreas | gemcitabine
- Lab models: 36 (response = AUC); proximal pathways: 45; rank-1 Ridge pathway: Homologous recombination repair of replication independent double strand breaks.
- Drug-treated patients: n = 104, deaths = 51 (23 in predicted responders, 28 in predicted non-responders); log-rank p = 0.223; HR 1.41 (95% CI 0.81-2.46).
- Untreated patients (same split): n = 61, deaths = 35; p = 0.897; HR 1.05 (95% CI 0.54-2.05).
- Other algorithms (rank-1): SVR p = 0.033, LinearRegression p = 0.151. Best of all 30 settings: p = 0.013 (SVR k=5); settings with p < 0.05: 5 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.296. Response shuffle (1,000 runs): rank-1 share = 0.805; best-of-30 share = 0.388 (best setting HR 2.01, expected direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.37 (p = 0.025), SVR 0.47, LinearRegression 0.28.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Some of the other 29 model settings reach p < 0.05 (5 of 30), but 38.8% of shuffled runs find a best p at least this small, so picking the best setting is not enough on its own.

#### TCGA pancreas | oxaliplatin
- Lab models: 36 (response = AUC); proximal pathways: 47; rank-1 Ridge pathway: Synthesis of substrates in n glycan biosythesis.
- Drug-treated patients: n = 21, deaths = 11 (8 in predicted responders, 3 in predicted non-responders); log-rank p = 0.223; HR 0.44 (95% CI 0.11-1.71).
- Untreated patients (same split): n = 61, deaths = 35; p = 0.072; HR 0.53 (95% CI 0.26-1.07).
- Other algorithms (rank-1): SVR p = 1.000, LinearRegression p = 0.636. Best of all 30 settings: p = 0.201 (LinearRegression k=10); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.160. Response shuffle (1,000 runs): rank-1 share = 0.135; best-of-30 share = 0.645 (best setting HR 0.42, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.14 (p = 0.399), SVR -0.40, LinearRegression 0.10.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Only 11 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.14).

#### TCGA pancreas | sn-38  (12-14 treated patients: borderline size)
- Lab models: 36 (response = AUC); proximal pathways: 78; rank-1 Ridge pathway: Regulation of the fanconi anemia pathway.
- Drug-treated patients: n = 13, deaths = 5 (3 in predicted responders, 2 in predicted non-responders); log-rank p = 0.570; HR 0.59 (95% CI 0.10-3.63).
- Untreated patients (same split): n = 61, deaths = 35; p = 0.186; HR 1.58 (95% CI 0.80-3.11).
- Other algorithms (rank-1): SVR p = 0.680, LinearRegression p = 0.570. Best of all 30 settings: p = 0.023 (Ridge k=9); settings with p < 0.05: 4 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.610. Response shuffle (1,000 runs): rank-1 share = 0.626; best-of-30 share = 0.189 (best setting HR not estimable, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.02 (p = 0.915), SVR 0.19, LinearRegression 0.00.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Some of the other 29 model settings reach p < 0.05 (4 of 30), but 18.9% of shuffled runs find a best p at least this small, so picking the best setting is not enough on its own. Only 5 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.02).

### Pancreatic organoids (Shi 2022)  (`PAAD_shi/`)

#### TCGA pancreas | fluorouracil
- Lab models: 39 (response = AUC); proximal pathways: 38; rank-1 Ridge pathway: Inhibition of replication initiation of damaged dna by rb1 e2f1.
- Drug-treated patients: n = 37, deaths = 18 (10 in predicted responders, 8 in predicted non-responders); log-rank p = 0.958; HR 0.98 (95% CI 0.38-2.48).
- Untreated patients (same split): n = 61, deaths = 35; p = 0.412; HR 0.75 (95% CI 0.39-1.48).
- Other algorithms (rank-1): SVR p = 0.696, LinearRegression p = 0.400. Best of all 30 settings: p = 0.097 (Ridge k=3); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.947. Response shuffle (1,000 runs): rank-1 share = 0.945; best-of-30 share = 0.533 (best setting HR 2.24, expected direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.01 (p = 0.928), SVR 0.08, LinearRegression 0.47.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Inside the lab models the pathway model does not predict drug response (leave-one-out rho = -0.01, negative).

#### TCGA pancreas | gemcitabine
- Lab models: 39 (response = AUC); proximal pathways: 45; rank-1 Ridge pathway: G1 s specific transcription.
- Drug-treated patients: n = 104, deaths = 51 (32 in predicted responders, 19 in predicted non-responders); log-rank p = 0.023; HR 0.52 (95% CI 0.30-0.92).
- Untreated patients (same split): n = 61, deaths = 35; p = 0.017; HR 0.45 (95% CI 0.23-0.88).
- Other algorithms (rank-1): SVR p = 0.146, LinearRegression p = 0.018. Best of all 30 settings: p = 0.005 (Ridge k=8); settings with p < 0.05: 6 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.022. Response shuffle (1,000 runs): rank-1 share = 0.100; best-of-30 share = 0.146 (best setting HR 0.44, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.08 (p = 0.643), SVR 0.04, LinearRegression 0.09.
- What I think: The pathway splits drug-treated patients into groups with different survival, but the predicted responders are the ones who die more often (the opposite direction). In the response-shuffle test 10.0% of shuffled runs did as well, and 2.2% of random pathway sets did as well. The same split also separates untreated patients, so part of the signal may be general prognosis rather than drug response. Inside the lab models the pathway model does not predict drug response (leave-one-out rho = -0.08, negative).

#### TCGA pancreas | irinotecan  (12-14 treated patients: borderline size)
- Lab models: 39 (response = AUC); proximal pathways: 78; rank-1 Ridge pathway: Fanconi anemia pathway.
- Drug-treated patients: n = 13, deaths = 5 (4 in predicted responders, 1 in predicted non-responders); log-rank p = 0.302; HR 0.33 (95% CI 0.04-3.00).
- Untreated patients (same split): n = 61, deaths = 35; p = 0.021; HR 2.22 (95% CI 1.11-4.43).
- Other algorithms (rank-1): SVR p = 0.044, LinearRegression p = 0.302. Best of all 30 settings: p = 0.030 (SVR k=2); settings with p < 0.05: 3 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.310. Response shuffle (1,000 runs): rank-1 share = 0.351; best-of-30 share = 0.305 (best setting HR not estimable, expected direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.27 (p = 0.102), SVR -0.08, LinearRegression -0.25.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. The same split also separates untreated patients, so part of the signal may be general prognosis rather than drug response. Some of the other 29 model settings reach p < 0.05 (3 of 30), but 30.5% of shuffled runs find a best p at least this small, so picking the best setting is not enough on its own. Only 5 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model does not predict drug response (leave-one-out rho = -0.27, negative).

#### TCGA pancreas | oxaliplatin
- Lab models: 39 (response = AUC); proximal pathways: 46; rank-1 Ridge pathway: Intrinsic pathway.
- Drug-treated patients: n = 21, deaths = 11 (7 in predicted responders, 4 in predicted non-responders); log-rank p = 0.497; HR 0.65 (95% CI 0.18-2.30).
- Untreated patients (same split): n = 61, deaths = 35; p = 0.242; HR 1.49 (95% CI 0.76-2.91).
- Other algorithms (rank-1): SVR p = 0.564, LinearRegression p = 0.460. Best of all 30 settings: p = 0.009 (LinearRegression k=10); settings with p < 0.05: 1 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.486. Response shuffle (1,000 runs): rank-1 share = 0.277; best-of-30 share = 0.027 (best setting HR 5.43, expected direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.35 (p = 0.028), SVR -0.07, LinearRegression -0.18.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Some of the other 29 model settings reach p < 0.05 (1 of 30). The best setting (LinearRegression k=10, p = 0.009) beats 97.3% of shuffled runs and points in the expected direction (predicted non-responders do worse). The rank-1 setting the paper's method would use does not reach p < 0.05. Only 11 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model does not predict drug response (leave-one-out rho = -0.35, negative).

### Liver cancer organoids (Ji 2023)  (`LIHC_ji/`)

#### TCGA liver | sorafenib
- Lab models: 61 (response = AUC); proximal pathways: 22; rank-1 Ridge pathway: Signaling by activated point mutants of fgfr1.
- Drug-treated patients: n = 30, deaths = 15 (8 in predicted responders, 7 in predicted non-responders); log-rank p = 0.188; HR 0.48 (95% CI 0.15-1.47).
- Untreated patients (same split): n = 328, deaths = 108; p = 0.129; HR 1.34 (95% CI 0.92-1.96).
- Other algorithms (rank-1): SVR p = 0.171, LinearRegression p = 0.188. Best of all 30 settings: p = 0.013 (LinearRegression k=8); settings with p < 0.05: 4 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.397. Response shuffle (1,000 runs): rank-1 share = 0.583; best-of-30 share = 0.401 (best setting HR 0.24, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.16 (p = 0.215), SVR -0.38, LinearRegression -0.16.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Some of the other 29 model settings reach p < 0.05 (4 of 30), but 40.1% of shuffled runs find a best p at least this small, so picking the best setting is not enough on its own. Inside the lab models the pathway model does not predict drug response (leave-one-out rho = -0.16, negative).

### Liver cancer organoids (Broutier 2017)  (`LIHC_broutier/`)

#### TCGA liver | sorafenib
- Lab models: 5 (response = AUC); proximal pathways: 15; rank-1 Ridge pathway: Reversible hydration of carbon dioxide.
- Drug-treated patients: n = 30, deaths = 15 (8 in predicted responders, 7 in predicted non-responders); log-rank p = 0.139; HR 0.43 (95% CI 0.14-1.35).
- Untreated patients (same split): n = 328, deaths = 108; p = 0.374; HR 1.19 (95% CI 0.81-1.73).
- Other algorithms (rank-1): SVR p = 0.879, LinearRegression p = 0.139. Best of all 30 settings: p = 0.093 (SVR k=4); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.364. Response shuffle (1,000 runs): rank-1 share = 0.521; best-of-30 share = 0.828 (best setting HR 0.38, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.70 (p = 0.188), SVR -0.70, LinearRegression -0.70.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated either. Inside the lab models the pathway model does not predict drug response (leave-one-out rho = -0.70, negative).

### Colorectal PDX (Isella 2017, cetuximab)  (`COAD_isella_pdx/`)

#### TCGA colon (current GDC data) | cetuximab
- Lab models: 192 (response = tumour volume change (3 wk)); proximal pathways: 41; rank-1 Ridge pathway: Signaling by notch3.
- Drug-treated patients: n = 15, deaths = 7 (4 in predicted responders, 3 in predicted non-responders); log-rank p = 0.648; HR 0.70 (95% CI 0.15-3.22).
- Untreated patients (same split): n = 302, deaths = 68; p = 0.046; HR 1.63 (95% CI 1.01-2.65).
- Other algorithms (rank-1): SVR p = 0.648, LinearRegression p = 0.648. Best of all 30 settings: p = 0.370 (LinearRegression k=7); settings with p < 0.05: 0 of 30.
- Rank-1 p bootstrap (10,000 random pathway sets of the same size; not the paper's top-rank statistic): share with p as small or smaller = 0.831. Response shuffle (1,000 runs): rank-1 share = 0.760; best-of-30 share = 1.000 (best setting HR 0.44, opposite direction).
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.50 (p = 8.2e-14), SVR 0.53, LinearRegression 0.50.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. The same split also separates untreated patients, so part of the signal may be general prognosis rather than drug response. Only 7 deaths among treated patients, so any estimate here is fragile.

### Combinations that could not be analysed

| set | cohort | drug | reason |
|---|---|---|---|
| Colorectal organoids (van de Wetering 2015) | READ | irinotecan | 7 treated patients with survival and expression (< 12) |
| Pancreatic organoids (Shi 2022) | PAAD | paclitaxel | 5 treated patients with survival and expression (< 12) |
| Liver cancer organoids (Ji 2023) | LIHC | doxorubicin | 3 treated patients with survival and expression (< 12) |
| Liver cancer organoids (Ji 2023) | LIHC | lenvatinib | no drug response or no proximity |
| Liver cancer organoids (Broutier 2017) | LIHC | doxorubicin | 3 treated patients with survival and expression (< 12) |
| Liver cancer organoids (Broutier 2017) | LIHC | fluorouracil | 1 treated patients with survival and expression (< 12) |
| Colorectal organoids (van de Wetering 2015) | COAD_authors | cetuximab | 9 treated patients with survival and expression (< 12) |
| Colorectal organoids (van de Wetering 2015) | READ | cetuximab | 0 treated patients with survival and expression (< 12) |
| Pancreatic organoids (Tiriac 2018) | PAAD | paclitaxel | 5 treated patients with survival and expression (< 12) |
| Pancreatic organoids (Shi 2022) | PAAD | mitomycin | no drug response or no proximity |
| Liver cancer organoids (Ji 2023) | LIHC | gemcitabine | 7 treated patients with survival and expression (< 12) |
| Liver cancer organoids (Ji 2023) | LIHC | fluorouracil | 1 treated patients with survival and expression (< 12) |
| Liver cancer organoids (Ji 2023) | LIHC | regorafenib | no drug response or no proximity |
| Liver cancer organoids (Broutier 2017) | LIHC | gemcitabine | 7 treated patients with survival and expression (< 12) |
| Liver cancer organoids (Broutier 2017) | LIHC | cisplatin | 4 treated patients with survival and expression (< 12) |

## Overall picture (all sets)
- 23 set x drug x cohort analyses ran (`tables/cross_set_summary.tsv`, `figures/cross_set_forest`); 15 combinations
  had too few treated patients or no drug data (table above).
- Treated-patient p < 0.05 at the paper's settings: colon 5-FU with the authors' files (p = 0.014, expected
  direction) and pancreas gemcitabine with Shi organoids (p = 0.023, but predicted responders died more, and
  untreated patients split the same way, p = 0.017). All other analyses give p between 0.09 and 0.96.
- Many untreated-patient splits reach p < 0.05 (e.g. COAD irinotecan 0.043, READ 5-FU 0.024, Isella cetuximab
  0.046). What I think: pathway scores often carry general prognosis; the untreated check is needed for every claim.
- Pancreas: Tiriac and Shi organoids, trained on the same drugs, pick different rank-1 pathways and give different
  results (gemcitabine p = 0.223 vs 0.023). Liver sorafenib: Ji p = 0.188, Broutier p = 0.139, 30 patients each,
  15 deaths, both HR < 1. What I think: these sets give honest, testable numbers, but patient counts are modest.
- Best of 30 settings (audit fix I03, `tables/cross_set_best_direction.tsv`, made by
  `nar/15_fixes/code/f02_best_setting_direction.py`): for colon oxaliplatin (current data), colon + rectum oxaliplatin
  and Shi pancreas oxaliplatin, the best of the 30 settings beats 99.2%, 99.3% and 97.3% of response shuffles. For the two
  colon ones that best setting points the OPPOSITE way (HR 0.17 and 0.18: predicted responders do worse) and no
  setting reaches p < 0.05 in the expected direction. For Shi pancreas oxaliplatin the best setting (LinearRegression
  k = 10, p = 0.009) points the expected way (HR 5.4), but the paper's rank-1 setting does not reach p < 0.05.
  What I think: a "best of 30" search can find strong splits, but only one of the three points the way the drug
  hypothesis predicts, and it is not the setting the method pre-specifies.
- Two different "bootstrap" numbers exist (audit fix I15). "Top-rank bootstrap (paper Fig 4e)" = share of random
  pathway sets whose best pathway reaches rank 1 (colon 0.0074, bladder 0.0437). "Rank-1 p bootstrap" = share of
  random sets whose rank-1 p is as small (colon 0.0072, bladder 0.0663). Quote the first when comparing with the paper.
- Survival sources (audit I14). This GI report uses TCGA-CDR survival and Efron ties; the BLCA/OV/LUAD/OTHER report
  (Novartis PDX sets) uses Xena GDC survival, a different drug-record parser and Breslow ties. So the same TCGA cohort
  can differ by a few patients: COAD cetuximab 15 treated / 7 deaths here vs 14 / 6 there; COAD untreated 302 vs 282;
  PAAD untreated 61 vs 60. Quote the GI numbers (TCGA-CDR) for COAD and PAAD; the PDX-set numbers were not re-run.
- Full cut-off sensitivity: `tables/report_numbers.md` (second table) and `figures/cross_set_cutoff_sensitivity`.

## Figures (PNG + PDF + caption.md each, 300 dpi)
Per set folder (`COAD_vandewetering/`, `PAAD_tiriac/`, `PAAD_shi/`, `LIHC_ji/`, `LIHC_broutier/`, `COAD_isella_pdx/`),
`figures/` holds, for each analysed cohort x drug:
- `Fig2a_Fig2d_Fig4e_SuppFig17_panels_<cohort>_<drug>`: the full story for one drug: proximal pathways, cut-off
  bars, pathway heatmap, leave-one-out predictions for 3 algorithms, top-k curve, survival curves for treated and
  untreated patients, bootstrap and shuffle histograms, number box. (Cetuximab version also = Supp Fig 11.)
- `Fig2b_SuppFig2_SuppFig3_coefficients_<cohort>_<drug>`: pathway weights from Ridge, SVR and LinearRegression and how
  they agree with each other.
- `SuppFig5_SourceData-multipathway_KMgrid_<cohort>_<drug>`: survival curves when the top 1 to 10 pathways are combined.
- `SuppFig10_size_vs_p_<cohort>_<drug>`: pathway size vs survival p, to show size does not drive the result.
- `Fig2d_main_composite_<set>`: one page with the treated and untreated survival curves for every drug of the set.
- `supp_proximity_heatmap`: proximity z of all pathways for all drugs of the set.
- `SuppFig1_PCA_models_vs_TCGA`: lab models vs patients before and after per-dataset scaling (batch check).

Colon only (`COAD_vandewetering/figures/`):
- `Fig2d_supp_composite_READ_and_pooled`: same composite for rectum and colon + rectum pooled.
- `Fig2d_regimen_split_5FU_vs_FOLFOX`: 5-FU patients split by whether they also got oxaliplatin or irinotecan.
- `Fig4b_SuppTable3_isogenic_HCT8_5FU`: top pathways in 5-FU-sensitive vs -resistant HCT-8 cells (authors' and recomputed scores).
- `Fig5a_cetuximab_BRAF_V600E`: predicted cetuximab resistance in BRAF V600E vs wild-type TCGA colon tumours (organoids and PDX).
- `SuppFig15_GDSC_cell_line_models_5FU`: 5-FU models trained on GDSC cell lines instead of organoids, applied to TCGA colon.

Cross-set (`figures/`):
- `cross_set_forest`: hazard ratios with 95% CI for treated and untreated patients, every set x drug.
- `cross_set_cutoff_sensitivity`: treated p for every set x drug at 4 proximity cut-offs.

## Tables
- `tables/cross_set_summary.tsv` (all numbers above), `tables/report_numbers.md`, `tables/report_sections.md`,
  `tables/cross_set_cutoff_sensitivity.tsv`, `tables/SuppTable1_clinical_features_by_predicted_group.tsv`.
- Per set `tables/`: all 30 tests, coefficients, leave-one-out predictions, proximity, cut-off sensitivity,
  bootstrap and shuffle p-value lists (gz), bootstrap summaries (Supp Table 4 analogue).
- Colon: `Fig4b_SuppTable3_isogenic_HCT8_5FU.tsv`, `Fig5a_cetuximab_BRAF_V600E.tsv`,
  `SuppFig15_GDSC_cell_line_models_5FU.tsv`, `regimen_split_5FU_FOLFOX.tsv`, `proximity_all_drugs.tsv`.

## Checklist IDs covered (nar/00_checklist/checklist.tsv, search_token)
Fig2a, Fig2b, Fig2d (left and right), Fig4b, Fig4e, Fig5a, SuppFig1, SuppFig2a, SuppFig2b, SuppFig3, SuppFig5a,
SourceData-multipathway, SuppFig10a, SuppFig11, SuppFig15a, SuppFig15b, SuppFig17a (as leave-one-out prediction
inside the organoids, not the 60/10/30 split), SuppTable1, SuppTable3, SuppTable4, Text-method-params (paper
cut-off, 10,000 bootstraps). Not in this assignment: bladder items, schematics (Fig1, Fig3a, Fig4a, Fig4d),
Fig3b-f comparison methods, SuppFig6-9, SuppFig12-14, SuppFig16, SuppFig18.

## Data versions (`_data/raw/`, bytes)
STRING v12.0 links 83,164,437 and info 1,970,090; TCGA STAR FPKM-UQ: COAD 128,114,149, READ 44,400,796, PAAD
49,309,161, LIHC 100,179,162; TCGA-CDR survival (Survival_SupplementalTable_S1_20171025) 2,419,504; GDC biotab drug
files COAD 237,155, PAAD 111,466; GSE81005 series matrix 5,104,292, GSE81006 26,244,390; GDSC1 fitted dose response
27Oct23 29,353,210; Cell Model Passports rnaseq_all_20220624.zip 940,707,972 and model_list_20230110 874,186;
derived gdsc_5FU_expression.tsv.gz 98,012,954 (moved to `_data/raw/`). Lab-model responses from
`data/curated/preclinical/`; authors' NES from MOESM3 (`nar/_cache/`). Seeds: bootstrap 0, shuffle 1.

## Code (`code/`)
common.py (loaders), fastridge.py, engine.py (models + survival; logrank now returns NaN when the variance is zero),
run_all.py (all combos -> `_data/results/*.pkl`), plots.py, extras.py, iso.py, report_tables.py, report_sections.py.
Re-run: `python run_all.py 0 2 & python run_all.py 1 2; python plots.py; python extras.py; python report_tables.py;
python report_sections.py`.
