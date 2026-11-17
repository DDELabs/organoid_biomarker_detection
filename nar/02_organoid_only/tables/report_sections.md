### Colorectal organoids (van de Wetering 2015)  (`COAD_vandewetering/`)

#### TCGA colon, authors' 2019 files (exact paper inputs) | fluorouracil
- Lab models: 19 (response = ln IC50); proximal pathways: 37; rank-1 Ridge pathway: Activation of bh3 only proteins.
- Drug-treated patients: n = 114, deaths = 8 (0 in predicted responders, 8 in predicted non-responders); log-rank p = 0.014; HR not estimable (no deaths in one group).
- Untreated patients (same split): n = 298, deaths = 41; p = 0.161; HR 0.64 (95% CI 0.34-1.20).
- Other algorithms (rank-1): SVR p = 0.014, LinearRegression p = 0.014. Best of all 30 settings: p = 0.014 (Ridge k=1); settings with p < 0.05: 3 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.007. Response shuffle (1,000 runs): rank-1 share = 0.034; best-of-30 share = 0.164.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.44 (p = 0.063), SVR 0.47, LinearRegression 0.49.
- What I think: The pathway splits drug-treated patients into groups with different survival, and the predicted non-responders die more often, as expected. In the response-shuffle test 3.4% of shuffled runs did as well, and 0.7% of random pathway sets did as well. Untreated patients are not separated, which fits a drug-specific signal if one exists. Only 8 deaths among treated patients, so any estimate here is fragile.

#### TCGA colon, authors' 2019 files (exact paper inputs) | irinotecan
- Lab models: 19 (response = ln IC50); proximal pathways: 76; rank-1 Ridge pathway: Hyaluronan uptake and degradation.
- Drug-treated patients: n = 22, deaths = 3 (3 in predicted responders, 0 in predicted non-responders); log-rank p = 0.182; HR not estimable (no deaths in one group).
- Untreated patients (same split): n = 298, deaths = 41; p = 0.718; HR 0.89 (95% CI 0.48-1.67).
- Other algorithms (rank-1): SVR p = 0.182, LinearRegression p = 0.182. Best of all 30 settings: p = 0.083 (Ridge k=3); settings with p < 0.05: 0 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.251. Response shuffle (1,000 runs): rank-1 share = 0.317; best-of-30 share = 0.495.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.21 (p = 0.379), SVR -0.22, LinearRegression -0.21.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated, which fits a drug-specific signal if one exists. Only 3 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = -0.21).

#### TCGA colon, authors' 2019 files (exact paper inputs) | oxaliplatin
- Lab models: 19 (response = ln IC50); proximal pathways: 43; rank-1 Ridge pathway: Bile acid and bile salt metabolism.
- Drug-treated patients: n = 100, deaths = 7 (2 in predicted responders, 5 in predicted non-responders); log-rank p = 0.844; HR 1.18 (95% CI 0.23-6.11).
- Untreated patients (same split): n = 298, deaths = 41; p = 0.660; HR 0.87 (95% CI 0.47-1.61).
- Other algorithms (rank-1): SVR p = 0.327, LinearRegression p = 0.844. Best of all 30 settings: p = 0.127 (LinearRegression k=5); settings with p < 0.05: 0 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.844. Response shuffle (1,000 runs): rank-1 share = 0.946; best-of-30 share = 0.569.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.15 (p = 0.533), SVR 0.16, LinearRegression 0.15.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated, which fits a drug-specific signal if one exists. Only 7 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.15).

#### TCGA colon (current GDC data) | cetuximab
- Lab models: 19 (response = ln IC50); proximal pathways: 32; rank-1 Ridge pathway: Gastrin creb signalling pathway via pkc and mapk.
- Drug-treated patients: n = 15, deaths = 7 (4 in predicted responders, 3 in predicted non-responders); log-rank p = 0.329; HR 0.44 (95% CI 0.08-2.42).
- Untreated patients (same split): n = 302, deaths = 68; p = 0.193; HR 1.37 (95% CI 0.85-2.21).
- Other algorithms (rank-1): SVR p = 0.329, LinearRegression p = 0.329. Best of all 30 settings: p = 0.018 (Ridge k=10); settings with p < 0.05: 3 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.403. Response shuffle (1,000 runs): rank-1 share = 0.539; best-of-30 share = 0.311.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.40 (p = 0.090), SVR -0.33, LinearRegression -0.32.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated, which fits a drug-specific signal if one exists. Some of the other 29 model settings reach p < 0.05 (3 of 30), but 31.1% of shuffled runs find a best p at least this small, so picking the best setting is not enough on its own. Only 7 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = -0.40).

#### TCGA colon (current GDC data) | fluorouracil
- Lab models: 19 (response = ln IC50); proximal pathways: 30; rank-1 Ridge pathway: Purine ribonucleoside monophosphate biosynthesis.
- Drug-treated patients: n = 129, deaths = 25 (12 in predicted responders, 13 in predicted non-responders); log-rank p = 0.866; HR 1.07 (95% CI 0.49-2.35).
- Untreated patients (same split): n = 302, deaths = 68; p = 0.419; HR 1.22 (95% CI 0.75-1.96).
- Other algorithms (rank-1): SVR p = 0.866, LinearRegression p = 0.866. Best of all 30 settings: p = 0.179 (SVR k=4); settings with p < 0.05: 0 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.885. Response shuffle (1,000 runs): rank-1 share = 0.941; best-of-30 share = 0.790.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.58 (p = 0.009), SVR 0.61, LinearRegression 0.53.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated, which fits a drug-specific signal if one exists.

#### TCGA colon (current GDC data) | irinotecan
- Lab models: 19 (response = ln IC50); proximal pathways: 66; rank-1 Ridge pathway: Hyaluronan uptake and degradation.
- Drug-treated patients: n = 32, deaths = 15 (9 in predicted responders, 6 in predicted non-responders); log-rank p = 0.257; HR 0.55 (95% CI 0.19-1.56).
- Untreated patients (same split): n = 302, deaths = 68; p = 0.043; HR 0.60 (95% CI 0.36-0.99).
- Other algorithms (rank-1): SVR p = 0.257, LinearRegression p = 0.257. Best of all 30 settings: p = 0.257 (Ridge k=1); settings with p < 0.05: 0 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.199. Response shuffle (1,000 runs): rank-1 share = 0.238; best-of-30 share = 0.787.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.22 (p = 0.355), SVR -0.21, LinearRegression -0.19.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. The same split also separates untreated patients, so part of the signal may be general prognosis rather than drug response. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = -0.22).

#### TCGA colon (current GDC data) | oxaliplatin
- Lab models: 19 (response = ln IC50); proximal pathways: 36; rank-1 Ridge pathway: Metabolism of polyamines.
- Drug-treated patients: n = 112, deaths = 21 (12 in predicted responders, 9 in predicted non-responders); log-rank p = 0.440; HR 0.71 (95% CI 0.30-1.70).
- Untreated patients (same split): n = 302, deaths = 68; p = 0.470; HR 1.19 (95% CI 0.74-1.93).
- Other algorithms (rank-1): SVR p = 0.440, LinearRegression p = 0.880. Best of all 30 settings: p = 0.001 (Ridge k=4); settings with p < 0.05: 14 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.511. Response shuffle (1,000 runs): rank-1 share = 0.448; best-of-30 share = 0.008.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.09 (p = 0.700), SVR 0.06, LinearRegression 0.06.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated, which fits a drug-specific signal if one exists. Some of the other 29 model settings reach p < 0.05 (14 of 30), but 0.8% of shuffled runs find a best p at least this small, so picking the best setting is not enough on its own. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.09).

#### TCGA rectum (current GDC data) | fluorouracil
- Lab models: 19 (response = ln IC50); proximal pathways: 30; rank-1 Ridge pathway: Purine ribonucleoside monophosphate biosynthesis.
- Drug-treated patients: n = 59, deaths = 5 (2 in predicted responders, 3 in predicted non-responders); log-rank p = 0.409; HR 2.12 (95% CI 0.34-13.17).
- Untreated patients (same split): n = 89, deaths = 18; p = 0.024; HR 3.34 (95% CI 1.09-10.20).
- Other algorithms (rank-1): SVR p = 0.409, LinearRegression p = 0.409. Best of all 30 settings: p = 0.062 (SVR k=6); settings with p < 0.05: 0 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.398. Response shuffle (1,000 runs): rank-1 share = 0.219; best-of-30 share = 0.536.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.58 (p = 0.009), SVR 0.61, LinearRegression 0.53.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. The same split also separates untreated patients, so part of the signal may be general prognosis rather than drug response. Only 5 deaths among treated patients, so any estimate here is fragile.

#### TCGA rectum (current GDC data) | oxaliplatin
- Lab models: 19 (response = ln IC50); proximal pathways: 36; rank-1 Ridge pathway: Metabolism of polyamines.
- Drug-treated patients: n = 42, deaths = 5 (1 in predicted responders, 4 in predicted non-responders); log-rank p = 0.093; HR 5.37 (95% CI 0.60-48.34).
- Untreated patients (same split): n = 89, deaths = 18; p = 0.326; HR 0.62 (95% CI 0.23-1.65).
- Other algorithms (rank-1): SVR p = 0.093, LinearRegression p = 0.133. Best of all 30 settings: p = 0.093 (Ridge k=1); settings with p < 0.05: 0 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.084. Response shuffle (1,000 runs): rank-1 share = 0.126; best-of-30 share = 0.391.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.09 (p = 0.700), SVR 0.06, LinearRegression 0.06.
- What I think: There is a weak trend in drug-treated patients. Untreated patients are not separated, which fits a drug-specific signal if one exists. Only 5 deaths among treated patients, so any estimate here is fragile. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.09).

#### TCGA colon + rectum pooled (current GDC data) | fluorouracil
- Lab models: 19 (response = ln IC50); proximal pathways: 30; rank-1 Ridge pathway: Purine ribonucleoside monophosphate biosynthesis.
- Drug-treated patients: n = 188, deaths = 30 (14 in predicted responders, 16 in predicted non-responders); log-rank p = 0.603; HR 1.21 (95% CI 0.59-2.48).
- Untreated patients (same split): n = 391, deaths = 86; p = 0.104; HR 1.42 (95% CI 0.93-2.19).
- Other algorithms (rank-1): SVR p = 0.603, LinearRegression p = 0.603. Best of all 30 settings: p = 0.343 (SVR k=4); settings with p < 0.05: 0 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.564. Response shuffle (1,000 runs): rank-1 share = 0.775; best-of-30 share = 0.978.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.58 (p = 0.009), SVR 0.61, LinearRegression 0.53.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated, which fits a drug-specific signal if one exists.

#### TCGA colon + rectum pooled (current GDC data) | irinotecan
- Lab models: 19 (response = ln IC50); proximal pathways: 66; rank-1 Ridge pathway: Hyaluronan uptake and degradation.
- Drug-treated patients: n = 39, deaths = 17 (10 in predicted responders, 7 in predicted non-responders); log-rank p = 0.392; HR 0.65 (95% CI 0.25-1.74).
- Untreated patients (same split): n = 391, deaths = 86; p = 0.028; HR 0.61 (95% CI 0.39-0.95).
- Other algorithms (rank-1): SVR p = 0.392, LinearRegression p = 0.392. Best of all 30 settings: p = 0.392 (Ridge k=1); settings with p < 0.05: 0 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.439. Response shuffle (1,000 runs): rank-1 share = 0.703; best-of-30 share = 0.983.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge -0.22 (p = 0.355), SVR -0.21, LinearRegression -0.19.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. The same split also separates untreated patients, so part of the signal may be general prognosis rather than drug response. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = -0.22).

#### TCGA colon + rectum pooled (current GDC data) | oxaliplatin
- Lab models: 19 (response = ln IC50); proximal pathways: 36; rank-1 Ridge pathway: Metabolism of polyamines.
- Drug-treated patients: n = 154, deaths = 26 (14 in predicted responders, 12 in predicted non-responders); log-rank p = 0.695; HR 0.86 (95% CI 0.40-1.86).
- Untreated patients (same split): n = 391, deaths = 86; p = 0.975; HR 0.99 (95% CI 0.65-1.52).
- Other algorithms (rank-1): SVR p = 0.695, LinearRegression p = 0.911. Best of all 30 settings: p = 0.00046 (LinearRegression k=6); settings with p < 0.05: 21 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.802. Response shuffle (1,000 runs): rank-1 share = 0.822; best-of-30 share = 0.007.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.09 (p = 0.700), SVR 0.06, LinearRegression 0.06.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated, which fits a drug-specific signal if one exists. Some of the other 29 model settings reach p < 0.05 (21 of 30), but 0.7% of shuffled runs find a best p at least this small, so picking the best setting is not enough on its own. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.09).

### Pancreatic organoids (Tiriac 2018)  (`PAAD_tiriac/`)

#### TCGA pancreas | fluorouracil
- Lab models: 36 (response = AUC); proximal pathways: 38; rank-1 Ridge pathway: Tetrahydrobiopterin bh4 synthesis recycling salvage and regulation.
- Drug-treated patients: n = 37, deaths = 18 (10 in predicted responders, 8 in predicted non-responders); log-rank p = 0.932; HR 0.96 (95% CI 0.38-2.44).
- Untreated patients (same split): n = 61, deaths = 35; p = 0.281; HR 0.69 (95% CI 0.34-1.37).
- Other algorithms (rank-1): SVR p = 0.420, LinearRegression p = 0.793. Best of all 30 settings: p = 0.277 (Ridge k=4); settings with p < 0.05: 0 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.926. Response shuffle (1,000 runs): rank-1 share = 0.873; best-of-30 share = 0.913.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.01 (p = 0.948), SVR -0.04, LinearRegression -0.19.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated, which fits a drug-specific signal if one exists. Inside the lab models the pathway model predicts drug response weakly (leave-one-out rho = 0.01).

#### TCGA pancreas | gemcitabine
- Lab models: 36 (response = AUC); proximal pathways: 45; rank-1 Ridge pathway: Homologous recombination repair of replication independent double strand breaks.
- Drug-treated patients: n = 104, deaths = 51 (23 in predicted responders, 28 in predicted non-responders); log-rank p = 0.223; HR 1.41 (95% CI 0.81-2.46).
- Untreated patients (same split): n = 61, deaths = 35; p = 0.897; HR 1.05 (95% CI 0.54-2.05).
- Other algorithms (rank-1): SVR p = 0.033, LinearRegression p = 0.151. Best of all 30 settings: p = 0.013 (SVR k=5); settings with p < 0.05: 5 of 30.
- Bootstrap (10,000 random pathway sets of the same size): share with p as small or smaller = 0.296. Response shuffle (1,000 runs): rank-1 share = 0.805; best-of-30 share = 0.388.
- Inside the lab models (leave-one-out Spearman rho, predicted vs observed): Ridge 0.37 (p = 0.025), SVR 0.47, LinearRegression 0.28.
- What I think: With this pathway, survival of the two predicted groups of treated patients is similar. Untreated patients are not separated, which fits a drug-specific signal if one exists. Some of the other 29 model settings reach p < 0.05 (5 of 30), but 38.8% of shuffled runs find a best p at least this small, so picking the best setting is not enough on its own.

