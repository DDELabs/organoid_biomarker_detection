## Cross-set table (paper's cut-off z <= -1.2816; rank-1 Ridge pathway; median split; overall survival)

* = 12-14 treated patients (below the ~15 target; kept and flagged).

| folder | cohort | drug | models | proximal | rank-1 pathway (Ridge) | treated n (deaths: resp vs non-resp) | HR treated (95% CI) | p treated | untreated n (deaths) | HR untreated | p untreated | best of 30 tests | bootstrap share | shuffle share rank-1 | shuffle share best-of-30 | LOO rho Ridge | What I think |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| COAD_vandewetering | COAD_authors | fluorouracil | 19 | 37 | Activation of bh3 only proteins | 114 (8: 0 vs 8) | not estimable | 0.014 | 298 (41) | 0.64 (0.34-1.20) | 0.161 | 0.014 (Ridge k=1) | 0.007 | 0.034 | 0.164 | 0.44 | the split separates treated patients in the expected direction; about 3% of shuffled-response runs did as well; untreated patients are not separated; only 8 deaths, so the estimate is fragile. |
| COAD_vandewetering | COAD_authors | irinotecan | 19 | 76 | Hyaluronan uptake and degradation | 22 (3: 3 vs 0) | not estimable | 0.182 | 298 (41) | 0.89 (0.48-1.67) | 0.718 | 0.083 (Ridge k=3) | 0.251 | 0.317 | 0.495 | -0.21 | treated patients are not separated by this pathway; untreated patients are not separated; only 3 deaths, so the estimate is fragile. |
| COAD_vandewetering | COAD_authors | oxaliplatin | 19 | 43 | Bile acid and bile salt metabolism | 100 (7: 2 vs 5) | 1.18 (0.23-6.11) | 0.844 | 298 (41) | 0.87 (0.47-1.61) | 0.660 | 0.127 (LinearRegression k=5) | 0.844 | 0.946 | 0.569 | 0.15 | treated patients are not separated by this pathway; untreated patients are not separated; only 7 deaths, so the estimate is fragile. |
| COAD_vandewetering | COAD | cetuximab | 19 | 32 | Gastrin creb signalling pathway via pkc and mapk | 15 (7: 4 vs 3) | 0.44 (0.08-2.42) | 0.329 | 302 (68) | 1.37 (0.85-2.21) | 0.193 | 0.018 (Ridge k=10) | 0.403 | 0.539 | 0.311 | -0.40 | treated patients are not separated by this pathway; untreated patients are not separated; only 7 deaths, so the estimate is fragile. |
| COAD_vandewetering | COAD | fluorouracil | 19 | 30 | Purine ribonucleoside monophosphate biosynthesis | 129 (25: 12 vs 13) | 1.07 (0.49-2.35) | 0.866 | 302 (68) | 1.22 (0.75-1.96) | 0.419 | 0.179 (SVR k=4) | 0.885 | 0.941 | 0.790 | 0.58 | treated patients are not separated by this pathway; untreated patients are not separated. |
| COAD_vandewetering | COAD | irinotecan | 19 | 66 | Hyaluronan uptake and degradation | 32 (15: 9 vs 6) | 0.55 (0.19-1.56) | 0.257 | 302 (68) | 0.60 (0.36-0.99) | 0.043 | 0.257 (Ridge k=1) | 0.199 | 0.238 | 0.787 | -0.22 | treated patients are not separated by this pathway; untreated patients are also separated, so part of it may be general prognosis. |
| COAD_vandewetering | COAD | oxaliplatin | 19 | 36 | Metabolism of polyamines | 112 (21: 12 vs 9) | 0.71 (0.30-1.70) | 0.440 | 302 (68) | 1.19 (0.74-1.93) | 0.470 | 0.001 (Ridge k=4) | 0.511 | 0.448 | 0.008 | 0.09 | treated patients are not separated by this pathway; untreated patients are not separated. |
| COAD_vandewetering | READ | fluorouracil | 19 | 30 | Purine ribonucleoside monophosphate biosynthesis | 59 (5: 2 vs 3) | 2.12 (0.34-13.17) | 0.409 | 89 (18) | 3.34 (1.09-10.20) | 0.024 | 0.062 (SVR k=6) | 0.398 | 0.219 | 0.536 | 0.58 | treated patients are not separated by this pathway; untreated patients are also separated, so part of it may be general prognosis; only 5 deaths, so the estimate is fragile. |
| COAD_vandewetering | READ | oxaliplatin | 19 | 36 | Metabolism of polyamines | 42 (5: 1 vs 4) | 5.37 (0.60-48.34) | 0.093 | 89 (18) | 0.62 (0.23-1.65) | 0.326 | 0.093 (Ridge k=1) | 0.084 | 0.126 | 0.391 | 0.09 | a weak trend in treated patients; untreated patients are not separated; only 5 deaths, so the estimate is fragile. |
| COAD_vandewetering | COADREAD | fluorouracil | 19 | 30 | Purine ribonucleoside monophosphate biosynthesis | 188 (30: 14 vs 16) | 1.21 (0.59-2.48) | 0.603 | 391 (86) | 1.42 (0.93-2.19) | 0.104 | 0.343 (SVR k=4) | 0.564 | 0.775 | 0.978 | 0.58 | treated patients are not separated by this pathway; untreated patients are not separated. |
| COAD_vandewetering | COADREAD | irinotecan | 19 | 66 | Hyaluronan uptake and degradation | 39 (17: 10 vs 7) | 0.65 (0.25-1.74) | 0.392 | 391 (86) | 0.61 (0.39-0.95) | 0.028 | 0.392 (Ridge k=1) | 0.439 | 0.703 | 0.983 | -0.22 | treated patients are not separated by this pathway; untreated patients are also separated, so part of it may be general prognosis. |
| COAD_vandewetering | COADREAD | oxaliplatin | 19 | 36 | Metabolism of polyamines | 154 (26: 14 vs 12) | 0.86 (0.40-1.86) | 0.695 | 391 (86) | 0.99 (0.65-1.52) | 0.975 | 0.00046 (LinearRegression k=6) | 0.802 | 0.822 | 0.007 | 0.09 | treated patients are not separated by this pathway; untreated patients are not separated. |
| PAAD_tiriac | PAAD | fluorouracil | 36 | 38 | Tetrahydrobiopterin bh4 synthesis recycling salva… | 37 (18: 10 vs 8) | 0.96 (0.38-2.44) | 0.932 | 61 (35) | 0.69 (0.34-1.37) | 0.281 | 0.277 (Ridge k=4) | 0.926 | 0.873 | 0.913 | 0.01 | treated patients are not separated by this pathway; untreated patients are not separated. |
| PAAD_tiriac | PAAD | gemcitabine | 36 | 45 | Homologous recombination repair of replication in… | 104 (51: 23 vs 28) | 1.41 (0.81-2.46) | 0.223 | 61 (35) | 1.05 (0.54-2.05) | 0.897 | 0.013 (SVR k=5) | 0.296 | 0.805 | 0.388 | 0.37 | treated patients are not separated by this pathway; untreated patients are not separated. |
| PAAD_tiriac | PAAD | oxaliplatin | 36 | 47 | Synthesis of substrates in n glycan biosythesis | 21 (11: 8 vs 3) | 0.44 (0.11-1.71) | 0.223 | 61 (35) | 0.53 (0.26-1.07) | 0.072 | 0.201 (LinearRegression k=10) | 0.160 | 0.135 | 0.645 | 0.14 | treated patients are not separated by this pathway; untreated patients are not separated; only 11 deaths, so the estimate is fragile. |
| PAAD_tiriac | PAAD | sn-38 * | 36 | 78 | Regulation of the fanconi anemia pathway | 13 (5: 3 vs 2) | 0.59 (0.10-3.63) | 0.570 | 61 (35) | 1.58 (0.80-3.11) | 0.186 | 0.023 (Ridge k=9) | 0.610 | 0.626 | 0.189 | 0.02 | treated patients are not separated by this pathway; untreated patients are not separated; only 5 deaths, so the estimate is fragile. |
| PAAD_shi | PAAD | fluorouracil | 39 | 38 | Inhibition of replication initiation of damaged d… | 37 (18: 10 vs 8) | 0.98 (0.38-2.48) | 0.958 | 61 (35) | 0.75 (0.39-1.48) | 0.412 | 0.097 (Ridge k=3) | 0.947 | 0.945 | 0.533 | -0.01 | treated patients are not separated by this pathway; untreated patients are not separated. |
| PAAD_shi | PAAD | gemcitabine | 39 | 45 | G1 s specific transcription | 104 (51: 32 vs 19) | 0.52 (0.30-0.92) | 0.023 | 61 (35) | 0.45 (0.23-0.88) | 0.017 | 0.005 (Ridge k=8) | 0.022 | 0.100 | 0.146 | -0.08 | the split separates treated patients, but predicted responders did worse (opposite direction); about 10% of shuffled-response runs did as well; untreated patients are also separated, so part of it may be general prognosis. |
| PAAD_shi | PAAD | irinotecan * | 39 | 78 | Fanconi anemia pathway | 13 (5: 4 vs 1) | 0.33 (0.04-3.00) | 0.302 | 61 (35) | 2.22 (1.11-4.43) | 0.021 | 0.030 (SVR k=2) | 0.310 | 0.351 | 0.305 | -0.27 | treated patients are not separated by this pathway; untreated patients are also separated, so part of it may be general prognosis; only 5 deaths, so the estimate is fragile. |
| PAAD_shi | PAAD | oxaliplatin | 39 | 46 | Intrinsic pathway | 21 (11: 7 vs 4) | 0.65 (0.18-2.30) | 0.497 | 61 (35) | 1.49 (0.76-2.91) | 0.242 | 0.009 (LinearRegression k=10) | 0.486 | 0.277 | 0.027 | -0.35 | treated patients are not separated by this pathway; untreated patients are not separated; only 11 deaths, so the estimate is fragile. |
| LIHC_ji | LIHC | sorafenib | 61 | 22 | Signaling by activated point mutants of fgfr1 | 30 (15: 8 vs 7) | 0.48 (0.15-1.47) | 0.188 | 328 (108) | 1.34 (0.92-1.96) | 0.129 | 0.013 (LinearRegression k=8) | 0.397 | 0.583 | 0.401 | -0.16 | treated patients are not separated by this pathway; untreated patients are not separated. |
| LIHC_broutier | LIHC | sorafenib | 5 | 15 | Reversible hydration of carbon dioxide | 30 (15: 8 vs 7) | 0.43 (0.14-1.35) | 0.139 | 328 (108) | 1.19 (0.81-1.73) | 0.374 | 0.093 (SVR k=4) | 0.364 | 0.521 | 0.828 | -0.70 | treated patients are not separated by this pathway; untreated patients are not separated. |
| COAD_isella_pdx | COAD | cetuximab | 192 | 41 | Signaling by notch3 | 15 (7: 4 vs 3) | 0.70 (0.15-3.22) | 0.648 | 302 (68) | 1.63 (1.01-2.65) | 0.046 | 0.370 (LinearRegression k=7) | 0.831 | 0.760 | 1.000 | 0.50 | treated patients are not separated by this pathway; untreated patients are also separated, so part of it may be general prognosis; only 7 deaths, so the estimate is fragile. |

## Cut-off sensitivity (rank-1 Ridge, treated p; number of proximal pathways in brackets)

| folder | cohort | drug | -1.1824 | -1.2816 | -1.6449 | -2.19254 |
|---|---|---|---|---|---|---|
| COAD_isella_pdx | COAD | cetuximab | 0.648 (49) | 0.648 (41) | 0.648 (24) | 0.020 (13) |
| COAD_vandewetering | COADREAD | fluorouracil | 0.603 (32) | 0.603 (30) | 0.470 (19) | 0.470 (13) |
| COAD_vandewetering | COADREAD | irinotecan | 0.392 (73) | 0.392 (66) | 0.392 (44) | 0.168 (10) |
| COAD_vandewetering | COADREAD | oxaliplatin | 0.911 (40) | 0.695 (36) | 0.817 (24) | 0.457 (10) |
| COAD_vandewetering | COAD | cetuximab | 0.928 (41) | 0.329 (32) | 0.329 (20) | 0.003 (11) |
| COAD_vandewetering | COAD | fluorouracil | 0.866 (32) | 0.866 (30) | 0.255 (19) | 0.255 (13) |
| COAD_vandewetering | COAD | irinotecan | 0.257 (73) | 0.257 (66) | 0.257 (44) | 0.760 (10) |
| COAD_vandewetering | COAD | oxaliplatin | 0.880 (40) | 0.440 (36) | 0.346 (24) | 0.165 (10) |
| COAD_vandewetering | COAD_authors | fluorouracil | 0.014 (40) | 0.014 (37) | 0.303 (23) | 0.303 (14) |
| COAD_vandewetering | COAD_authors | irinotecan | 0.182 (83) | 0.182 (76) | 0.763 (49) | 0.564 (13) |
| COAD_vandewetering | COAD_authors | oxaliplatin | 0.761 (48) | 0.844 (43) | 0.844 (30) | 0.844 (14) |
| COAD_vandewetering | READ | fluorouracil | 0.409 (32) | 0.409 (30) | 0.802 (19) | 0.802 (13) |
| COAD_vandewetering | READ | oxaliplatin | 0.133 (40) | 0.093 (36) | 0.243 (24) | 0.202 (10) |
| LIHC_broutier | LIHC | sorafenib | 0.879 (20) | 0.139 (15) | 0.139 (6) | 0 pathways |
| LIHC_ji | LIHC | sorafenib | 0.188 (27) | 0.188 (22) | 0.778 (12) | 0 pathways |
| PAAD_shi | PAAD | fluorouracil | 0.563 (42) | 0.958 (38) | 0.784 (23) | 0.860 (14) |
| PAAD_shi | PAAD | gemcitabine | 0.461 (51) | 0.023 (45) | 0.023 (31) | 0.023 (14) |
| PAAD_shi | PAAD | irinotecan | 0.302 (85) | 0.302 (78) | 0.302 (51) | 0.959 (15) |
| PAAD_shi | PAAD | oxaliplatin | 0.569 (51) | 0.497 (46) | 0.226 (33) | 0.636 (16) |
| PAAD_tiriac | PAAD | fluorouracil | 0.932 (42) | 0.932 (38) | 0.793 (23) | 0.793 (14) |
| PAAD_tiriac | PAAD | gemcitabine | 0.043 (51) | 0.223 (45) | 0.151 (31) | 0.074 (14) |
| PAAD_tiriac | PAAD | oxaliplatin | 0.506 (53) | 0.223 (47) | 0.483 (33) | 0.562 (16) |
| PAAD_tiriac | PAAD | sn-38 | 0.570 (85) | 0.570 (78) | 0.570 (51) | 0.302 (15) |
