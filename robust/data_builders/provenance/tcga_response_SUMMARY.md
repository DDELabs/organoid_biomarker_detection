# TCGA pan-cancer drug-response labels

Built by `robust/obd/tcga_response.py::build()`; load with `load_response_labels()`, `load_cdr()`,
`drug_label_matrix()`.

* Drug records (rows, one per patient x drug x regimen record): **12,423**
  (4,230 patients, 32 projects, 442 distinct drug names,
  97.5% of rows resolved to a DrugBank/curated name)
* Rows with a RECIST-like label (CR/PR/SD/PD): **3,832** (1,650 patients);
  responders (CR/PR) 2,269 (59.2%);
  setting_proxy of labelled rows: 62.2% early (started <=180 d after diagnosis, no prior
  progression; adjuvant/neoadjuvant/first-line), 28.2% post_progression (started after the CDR PFI event).
* `recist` harmonised from `treatment_best_response`: Complete Response->CR, Partial Response->PR,
  Stable Disease->SD, Clinical/Radiographic Progressive Disease->PD; `responder` = 1 for CR/PR, 0 for SD/PD.
* `n_drugs_in_regimen`: distinct drugs with the same patient and start day (same biotab record if start missing).
* `regimen_indication`: `regimen_indication` or `therapy_regimen` column (ADJUVANT, PROGRESSION, RECURRENCE,
  PALLIATIVE, PRIMARY, OTHER...); absent for 43% of rows and, in these
  biotab versions, NEVER present on a row that carries a response (the form versions with `therapy_regimen`
  dropped `treatment_best_response`). Use `setting_proxy` instead.
* `setting_proxy`: post_progression (start >= CDR PFI event time - 14 d, or indication PROGRESSION/RECURRENCE),
  early (start <= 180 d from diagnosis, or indication ADJUVANT), late_no_event, unknown (no start day).
* `in_ding2016`: (patient, drug) also present in the Ding, Zu & Gu 2016 curated table (2,484 pairs);
  RECIST agreement on 2,425 shared labelled pairs: 99.3%.
* LAML: TCGA released no clinical_drug biotab table (no labels). DLBC's table has no response field.
* survival_cdr.tsv.gz: TCGA-CDR (Liu et al. Cell 2018) for 11,081 patients.

## Top 40 drugs by patients with a RECIST label

Patient-level label = earliest labelled record of the drug for that patient. `post_progression` restricts to
records with setting_proxy == post_progression (treatment of measurable recurrent/progressive disease).

| drug | n_patients | n_responders | responder_rate | n_patients_post_progression | responder_rate_post_progression | n_projects | projects |
|---|---|---|---|---|---|---|---|
| CISPLATIN | 433 | 335 | 0.77 | 69 | 0.57 | 21 | CESC:67,BLCA:65,LUAD:52,TGCT:50,HNSC:50,STAD:39,LUSC:32,MESO:21,ESCA:17,UCEC:9,SKCM:7,UCS:6,THYM:4,SARC:3,CHOL:3,OV:2,LIHC:2,PAAD:1,ACC:1,BRCA:1,KIRP:1 |
| FLUOROURACIL | 285 | 204 | 0.72 | 43 | 0.26 | 10 | STAD:93,COAD:65,BRCA:51,READ:34,PAAD:20,ESCA:11,CESC:6,HNSC:3,BLCA:1,LIHC:1 |
| CARBOPLATIN | 263 | 172 | 0.65 | 64 | 0.28 | 20 | UCEC:54,LUAD:43,UCS:29,HNSC:28,LUSC:28,BLCA:19,CESC:11,BRCA:9,MESO:9,OV:7,ESCA:6,SKCM:5,TGCT:3,STAD:3,LGG:3,KIRP:2,SARC:1,ACC:1,PAAD:1,PRAD:1 |
| PACLITAXEL | 250 | 177 | 0.71 | 39 | 0.31 | 14 | BRCA:64,UCEC:53,UCS:30,LUAD:27,HNSC:22,CESC:12,BLCA:10,LUSC:10,OV:6,ESCA:6,STAD:4,SKCM:3,PAAD:2,TGCT:1 |
| GEMCITABINE | 238 | 110 | 0.46 | 85 | 0.15 | 18 | BLCA:74,PAAD:70,SARC:34,LUSC:17,LUAD:9,CHOL:8,MESO:6,LIHC:4,BRCA:4,ESCA:2,CESC:2,SKCM:2,KIRP:1,UCS:1,OV:1,PCPG:1,HNSC:1,UCEC:1 |
| CYCLOPHOSPHAMIDE | 174 | 160 | 0.92 | 7 | 0.43 | 11 | BRCA:162,PCPG:3,PAAD:1,LIHC:1,CESC:1,UCS:1,STAD:1,SKCM:1,KIRP:1,LGG:1,THYM:1 |
| DOXORUBICIN | 171 | 125 | 0.73 | 42 | 0.29 | 15 | BRCA:97,SARC:31,BLCA:14,UCEC:12,LIHC:2,STAD:2,CESC:2,UCS:2,THCA:2,MESO:2,PAAD:1,ACC:1,OV:1,KIRP:1,THYM:1 |
| DOCETAXEL | 161 | 109 | 0.68 | 46 | 0.13 | 13 | BRCA:73,SARC:31,LUSC:13,LUAD:12,HNSC:9,UCEC:4,STAD:4,BLCA:4,UCS:3,PRAD:3,OV:2,ESCA:2,PAAD:1 |
| TEMOZOLOMIDE | 148 | 25 | 0.17 | 40 | 0.05 | 6 | LGG:131,GBM:12,LIHC:2,PCPG:1,LUSC:1,SKCM:1 |
| LEUCOVORIN | 132 | 90 | 0.68 | 28 | 0.25 | 5 | COAD:62,STAD:36,READ:25,PAAD:6,ESCA:3 |
| OXALIPLATIN | 112 | 70 | 0.62 | 27 | 0.26 | 8 | COAD:52,STAD:23,READ:20,PAAD:10,ESCA:4,BLCA:1,LUSC:1,LIHC:1 |
| ETOPOSIDE | 108 | 87 | 0.81 | 33 | 0.61 | 11 | TGCT:52,STAD:17,LUAD:15,LUSC:8,LGG:7,BLCA:3,THYM:2,ACC:1,SARC:1,READ:1,ESCA:1 |
| CAPECITABINE | 95 | 61 | 0.64 | 22 | 0.27 | 9 | STAD:33,ESCA:22,COAD:19,PAAD:7,BRCA:7,READ:3,HNSC:2,CHOL:1,PCPG:1 |
| BEVACIZUMAB | 76 | 20 | 0.26 | 48 | 0.17 | 15 | LGG:24,COAD:20,LUAD:9,READ:5,GBM:4,BRCA:4,SARC:2,BLCA:1,SKCM:1,CESC:1,MESO:1,UCS:1,ACC:1,HNSC:1,KIRC:1 |
| PEMETREXED | 58 | 27 | 0.47 | 20 | 0.35 | 7 | LUAD:29,MESO:24,KIRP:1,HNSC:1,LUSC:1,SARC:1,THYM:1 |
| BLEOMYCIN | 54 | 50 | 0.93 | 16 | 0.94 | 3 | TGCT:51,CESC:2,SKCM:1 |
| EPIRUBICIN | 52 | 44 | 0.85 | 4 | 0.25 | 3 | STAD:26,BRCA:25,ESCA:1 |
| VINORELBINE | 42 | 31 | 0.74 | 9 | 0.22 | 4 | LUSC:21,LUAD:16,BRCA:4,MESO:1 |
| IRINOTECAN | 40 | 11 | 0.28 | 23 | 0.22 | 6 | COAD:19,LGG:8,PAAD:5,READ:4,STAD:3,LUAD:1 |
| DACARBAZINE | 38 | 14 | 0.37 | 25 | 0.28 | 4 | SKCM:27,SARC:7,PCPG:3,STAD:1 |
| IFOSFAMIDE | 33 | 18 | 0.55 | 19 | 0.37 | 4 | SARC:16,UCS:8,BLCA:7,TGCT:2 |
| TAMOXIFEN | 32 | 20 | 0.62 | 7 | 0.0 | 4 | BRCA:25,SKCM:4,SARC:2,LGG:1 |
| CETUXIMAB | 27 | 13 | 0.48 | 10 | 0.2 | 3 | HNSC:20,COAD:6,STAD:1 |
| LOMUSTINE | 27 | 5 | 0.19 | 15 | 0.0 | 3 | LGG:23,SKCM:2,GBM:2 |
| LEUPROLIDE | 24 | 18 | 0.75 | 8 | 0.38 | 3 | PRAD:21,BRCA:2,SKCM:1 |
| SORAFENIB | 23 | 4 | 0.17 | 20 | 0.1 | 5 | LIHC:18,KICH:2,KIRC:1,SARC:1,ACC:1 |
| ANASTROZOLE | 22 | 17 | 0.77 | 2 | 0.0 | 1 | BRCA:22 |
| METHOTREXATE | 20 | 13 | 0.65 | 4 | 0.0 | 4 | BLCA:10,BRCA:8,STAD:1,HNSC:1 |
| VINCRISTINE | 19 | 9 | 0.47 | 8 | 0.38 | 7 | LGG:10,PCPG:3,SKCM:2,BRCA:1,STAD:1,CESC:1,KIRP:1 |
| INTERFERON ALFA-2B, RECOMBINANT | 19 | 11 | 0.58 | 4 | 0.5 | 2 | SKCM:18,KICH:1 |
| BICALUTAMIDE | 18 | 15 | 0.83 | 4 | 0.5 | 1 | PRAD:18 |
| TRASTUZUMAB | 18 | 16 | 0.89 | 2 | 0.0 | 2 | BRCA:17,SKCM:1 |
| VINBLASTINE | 17 | 11 | 0.65 | 3 | 0.33 | 4 | BLCA:10,SKCM:5,LUAD:1,BRCA:1 |
| ERLOTINIB | 15 | 2 | 0.13 | 11 | 0.0 | 6 | LUAD:7,LUSC:4,BLCA:1,PAAD:1,MESO:1,LGG:1 |
| IPILIMUMAB | 15 | 5 | 0.33 | 11 | 0.36 | 1 | SKCM:15 |
| PROCARBAZINE | 13 | 3 | 0.23 | 6 | 0.0 | 2 | LGG:12,GBM:1 |
| CARMUSTINE | 11 | 0 | 0.0 | 9 | 0.0 | 2 | LGG:10,SKCM:1 |
| SUNITINIB | 11 | 0 | 0.0 | 5 | 0.0 | 6 | KIRC:4,KIRP:3,LIHC:1,ACC:1,MESO:1,KICH:1 |
| LETROZOLE | 10 | 6 | 0.6 | 4 | 0.0 | 2 | BRCA:9,SARC:1 |
| PAZOPANIB | 8 | 1 | 0.12 | 6 | 0.0 | 3 | SARC:4,KIRP:2,KIRC:2 |

## Per project

| project | drug_records | patients_treated | labelled_records | labelled_patients | responder_records | cdr_patients |
|---|---|---|---|---|---|---|
| ACC | 18 | 5 | 13 | 3 | 0 | 92 |
| BLCA | 295 | 113 | 247 | 97 | 133 | 412 |
| BRCA | 2404 | 772 | 635 | 214 | 544 | 1097 |
| CESC | 198 | 140 | 114 | 80 | 83 | 307 |
| CHOL | 14 | 9 | 12 | 8 | 5 | 36 |
| COAD | 609 | 153 | 272 | 82 | 160 | 458 |
| DLBC | 252 | 43 | 0 | 0 | 0 | 48 |
| ESCA | 85 | 42 | 76 | 40 | 55 | 185 |
| GBM | 1464 | 443 | 26 | 14 | 2 | 588 |
| HNSC | 335 | 169 | 149 | 88 | 102 | 528 |
| KICH | 21 | 10 | 5 | 4 | 2 | 66 |
| KIRC | 167 | 84 | 16 | 12 | 1 | 536 |
| KIRP | 47 | 23 | 17 | 10 | 7 | 291 |
| LAML | 0 | 0 | 0 | 0 | 0 | 200 |
| LGG | 707 | 285 | 314 | 142 | 38 | 515 |
| LIHC | 63 | 37 | 40 | 24 | 9 | 377 |
| LUAD | 445 | 169 | 241 | 101 | 144 | 519 |
| LUSC | 350 | 135 | 143 | 63 | 94 | 504 |
| MESO | 151 | 63 | 93 | 35 | 29 | 87 |
| OV | 2448 | 527 | 21 | 9 | 17 | 584 |
| PAAD | 247 | 116 | 144 | 79 | 50 | 185 |
| PCPG | 12 | 4 | 12 | 4 | 3 | 179 |
| PRAD | 130 | 72 | 61 | 40 | 45 | 498 |
| READ | 231 | 76 | 106 | 41 | 67 | 166 |
| SARC | 205 | 70 | 165 | 59 | 54 | 261 |
| SKCM | 231 | 129 | 156 | 85 | 71 | 470 |
| STAD | 466 | 188 | 339 | 140 | 208 | 443 |
| TGCT | 171 | 61 | 163 | 59 | 159 | 134 |
| THCA | 56 | 43 | 12 | 12 | 11 | 507 |
| THYM | 16 | 8 | 9 | 4 | 4 | 124 |
| UCEC | 486 | 198 | 144 | 65 | 114 | 547 |
| UCS | 95 | 39 | 87 | 36 | 58 | 57 |
| UVM | 4 | 4 | 0 | 0 | 0 | 80 |

## SOURCE

Mirrors of the TCGA BCR biotab `nationwidechildrens.org_clinical_drug_<project>.txt` (GDC legacy archive,
data frozen ~2016; GDC API blocked from the build environment). kemplab/FBA-pipeline covers 31 projects,
nikcheerla/mirnanalyze provides ACC. Shicheng-Guo/HowtoBook `TCGA/drug_response/mutation/tcga.TCGA-*.drug.txt`
(TCGAbiolinks-prepared) was checked as a cross-mirror (same records, ~1% fewer) but not used.

| file | url | sha256 |
|---|---|---|
| raw/nationwidechildrens.org_clinical_drug_acc.txt | https://raw.githubusercontent.com/nikcheerla/mirnanalyze/master/Scripts/cancer_clinical_data/acc/nationwidechildrens.org_clinical_drug_acc.txt | `ea274f34a4cc16b15097bf675082272ad784458f8036c2a0c11bb43596d8ccb5` |
| raw/nationwidechildrens.org_clinical_drug_blca.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_blca.txt | `03ee4d30349944837dd4da6dc2bc987b816adf28630eb7d9d7efe2459d632f11` |
| raw/nationwidechildrens.org_clinical_drug_brca.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_brca.txt | `10d58859f43895eeb6ba5f14ec3b3ea3e670f9174a99fdce45588972836d804e` |
| raw/Survival_SupplementalTable_S1_20171025_xena_sp.tsv | https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp | `a5e704158bb5c51cded8a368accd999dfb259d428e88bbb0c4386078c5df9617` |
| raw/nationwidechildrens.org_clinical_drug_cesc.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_cesc.txt | `0ab51fa0b2820d7348c9497e5bd2a2cb5420a6c05448d662222b24efb929e24f` |
| raw/nationwidechildrens.org_clinical_drug_chol.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_chol.txt | `08de11866a463e1c7ebfebf729d1eee2b50c7bbd4c37bb0924d0cedb3cca7953` |
| raw/nationwidechildrens.org_clinical_drug_coad.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_coad.txt | `4c063a87f615625d0768e1536c70051d828c6ac785a62403972887e460608f3a` |
| raw/ding2016_vaen_drug_response.txt | https://raw.githubusercontent.com/bsml320/VAEN/master/DATA/response/drug_response.txt | `fc564687c4e0af55ce84e1f03e69ed470cca976a606f6e4fe6f5b00a877d9382` |
| raw/nationwidechildrens.org_clinical_drug_dlbc.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_dlbc.txt | `c2aaf29f197ce593f6f899fd934800f234457b2a23efb8a5dd1909aa2c3a937d` |
| raw/nationwidechildrens.org_clinical_drug_esca.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_esca.txt | `ad88d2f9e2bc3ed1fe5855f51c11d60faa202fa0a9245ed32e90c33a5809652a` |
| raw/nationwidechildrens.org_clinical_drug_gbm.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_gbm.txt | `524c183ab1c49fd551c5841a42d1f33659de8f1eb88453caddce4e92f35a4065` |
| raw/nationwidechildrens.org_clinical_drug_hnsc.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_hnsc.txt | `6662f2270037ef9bd85146cc1c27920fd7b750c2f6fa323016bf0587db0bed26` |
| raw/nationwidechildrens.org_clinical_drug_kich.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_kich.txt | `8cac65f349522fe7bd01961bf8bd9f057d8e75d46966b00a0dcd0ead0c416437` |
| raw/nationwidechildrens.org_clinical_drug_kirc.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_kirc.txt | `93448dfd47abeda804429eb30db0dcd276df77bbde90ceb5f8392dcfa446ad07` |
| raw/nationwidechildrens.org_clinical_drug_kirp.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_kirp.txt | `f60498d19a26d3f29326e05e3c1c1d7a6a64e492353d5992068db6bebf77a823` |
| raw/nationwidechildrens.org_clinical_drug_lgg.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_lgg.txt | `bca8d620aa5af0621042d5e7f99f277528b5a71948267b32e2342aa482c17bcd` |
| raw/nationwidechildrens.org_clinical_drug_lihc.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_lihc.txt | `7e15d3321d3913158e945cc87fb9c3ebd77ec90ce406733d40a9e862cebfe100` |
| raw/nationwidechildrens.org_clinical_drug_luad.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_luad.txt | `9f6977eae24d77017b8f3bec1d7564349e75334827dc67d8bef9e88820e640f2` |
| raw/nationwidechildrens.org_clinical_drug_lusc.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_lusc.txt | `48e33101472eeebad8ec59269a76e68f172959eed8d66144e86dbd37d6b69b45` |
| raw/nationwidechildrens.org_clinical_drug_meso.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_meso.txt | `28f8c5c12768cbd111c81c4965eca933bcff16bb46f8905211ec142e61c9c9ec` |
| raw/nationwidechildrens.org_clinical_drug_ov.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_ov.txt | `c5c53d0b6dd6176c22747b3851c206ad80c39af3c008aea5dba799d97e0eecb3` |
| raw/nationwidechildrens.org_clinical_drug_paad.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_paad.txt | `a8cfff2d4f9d59da2f79b3d8c66d69d356dfe89a1a550e574d8459786a657c82` |
| raw/nationwidechildrens.org_clinical_drug_pcpg.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_pcpg.txt | `9bf6d6f47dac464eb6fd006051a6c2aedc6ddc43c6a87d95611b65d92698fbf9` |
| raw/nationwidechildrens.org_clinical_drug_prad.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_prad.txt | `4c3d417d2ed8086032727e6862a075dc7b7a38ed4635f3fae1298a1fe84cc180` |
| raw/nationwidechildrens.org_clinical_drug_read.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_read.txt | `33e3be6f805a73df5b0d17f0082598d9e5e369911d99fd5a4b835b26991f5292` |
| raw/nationwidechildrens.org_clinical_drug_sarc.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_sarc.txt | `70541a110d72120174a704d7bf4412c471d2aa9f44840a29f78909418379880f` |
| raw/nationwidechildrens.org_clinical_drug_skcm.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_skcm.txt | `5392c85e8a5021ea5afa8c8c3cc7ac8c6548a6cd94a1e85fa4c5f02d2b138019` |
| raw/nationwidechildrens.org_clinical_drug_stad.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_stad.txt | `32d788bc15b596b1c39ad256cea3f96d890804e790fe98167ff1542dca4e0dfb` |
| raw/nationwidechildrens.org_clinical_drug_tgct.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_tgct.txt | `ffcaa995613732c5fe4d746c2755b7943086a5ecd5b12883c7e9dd919bf15bf6` |
| raw/nationwidechildrens.org_clinical_drug_thca.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_thca.txt | `c02263133d93e99614078bb2062f459f967e35e5be4ea245330ef5e0cf7b73ea` |
| raw/nationwidechildrens.org_clinical_drug_thym.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_thym.txt | `5cbade746d1e029907cf430f15542d9ca91776742d8960a3f267830a19558c8f` |
| raw/nationwidechildrens.org_clinical_drug_ucec.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_ucec.txt | `810a80878a0a39a5a0a2af128f6f5344c9d236a4b2fa9312238e4c8b069507ab` |
| raw/nationwidechildrens.org_clinical_drug_ucs.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_ucs.txt | `b171e529e3f4dcb916f0c020ea38568a9cc465bb175dcc275578d615943d4e8a` |
| raw/nationwidechildrens.org_clinical_drug_uvm.txt | https://raw.githubusercontent.com/kemplab/FBA-pipeline/master/Code%20%2B%20Models/data/clinical/_data_/input/TCGA/nationwidechildrens.org_clinical_drug_uvm.txt | `718abe00694cb6d72440c3e47e81ae74a76ee0394d8ff81714f2acd5bdaec50d` |
| drug_response_long.tsv.gz | derived (obd.tcga_response.build) | `09ec096a63a15a51bb51eca25ca9cd840e9a778b9eaa230ee8e784aae6099402` |
| survival_cdr.tsv.gz | derived (obd.tcga_response.build) | `6f31a9d29d7e912a3bed19b55c6638851b2aef66b78a81952dfe1e987483b5cb` |
