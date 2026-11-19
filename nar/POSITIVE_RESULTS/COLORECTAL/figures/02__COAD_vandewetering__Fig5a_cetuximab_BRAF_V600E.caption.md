# Fig. 5a: cetuximab model vs BRAF V600E in TCGA colon tumours

BRAF V600E is a known marker of cetuximab resistance. A Ridge model was trained on cetuximab response of the lab
models using the cetuximab-proximal pathways; its rank-1 pathway gives each TCGA-COAD tumour a predicted-resistance
score (coefficient x scaled pathway score; higher = more resistant). BRAF V600E status comes from the GDC
open-access somatic mutation calls (cases with whole-exome mutation data). The test is the paper's one-sided
Mann-Whitney U (are V600E tumours predicted more resistant?). Paper: 47 V600E vs 328 wild-type, p = 0.037, using
organoids without the 3 hypermutated lines.

- organoids, authors' NES, 3 hypermutated removed (paper): 16 models; pathway Gastrin creb signalling pathway via pkc and mapk; V600E n = 49, wild-type n = 369; one-sided Mann-Whitney p = 0.031
- organoids, authors' NES, all 19: 19 models; pathway Gastrin creb signalling pathway via pkc and mapk; V600E n = 49, wild-type n = 369; one-sided Mann-Whitney p = 0.031
- organoids, current data, 3 hypermutated removed: 16 models; pathway Gastrin creb signalling pathway via pkc and mapk; V600E n = 50, wild-type n = 377; one-sided Mann-Whitney p = 0.00072
- organoids, current data, all: 19 models; pathway Gastrin creb signalling pathway via pkc and mapk; V600E n = 50, wild-type n = 377; one-sided Mann-Whitney p = 0.00072
- PDX (Isella), current data: 192 models; pathway Signaling by notch3; V600E n = 50, wild-type n = 377; one-sided Mann-Whitney p = 0.266
