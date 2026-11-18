# What changed between the 2019 colon 5-FU result and today's data

The paper split 5-FU-treated TCGA colon patients by an organoid-trained model and got log-rank p = 0.014. On today's
data the same pipeline gives p = 0.866. Here we change one input at a time to see which change moves the result.
Orange = the paper's method run fresh at each step (the model re-selects its top pathway). Blue = the paper's own
pathway, "Activation of BH3-only proteins", kept fixed (frozen) with the direction learned in 2019.

- **a** Cumulative path from the authors' 2019 inputs (left) to the current inputs (right). Each step adds one
  current input on top of the previous ones. Height = -log10 p (the dashed line is p = 0.05).
- **b** Same steps; the frozen BH3-only score used as a continuous predictor (Cox hazard ratio per standard deviation,
  95% CI). HR > 1 means patients predicted to respond poorly die sooner, the paper's direction.
- **c, d** The reverse view: start from the current setup and put back one of the authors' inputs at a time.

Steps of panel a:
- Step 0 (authors' setup (2019)): n = 114, deaths = 8; re-selected rank-1 = activation of bh3 only proteins, p = 0.014; frozen BH3-only p = 0.014, HR per SD 1.52 (0.94-2.46)
- Step 1 (+ only patients with current STAR data): n = 114, deaths = 8; re-selected rank-1 = activation of bh3 only proteins, p = 0.014; frozen BH3-only p = 0.014, HR per SD 1.52 (0.94-2.46)
- Step 2 (+ organoid scores recomputed): n = 114, deaths = 8; re-selected rank-1 = endosomal sorting complex required for transport escrt, p = 0.396; frozen BH3-only p = 0.014, HR per SD 1.52 (0.94-2.46)
- Step 3 (+ patient scores from current STAR FPKM-UQ): n = 114, deaths = 8; re-selected rank-1 = endosomal sorting complex required for transport escrt, p = 0.839; frozen BH3-only p = 0.319, HR per SD 1.63 (0.84-3.17)
- Step 4 (+ shared gene universe (02 pipeline)): n = 117, deaths = 9; re-selected rank-1 = endosomal sorting complex required for transport escrt, p = 0.228; frozen BH3-only p = 0.596, HR per SD 1.21 (0.60-2.43)
- Step 5 (+ organoid response from curated table): n = 117, deaths = 9; re-selected rank-1 = purine ribonucleoside monophosphate biosynthesis, p = 0.414; frozen BH3-only p = 0.596, HR per SD 1.21 (0.60-2.43)
- Step 6 (+ current drug-record parser): n = 129, deaths = 11; re-selected rank-1 = purine ribonucleoside monophosphate biosynthesis, p = 0.413; frozen BH3-only p = 0.774, HR per SD 1.17 (0.63-2.17)
- Step 7 (+ TCGA-CDR survival = current setup): n = 129, deaths = 25; re-selected rank-1 = purine ribonucleoside monophosphate biosynthesis, p = 0.866; frozen BH3-only p = 0.620, HR per SD 1.20 (0.82-1.75)
