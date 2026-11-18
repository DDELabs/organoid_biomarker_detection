**Benchmark heatmap: every model (rows) on every task (columns).** Numbers are averages over the 41 lab-model set x drug pairs
that have at least 15 drug-treated TCGA patients (cell-line columns: averages over drug x screen units of GDSC and CTRPv2).
Columns: (1) Spearman correlation between predicted and measured drug response in lab models left out of training
(leave-one-out for sets of up to 40 models, else 10-fold); (2) C-index minus 0.5 in treated TCGA patients (0 = coin toss;
positive = patients predicted resistant died sooner); (3) log hazard ratio per standard deviation of the predicted resistance
(positive = expected direction); (4) cell-blind and (5) drug-blind cross-validation on cell lines (Spearman). Row label colour:
blue = the paper's network method, grey = the paper's comparison methods, red = modern models, green = with the
lab-model-to-tumour alignment of nar/03. "n/a" = the model cannot be run on that task (for example, the paper's method needs
response data for the drug, so it has no drug-blind version). Values are clipped at +/-0.25 for colour only.
