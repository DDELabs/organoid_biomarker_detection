# Are the mutation checks special to the organoid pathways, and do they survive adjustment?

The paper checked its models against known mutations: bladder tumours with ERCC2 mutations should look more
cisplatin-sensitive, and colon tumours with BRAF V600E should look more cetuximab-resistant. Both reproduced. Here we
ask two questions. (1) Would any pathway do as well? (2) Is the result explained by tumour type features that travel
with the mutation (high tumour mutation burden = TMB for ERCC2; MSI-high status, right-sided colon and TMB for BRAF)?

- **a, b, d, e** Red line = observed p (higher = stronger). Blue = the same test with every Reactome pathway in both
  directions; grey = 2,000 runs of the whole pipeline with random pathway sets of the proximal size (Ridge on the
  organoids, rank-1 pathway, test). The titles give the share of each null that is at least as strong.
- **c, f** Difference in the (rank-normalised) predicted-resistance score between mutant and wild-type tumours, first
  without and then with adjustment (linear regression). For ERCC2 a negative difference is the expected direction;
  for BRAF V600E a positive difference.

- ERCC2, paper data: n = 36 mutant vs 367 wild type, p = 0.003; share of pathways (both directions) as small 0.017; random proximal-sized sets 0.063; TMB-adjusted p = 0.024
- ERCC2, current data: n = 36 mutant vs 366 wild type, p = 0.00055; share of pathways (both directions) as small 0.002; random proximal-sized sets 0.000; TMB-adjusted p = 0.008
- BRAF V600E, authors data (gastrin creb signalling pathway via pkc and mapk): n = 49 vs 369, p = 0.031; pathways 0.272; random sets 0.352; adjusted for MSI, side, stage and TMB p = 0.093; MSS tumours only p = 0.034
- BRAF V600E, current data (gastrin creb signalling pathway via pkc and mapk): n = 50 vs 377, p = 0.00072; pathways 0.163; random sets 0.177; adjusted for MSI, side, stage and TMB p = 0.026; MSS tumours only p = 0.007
