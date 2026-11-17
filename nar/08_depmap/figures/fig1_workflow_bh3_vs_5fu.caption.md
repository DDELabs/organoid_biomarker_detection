**Figure 1. Testing the organoid 5-FU biomarker in colorectal cancer cell lines.**
(a) Workflow. The 5-FU model was trained on 19 colorectal organoids only. Its strongest pathway is "Activation of
BH3-only proteins" (BH3-only proteins are pro-death switches such as BIM, PUMA and NOXA). In organoids, a higher
score went with a LOWER IC50, that is, more sensitivity to 5-FU. We then score the same pathway in cell lines from
four drug screens without changing anything.
(b) Each dot is one colorectal cell line. x = BH3-only activation score from RNA-seq (single-sample rank score).
y = 5-FU response (GDSC: natural log of IC50; CTRPv2 and PRISM: area under the dose-response curve). Higher y means
more resistant. The organoid model predicts a downward slope (negative rho). Spearman rho and p are shown; the
line is a least-squares fit for the eye only. Numbers: `tables/A_b_bh3_vs_5fu_colon.tsv`.
