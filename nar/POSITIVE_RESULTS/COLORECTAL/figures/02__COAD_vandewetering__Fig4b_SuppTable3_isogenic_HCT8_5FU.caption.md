# Fig. 4b (and Supp. Table 3): isogenic 5-FU-sensitive vs -resistant colorectal cells

HCT-8 colorectal cancer cells (5-FU-sensitive, n = 3) and HCT-8 cells made resistant by growing them in rising
5-FU doses (HCT-8/5-FU at 0, 24 and 48 h, n = 9), GEO GSE81008 (sub-series GSE81005), as in the paper. Each box
shows a pathway's ssGSEA score; dots are single arrays. p = two-sided two-sample Student's t-test (paper's test).
First column: the paper's biomarker (activation of BH3-only proteins). Other columns: the top Ridge pathways of the
model in that row. Row 1 uses the authors' own scores (Source Data sheet 'Figure 4b'); the other rows use scores
recomputed from GEO with this project's ssGSEA. "Direction matches" = the pathway is higher in sensitive cells
when a higher score predicts lower IC50 in organoids (negative Ridge coefficient).

- COADREAD: Activation of bh3 only proteins (Ridge rank 5): sensitive mean 3960, resistant mean 3816, p = 0.039; direction matches the organoid model: True
- COADREAD: Purine ribonucleoside monophosphate biosynthesis (Ridge rank 1): sensitive mean 6565, resistant mean 6195, p = 0.016; direction matches the organoid model: True
- COADREAD: Tetrahydrobiopterin bh4 synthesis recycling salvage and regulation (Ridge rank 2): sensitive mean 4991, resistant mean 4996, p = 0.980; direction matches the organoid model: False
- COADREAD: Hyaluronan uptake and degradation (Ridge rank 3): sensitive mean 3745, resistant mean 3319, p = 0.004; direction matches the organoid model: False
- COAD: Activation of bh3 only proteins (Ridge rank 5): sensitive mean 3960, resistant mean 3816, p = 0.039; direction matches the organoid model: True
- COAD: Purine ribonucleoside monophosphate biosynthesis (Ridge rank 1): sensitive mean 6565, resistant mean 6195, p = 0.016; direction matches the organoid model: True
- COAD: Tetrahydrobiopterin bh4 synthesis recycling salvage and regulation (Ridge rank 2): sensitive mean 4991, resistant mean 4996, p = 0.980; direction matches the organoid model: False
- COAD: Hyaluronan uptake and degradation (Ridge rank 3): sensitive mean 3745, resistant mean 3319, p = 0.004; direction matches the organoid model: False
- COAD_authors: Activation of bh3 only proteins (Ridge rank 1): sensitive mean 0.1922, resistant mean 0.1772, p = 0.008; direction matches the organoid model: True
- COAD_authors: Amino acid synthesis and interconversion transamination (Ridge rank 2): sensitive mean 0.4012, resistant mean 0.3629, p = 1.7e-05; direction matches the organoid model: True
- COAD_authors: Reversible hydration of carbon dioxide (Ridge rank 3): sensitive mean -0.05246, resistant mean -0.08346, p = 0.005; direction matches the organoid model: False
- READ: Activation of bh3 only proteins (Ridge rank 5): sensitive mean 3960, resistant mean 3816, p = 0.039; direction matches the organoid model: True
- READ: Purine ribonucleoside monophosphate biosynthesis (Ridge rank 1): sensitive mean 6565, resistant mean 6195, p = 0.016; direction matches the organoid model: True
- READ: Tetrahydrobiopterin bh4 synthesis recycling salvage and regulation (Ridge rank 2): sensitive mean 4991, resistant mean 4996, p = 0.980; direction matches the organoid model: False
- READ: Hyaluronan uptake and degradation (Ridge rank 3): sensitive mean 3745, resistant mean 3319, p = 0.004; direction matches the organoid model: False
