# 5-FU-treated colorectal patients split by regimen

The paper's "5-FU-treated" patients are everyone with a 5-FU record. Most of them also received oxaliplatin
(FOLFOX). Here the same predicted groups (rank-1 Ridge pathway, median split over all 5-FU-treated patients) are
shown separately for patients who also had oxaliplatin, patients who had neither oxaliplatin nor irinotecan
(closest to "5-FU alone", often with leucovorin or capecitabine), and patients with irinotecan but no oxaliplatin.

- COADREAD, with oxaliplatin (FOLFOX-type): n = 146, deaths = 22, p = 0.478, HR 1.35 (95% CI 0.58-3.14)
- COADREAD, no oxaliplatin or irinotecan record: n = 31, deaths = 1, p = 0.405, HR not estimable (0 deaths in one group)
- COADREAD, with irinotecan, no oxaliplatin (FOLFIRI-type): n = 11, deaths = 7, p = 0.980, HR 0.98 (95% CI 0.19-5.03)
- COAD, with oxaliplatin (FOLFOX-type): n = 106, deaths = 18, p = 0.583, HR 1.30 (95% CI 0.51-3.29)
- COAD, no oxaliplatin or irinotecan record: n = 12, deaths = 0, p = 1.000, HR not estimable (0 deaths in one group)
- COAD, with irinotecan, no oxaliplatin (FOLFIRI-type): n = 11, deaths = 7, p = 0.980, HR 0.98 (95% CI 0.19-5.03)
- COAD_authors, with oxaliplatin (FOLFOX-type): n = 94, deaths = 6, p = 0.036, HR not estimable (0 deaths in one group)
- COAD_authors, no oxaliplatin or irinotecan record: n = 11, deaths = 0, p = 1.000, HR not estimable (0 deaths in one group)
- COAD_authors, with irinotecan, no oxaliplatin (FOLFIRI-type): n = 9, deaths = 2, p = NA, HR not estimable (0 deaths in one group)
- READ, with oxaliplatin (FOLFOX-type): n = 40, deaths = 4, p = 0.502, HR 1.97 (95% CI 0.26-14.89)
- READ, no oxaliplatin or irinotecan record: n = 19, deaths = 1, p = 0.480, HR not estimable (0 deaths in one group)
- READ, with irinotecan, no oxaliplatin (FOLFIRI-type): n = 0
