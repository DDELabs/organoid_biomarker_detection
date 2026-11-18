**Does a pathway score tell us who benefits from chemotherapy after lung surgery?**
Each row is one lung cohort that was kept out of training (the model was learned on the other cohorts).
**a**, the "interaction" hazard ratio (HR, how much the chemo effect changes per 1 standard deviation of the score; the dashed line at 1 means no change)
with 95% confidence interval. Blue: the score uses all 673 Reactome pathways. Orange: the score uses only the 47 pathways close to cisplatin's targets in the STRING network (the paper's shortlist).
The bottom row pools the tested cohorts. GSE14814 is the randomised trial, so it is the cleanest test.
**b**, the plain chemo-vs-no-chemo hazard ratio inside each half of the patients, split at the median predicted benefit (all-pathway score). If the score works, green should sit left of red.
n = patients, events = deaths in the held-out cohort.
