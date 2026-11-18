# Fig. 1a (redrawn): the whole method in one picture

Read from top to bottom. Blue = lab models. Orange = the protein network step. Green = the machine-learning step.
Pink = the patient test. Grey = extra checks.
1. Lab models (organoids, PDX or cell lines) give two things: gene activity (RNA) and how well the drug kills them.
2. Each gene set from Reactome (a "pathway") gets one activity score per sample (ssGSEA NES).
3. In the STRING protein network we measure how close every pathway is to the drug's target proteins.
   Pathways much closer than random (z <= -1.2816) are kept.
4. Three simple models learn drug response from the kept pathways. The pathway with the largest Ridge weight is the
   candidate biomarker.
5. TCGA patients who received the drug are scored with that pathway, split at the median, and their survival compared.
6. Checks: the same split in untreated patients, random pathway sets, drug-resistant cell lines and known mutations.
This is a drawing; it contains no data.
