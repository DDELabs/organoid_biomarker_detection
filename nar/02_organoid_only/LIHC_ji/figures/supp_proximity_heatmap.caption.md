# Proximity heatmap: Liver cancer organoids (Ji 2023)

Rows: Reactome pathways that pass the paper's proximity cut-off (z <= -1.2816) for at least one of the drugs
tested in this set (22 pathways; the 22 closest are shown). Columns: drugs. Colour: network
proximity z-score between the drug's protein targets and the pathway's genes in the STRING network (authors'
published table; purple/negative = closer than random). A dot marks a proximal pathway. Oxaliplatin has no protein
target, so cisplatin (same platinum class) is used. Proximal pathways per drug: sorafenib 22.
Full table: `tables/proximity_all_drugs.tsv`.
