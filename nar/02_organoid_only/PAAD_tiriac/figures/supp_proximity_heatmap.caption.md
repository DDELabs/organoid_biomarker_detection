# Proximity heatmap: Pancreatic organoids (Tiriac 2018)

Rows: Reactome pathways that pass the paper's proximity cut-off (z <= -1.2816) for at least one of the drugs
tested in this set (150 pathways; the 60 closest are shown). Columns: drugs. Colour: network
proximity z-score between the drug's protein targets and the pathway's genes in the STRING network (authors'
published table; purple/negative = closer than random). A dot marks a proximal pathway. Oxaliplatin has no protein
target, so cisplatin (same platinum class) is used. Proximal pathways per drug: fluorouracil 38, gemcitabine 45, oxaliplatin 47, irinotecan 78.
Full table: `tables/proximity_all_drugs.tsv`.
