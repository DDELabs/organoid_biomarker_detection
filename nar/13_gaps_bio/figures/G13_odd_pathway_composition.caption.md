# What is inside the "HIV" pathways among the frozen hypotheses?

Two of the 35 frozen cell-line hypotheses carry virus names (Reactome "APOBEC3G-mediated resistance to HIV-1
infection" and "HIV infection"). Reactome describes how the virus uses human proteins, so these sets are made of human
genes. Here we show which human genes they contain.

- **a** Share of genes from common protein families (by gene-name pattern).
- **b** The five canonical pathways (KEGG, Reactome, BioCarta, PID; MSigDB v6.1) that overlap each set most.
- **c** For all 35 frozen pathways: how the pathway score tracks tumour purity (ABSOLUTE) and an adaptive-immune pathway
  score, as the median correlation within each TCGA cancer type. Large outlined dots = the two "HIV" sets.

- REACTOME_APOBEC3G_MEDIATED_RESISTANCE_TO_HIV1_INFECTION: 12 genes; 
- REACTOME_HIV_INFECTION: 207 genes; 
- REACTOME_APOBEC3G_MEDIATED_RESISTANCE_TO_HIV1_INFECTION vs purity: median rho -0.23 (10th-90th percentile -0.41 to -0.01)
- REACTOME_APOBEC3G_MEDIATED_RESISTANCE_TO_HIV1_INFECTION vs adaptive immune score: median rho 0.48 (10th-90th percentile 0.27 to 0.67)
- REACTOME_HIV_INFECTION vs purity: median rho -0.06 (10th-90th percentile -0.32 to 0.25)
- REACTOME_HIV_INFECTION vs adaptive immune score: median rho 0.48 (10th-90th percentile 0.22 to 0.63)
