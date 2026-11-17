**Figure 5. ANCHOR Stage 0: pre-registered, outcome-free mechanism hypotheses.**
(a) Each row is a drug; each column a Reactome pathway that was frozen for at least one drug. Coloured cells are
frozen hypotheses: red = a HIGHER pathway score is expected to go with RESISTANCE, blue = with SENSITIVITY; the
number is how many of the three independent label-free sources supported it (1: expression association across
GDSC1/GDSC2, CTRPv2 and PRISM after removing cancer type and general drug sensitivity; 2: CRISPR co-dependency in
Project Score; 3: STRING network proximity to the drug's targets). Pale cells show the same pathway in another drug
where it was not frozen. Every frozen hypothesis also passed a specificity rule against all drugs and against drugs
of the same mechanism class. No patient data were used. (b) How many pathways passed each source and each step.
Numbers: `tables/frozen_hypotheses.tsv`, `tables/B_stage0_summary.tsv`, `tables/B_candidates_two_or_more_sources.tsv`.
