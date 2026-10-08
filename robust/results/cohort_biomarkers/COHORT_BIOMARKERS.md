# Pathway biomarkers in liver, kidney, glioma and melanoma (public treated cohorts)

Response cohorts: logistic score-test z (+ = higher in responders); `*` q < 0.1; stab = 200-bootstrap top-10 frequency. Predictive tests: pathway x treatment interaction z in Cox models (+ = treatment benefit increases with the pathway). Exploratory.

## liver (HCC) - Tace (GSE104580; n = 147, responders = 81)

| # | Response pathway | z | stab | Resistance pathway | z | stab |
|---|---|---|---|---|---|---|
| 1 | xenobiotics* | +5.6 | 0.57 | glucose transport* | -5.7 | 0.54 |
| 2 | adenylate cyclase activating pathway* | +5.6 | 0.49 | regulation of hypoxia inducible factor hif by oxygen* | -5.6 | 0.46 |
| 3 | adenylate cyclase inhibitory pathway* | +5.6 | 0.48 | shc1 events in erbb4 signaling* | -5.5 | 0.43 |
| 4 | ethanol oxidation* | +5.3 | 0.36 | oxygen dependent proline hydroxylation of hypoxia induc* | -5.4 | 0.36 |
| 5 | phase1 functionalization of compounds* | +5.3 | 0.34 | egfr downregulation* | -5.3 | 0.39 |

## kidney (ccRCC) - Sunitinib (EMTAB3267; n = 43, responders = 19)

| # | Response pathway | z | stab | Resistance pathway | z | stab |
|---|---|---|---|---|---|---|
| 1 | regulation of pyruvate dehydrogenase pdh complex | +2.6 | 0.35 | vitamin b5 pantothenate metabolism | -2.6 | 0.28 |
| 2 | norepinephrine neurotransmitter release cycle | +2.5 | 0.19 | acyl chain remodelling of pi | -2.5 | 0.23 |
| 3 | rap1 signalling | +2.5 | 0.18 | microrna mirna biogenesis | -2.5 | 0.28 |
| 4 | recycling of bile acids and salts | +2.4 | 0.14 | elevation of cytosolic ca2 levels | -2.5 | 0.20 |
| 5 | nuclear receptor transcription pathway | +2.4 | 0.12 | kinesins | -2.5 | 0.14 |

## kidney (ccRCC) - Nivolumab (BRAUN2020_CHECKMATE; n = 172, responders = 39)

| # | Response pathway | z | stab | Resistance pathway | z | stab |
|---|---|---|---|---|---|---|
| 1 | platelet adhesion to exposed collagen | +2.8 | 0.27 | regulatory rna pathways | -2.8 | 0.34 |
| 2 | the nlrp3 inflammasome | +2.7 | 0.28 | transport of vitamins nucleosides and related molecules | -2.8 | 0.32 |
| 3 | signaling by gpcr | +2.6 | 0.24 | microrna mirna biogenesis | -2.7 | 0.29 |
| 4 | gpcr ligand binding | +2.4 | 0.16 | peroxisomal lipid metabolism | -2.6 | 0.17 |
| 5 | gpcr downstream signaling | +2.4 | 0.12 | rora activates circadian expression | -2.5 | 0.28 |

## kidney (ccRCC) - Everolimus (BRAUN2020_CHECKMATE; n = 109, responders = 5)

| # | Response pathway | z | stab | Resistance pathway | z | stab |
|---|---|---|---|---|---|---|
| 1 | cgmp effects | +2.9 | 0.63 | highly calcium permeable postsynaptic nicotinic acetylc | -1.9 | 0.07 |
| 2 | rap1 signalling | +2.9 | 0.57 | androgen biosynthesis | -1.7 | 0.04 |
| 3 | nitric oxide stimulates guanylate cyclase | +2.8 | 0.62 | acetylcholine binding and downstream events | -1.6 | 0.01 |
| 4 | platelet homeostasis | +2.5 | 0.36 | recruitment of mitotic centrosome proteins and complexe | -1.5 | 0.09 |
| 5 | regulation of kit signaling | +2.3 | 0.33 | downregulation of erbb2 erbb3 signaling | -1.5 | 0.02 |

## melanoma - Nivolumab (GSE91061; n = 49, responders = 10)

| # | Response pathway | z | stab | Resistance pathway | z | stab |
|---|---|---|---|---|---|---|
| 1 | highly calcium permeable postsynaptic nicotinic acetylc | +2.4 | 0.33 | retrograde neurotrophin signalling | -2.6 | 0.28 |
| 2 | irak2 mediated activation of tak1 complex upon tlr7 8 o | +2.3 | 0.12 | post translational modification synthesis of gpi anchor | -2.4 | 0.14 |
| 3 | purine metabolism | +2.2 | 0.17 | metabolism of proteins | -2.4 | 0.14 |
| 4 | p2y receptors | +2.2 | 0.12 | synthesis of glycosylphosphatidylinositol gpi | -2.3 | 0.09 |
| 5 | rap1 signalling | +2.2 | 0.15 | hormone ligand binding receptors | -2.1 | 0.15 |

## melanoma - Pembrolizumab (GSE78220; n = 28, responders = 15)

| # | Response pathway | z | stab | Resistance pathway | z | stab |
|---|---|---|---|---|---|---|
| 1 | phospholipid metabolism | +2.6 | 0.15 | yap1 and wwtr1 taz stimulated gene expression | -2.9 | 0.38 |
| 2 | glycogen breakdown glycogenolysis | +2.6 | 0.26 | hs gag degradation | -2.9 | 0.37 |
| 3 | endosomal sorting complex required for transport escrt | +2.5 | 0.24 | g alpha s signalling events | -2.8 | 0.31 |
| 4 | synthesis of pa | +2.5 | 0.20 | hs gag biosynthesis | -2.7 | 0.20 |
| 5 | glycerophospholipid biosynthesis | +2.4 | 0.12 | class b 2 secretin family receptors | -2.6 | 0.17 |

## Predictive (treatment x pathway interaction) tests

**kidney CM-025 nivolumab vs everolimus (PFS)** (n = 250; min q = 0.70)

| Pathway | interaction z | q |
|---|---|---|
| nrif signals cell death from the nucleus | -2.85 | 0.70 |
| gaba synthesis release reuptake and degradation | -2.83 | 0.70 |
| signaling by tgf beta receptor complex | +2.66 | 0.70 |
| regulatory rna pathways | -2.42 | 0.70 |
| signaling by notch3 | -2.41 | 0.70 |
| signaling by notch2 | -2.38 | 0.70 |
| prostanoid ligand receptors | +2.36 | 0.70 |
| glucagon type ligand receptors | +2.35 | 0.70 |
| nucleotide excision repair | -2.33 | 0.70 |
| neurotransmitter release cycle | -2.33 | 0.70 |

**CGGA all glioma TMZ vs no TMZ (OS)** (n = 636; min q = 0.99)

| Pathway | interaction z | q |
|---|---|---|
| sulfur amino acid metabolism | +2.78 | 0.99 |
| regulation of insulin secretion by acetylcholine | -2.56 | 0.99 |
| tandem pore domain potassium channels | -2.21 | 0.99 |
| unblocking of nmda receptor glutamate binding and activation | -2.17 | 0.99 |
| platelet homeostasis | -2.12 | 0.99 |
| activation of nmda receptor upon glutamate binding and posts | -2.04 | 0.99 |
| creb phosphorylation through the activation of ras | -2.03 | 0.99 |
| g0 and early g1 | +1.99 | 0.99 |
| trafficking of glur2 containing ampa receptors | -1.99 | 0.99 |
| gaba b receptor activation | -1.95 | 0.99 |

**CGGA GBM (WHO IV) TMZ vs no TMZ (OS)** (n = 226; min q = 0.97)

| Pathway | interaction z | q |
|---|---|---|
| metabolism of vitamins and cofactors | +2.48 | 0.97 |
| the activation of arylsulfatases | +2.25 | 0.97 |
| cs ds degradation | +2.11 | 0.97 |
| metabolism of lipids and lipoproteins | +2.09 | 0.97 |
| sulfur amino acid metabolism | +2.01 | 0.97 |
| elevation of cytosolic ca2 levels | +1.95 | 0.97 |
| heparan sulfate heparin hs gag metabolism | +1.90 | 0.97 |
| fatty acyl coa biosynthesis | +1.88 | 0.97 |
| alpha linolenic acid ala metabolism | +1.86 | 0.97 |
| regulation of beta cell development | -1.85 | 0.97 |

**GSE7696 GBM RT+TMZ vs RT (OS)** (n = 70; min q = 0.42)

| Pathway | interaction z | q |
|---|---|---|
| gpcr ligand binding | +3.27 | 0.42 |
| class a1 rhodopsin like receptors | +3.11 | 0.42 |
| signaling by gpcr | +3.04 | 0.42 |
| synthesis of glycosylphosphatidylinositol gpi | -3.00 | 0.42 |
| amino acid transport across the plasma membrane | +2.86 | 0.42 |
| gpcr downstream signaling | +2.82 | 0.42 |
| amine ligand binding receptors | +2.77 | 0.42 |
| peptide ligand binding receptors | +2.74 | 0.42 |
| serotonin receptors | +2.64 | 0.42 |
| transport of inorganic cations anions and amino acids oligop | +2.63 | 0.42 |

**GBM meta (CGGA WHO IV + GSE7696), Stouffer** (n = 296; min q = 0.67)

| Pathway | interaction z | q |
|---|---|---|
| defensins | +2.82 | 0.67 |
| amino acid transport across the plasma membrane | +2.82 | 0.67 |
| gpcr ligand binding | +2.58 | 0.67 |
| transport of inorganic cations anions and amino acids oligop | +2.51 | 0.67 |
| il1 signaling | +2.46 | 0.67 |
| signaling by gpcr | +2.39 | 0.67 |
| lysosome vesicle biogenesis | +2.37 | 0.67 |
| class a1 rhodopsin like receptors | +2.34 | 0.67 |
| serotonin receptors | +2.29 | 0.67 |
| elevation of cytosolic ca2 levels | +2.27 | 0.67 |

