# Top-10 pathway biomarkers per cancer x drug - relaxed criteria

Pairs not covered by the standard analysis (TOP10_BY_CANCER.md). Tier A = RECIST (n >= 10, >= 3 per class); B = TCGA first-course treatment outcome among first-line recipients; C = progression-free interval among treated patients (Cox; positive z = longer PFI). For tier C, `pred` = treated minus untreated effect z (crude predictive check; untreated = no recorded systemic therapy). `*` = local q < 0.1. Combined z uses the ATLAS pan-cancer prior when the drug is in ATLAS. Exploratory.

## BLCA - bladder

**Carboplatin** - tier A (RECIST response), n = 18, responders = 8

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of rna | +1.8 | +3.5 | 0.00 | prostanoid ligand receptors | -1.9 | -3.6 | 0.06 |
| 2 | creb phosphorylation through the activation of camkii | +1.0 | +3.2 | 0.00 | nef mediated downregulation of mhc class i complex cell | -2.6 | -3.6 | 0.18 |
| 3 | signaling by fgfr3 mutants | +1.6 | +2.9 | 0.03 | effects of pip2 hydrolysis | -2.1 | -3.3 | 0.14 |
| 4 | regulation of beta cell development | +2.1 | +2.8 | 0.05 | rip mediated nfkb activation via dai | -1.6 | -3.2 | 0.00 |
| 5 | g0 and early g1 | +2.2 | +2.8 | 0.02 | semaphorin interactions | -2.2 | -3.2 | 0.03 |

**Doxorubicin** - tier A (RECIST response), n = 14, responders = 9

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | synthesis of pips at the early endosome membrane | +2.0 | +3.2 | 0.06 | phospholipase c mediated cascade | -1.7 | -3.0 | 0.00 |
| 2 | regulation of the fanconi anemia pathway | +0.9 | +3.0 | 0.00 | negative regulation of fgfr signaling | -2.5 | -2.7 | 0.20 |
| 3 | synthesis of pips at the late endosome membrane | +2.2 | +2.8 | 0.13 | amino acid synthesis and interconversion transamination | -1.4 | -2.6 | 0.04 |
| 4 | synthesis of pips at the golgi membrane | +1.9 | +2.8 | 0.19 | acyl chain remodelling of pc | -0.7 | -2.5 | 0.00 |
| 5 | ctla4 inhibitory signaling | +1.1 | +2.6 | 0.00 | fgfr ligand binding and activation | -2.4 | -2.5 | 0.16 |

**Methotrexate** - tier A (RECIST response), n = 10, responders = 5

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | activation of chaperones by atf6 alpha | +2.6 | +2.6 | 0.46 | amine derived hormones | -2.4 | -2.4 | 0.26 |
| 2 | synthesis of very long chain fatty acyl coas | +2.5 | +2.5 | 0.47 | pi 3k cascade | -2.2 | -2.2 | 0.11 |
| 3 | synthesis of pips at the late endosome membrane | +2.4 | +2.4 | 0.24 | pi3k cascade | -2.1 | -2.1 | 0.06 |
| 4 | fatty acyl coa biosynthesis | +2.4 | +2.4 | 0.27 | signaling by insulin receptor | -2.0 | -2.0 | 0.05 |
| 5 | ikk complex recruitment mediated by rip1 | +2.4 | +2.4 | 0.28 | insulin receptor signalling cascade | -2.0 | -2.0 | 0.01 |

**Vinblastine** - tier A (RECIST response), n = 10, responders = 5

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | activation of chaperones by atf6 alpha | +2.6 | +2.6 | 0.47 | amine derived hormones | -2.4 | -2.4 | 0.24 |
| 2 | synthesis of very long chain fatty acyl coas | +2.5 | +2.5 | 0.40 | pi 3k cascade | -2.2 | -2.2 | 0.12 |
| 3 | synthesis of pips at the late endosome membrane | +2.4 | +2.4 | 0.22 | pi3k cascade | -2.1 | -2.1 | 0.07 |
| 4 | fatty acyl coa biosynthesis | +2.4 | +2.4 | 0.31 | signaling by insulin receptor | -2.0 | -2.0 | 0.09 |
| 5 | ikk complex recruitment mediated by rip1 | +2.4 | +2.4 | 0.35 | insulin receptor signalling cascade | -2.0 | -2.0 | 0.01 |

## BRCA - breast

**Anastrozole** - tier A (RECIST response), n = 22, responders = 17

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | creation of c4 and c2 activators | +2.8 | +2.8 | 0.28 | opsins | -2.7 | -2.7 | 0.22 |
| 2 | grb2 sos provides linkage to mapk signaling for intergr | +2.8 | +2.8 | 0.35 | akt phosphorylates targets in the cytosol | -2.7 | -2.7 | 0.31 |
| 3 | innate immune system | +2.7 | +2.7 | 0.27 | digestion of dietary carbohydrate | -2.6 | -2.6 | 0.21 |
| 4 | il 7 signaling | +2.7 | +2.7 | 0.13 | regulation of gene expression in beta cells | -2.5 | -2.5 | 0.21 |
| 5 | p130cas linkage to mapk signaling for integrins | +2.7 | +2.7 | 0.21 | metabolism of proteins | -2.5 | -2.5 | 0.12 |

**Exemestane** - tier C (PFI among treated (Cox)), n = 48, events = 10

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | n glycan trimming in the er and calnexin calreticulin c | +3.5 | +3.5 | 0.70 | +4.1 | cgmp effects | -2.4 | -2.4 | 0.36 | -2.0 |
| 2 | calnexin calreticulin cycle | +3.4 | +3.4 | 0.65 | +3.9 | amine ligand binding receptors | -2.0 | -2.0 | 0.15 | -1.7 |
| 3 | sphingolipid metabolism | +2.6 | +2.6 | 0.29 | +2.7 | rna pol iii transcription | -2.0 | -2.0 | 0.06 | -1.7 |
| 4 | erks are inactivated | +2.5 | +2.5 | 0.18 | +2.6 | recycling of bile acids and salts | -1.9 | -1.9 | 0.14 | -1.9 |
| 5 | activation of chaperone genes by xbp1s | +2.4 | +2.4 | 0.19 | +2.7 | nitric oxide stimulates guanylate cyclase | -1.9 | -1.9 | 0.22 | -1.4 |

**Fluorouracil** - tier C (PFI among treated (Cox)), n = 101, events = 10

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | mrna decay by 3 to 5 exoribonuclease | +2.0 | +3.5 | 0.00 | +1.6 | chondroitin sulfate dermatan sulfate metabolism | -2.2 | -3.7 | 0.08 | -1.4 |
| 2 | metabolism of rna | +1.5 | +3.2 | 0.00 | +1.4 | gap junction degradation | -1.8 | -3.6 | 0.03 | -1.7 |
| 3 | alpha linolenic acid ala metabolism | +1.5 | +3.0 | 0.03 | +2.2 | other semaphorin interactions | -2.0 | -3.4 | 0.07 | -1.1 |
| 4 | destabilization of mrna by tristetraprolin ttp | +2.5 | +3.0 | 0.12 | +2.0 | signaling by notch | -3.3 | -3.3 | 0.53 | -2.9 |
| 5 | regulation of mrna stability by proteins that bind au r | +1.7 | +2.8 | 0.06 | +1.6 | cgmp effects | -2.6 | -3.1 | 0.21 | -2.2 |

**Letrozole** - tier C (PFI among treated (Cox)), n = 74, events = 12

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | nod1 2 signaling pathway | +3.1 | +3.1 | 0.37 | +2.2 | inhibition of replication initiation of damaged dna by  | -2.8 | -2.8 | 0.30 | -2.4 |
| 2 | nucleotide binding domain leucine rich repeat containin | +2.9 | +2.9 | 0.23 | +1.7 | enos activation and regulation | -2.7 | -2.7 | 0.23 | -2.4 |
| 3 | the nlrp3 inflammasome | +2.7 | +2.7 | 0.15 | +1.7 | pkb mediated events | -2.4 | -2.4 | 0.20 | -2.4 |
| 4 | il 7 signaling | +2.7 | +2.7 | 0.17 | +1.4 | akt phosphorylates targets in the cytosol | -2.4 | -2.4 | 0.09 | -2.9 |
| 5 | amino acid transport across the plasma membrane | +2.6 | +2.6 | 0.23 | +3.1 | olfactory signaling pathway | -2.3 | -2.3 | 0.17 | -1.9 |

## CESC - cervical

**Carboplatin** - tier A (RECIST response), n = 11, responders = 6

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | base free sugar phosphate removal via the single nucleo | +1.7 | +2.8 | 0.05 | prostanoid ligand receptors | -1.2 | -3.2 | 0.00 |
| 2 | transport of glucose and other sugars bile salts and or | +1.2 | +2.6 | 0.00 | phospholipase c mediated cascade | -1.2 | -3.0 | 0.00 |
| 3 | n glycan antennae elongation | +2.0 | +2.5 | 0.08 | dag and ip3 signaling | -2.1 | -2.9 | 0.17 |
| 4 | glucuronidation | +1.8 | +2.5 | 0.01 | acyl chain remodelling of pc | -1.3 | -2.9 | 0.01 |
| 5 | transport of organic anions | +0.9 | +2.3 | 0.02 | downregulation of smad2 3 smad4 transcriptional activit | -1.2 | -2.8 | 0.00 |

**Paclitaxel** - tier A (RECIST response), n = 12, responders = 7

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | phospholipid metabolism | +2.2 | +3.5 | 0.14 | prostanoid ligand receptors | -1.8 | -3.8 | 0.01 |
| 2 | sphingolipid metabolism | +2.2 | +2.7 | 0.20 | signaling by bmp | -2.5 | -3.6 | 0.26 |
| 3 | glucuronidation | +1.9 | +2.5 | 0.05 | phospholipase c mediated cascade | -1.6 | -3.3 | 0.00 |
| 4 | n glycan antennae elongation | +1.4 | +2.4 | 0.00 | integration of provirus | -2.4 | -3.3 | 0.32 |
| 5 | synthesis of bile acids and bile salts via 24 hydroxych | +0.7 | +2.4 | 0.00 | regulation of water balance by renal aquaporins | -0.6 | -3.1 | 0.00 |

## COAD - colon

**Capecitabine** - tier A (RECIST response), n = 19, responders = 12

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | alpha linolenic acid ala metabolism | +1.5 | +3.1 | 0.00 | pol switching | -2.2 | -3.6 | 0.08 |
| 2 | adherens junctions interactions | +1.9 | +2.9 | 0.05 | telomere maintenance | -2.3 | -3.5 | 0.08 |
| 3 | ctla4 inhibitory signaling | +1.3 | +2.9 | 0.01 | microrna mirna biogenesis | -2.4 | -3.3 | 0.11 |
| 4 | signaling by fgfr1 mutants | +1.2 | +2.9 | 0.04 | rna pol iii transcription termination | -2.0 | -3.2 | 0.11 |
| 5 | map kinase activation in tlr cascade | +1.7 | +2.8 | 0.05 | insulin receptor recycling | -3.0 | -3.1 | 0.52 |

**Irinotecan** - tier A (RECIST response), n = 19, responders = 6

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | glucose metabolism | +1.7 | +3.3 | 0.04 | mtorc1 mediated signalling | -3.0 | -3.6 | 0.65 |
| 2 | extrinsic pathway for apoptosis | +2.6 | +2.9 | 0.48 | metabolism of proteins | -1.0 | -3.0 | 0.02 |
| 3 | gluconeogenesis | +1.8 | +2.8 | 0.09 | integration of provirus | -1.7 | -2.8 | 0.02 |
| 4 | mrna decay by 3 to 5 exoribonuclease | +1.1 | +2.7 | 0.02 | regulation of water balance by renal aquaporins | -0.6 | -2.7 | 0.00 |
| 5 | inflammasomes | +1.6 | +2.7 | 0.03 | gap junction degradation | -0.8 | -2.7 | 0.02 |

## DLBC - DLBC

**Cyclophosphamide** - tier B (first-course outcome), n = 39, responders = 32

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | glucose metabolism | +1.6 | +3.1 | 0.00 | rip mediated nfkb activation via dai | -2.9 | -3.6 | 0.40 |
| 2 | resolution of ap sites via the multiple nucleotide patc | +2.4 | +3.1 | 0.02 | p38mapk events | -3.2 | -3.3 | 0.47 |
| 3 | gluconeogenesis | +2.0 | +3.0 | 0.03 | rig i mda5 mediated induction of ifn alpha beta pathway | -2.2 | -3.1 | 0.03 |
| 4 | base free sugar phosphate removal via the single nucleo | +1.9 | +2.9 | 0.01 | insulin synthesis and processing | -2.6 | -3.0 | 0.11 |
| 5 | conversion from apc c cdc20 to apc c cdh1 in late anaph | +1.5 | +2.5 | 0.00 | traf6 mediated nfkb activation | -3.4 | -3.0 | 0.64 |

**Doxorubicin** - tier B (first-course outcome), n = 38, responders = 31

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | base free sugar phosphate removal via the single nucleo | +1.8 | +3.2 | 0.00 | rip mediated nfkb activation via dai | -2.8 | -3.7 | 0.42 |
| 2 | gluconeogenesis | +2.2 | +3.2 | 0.02 | rig i mda5 mediated induction of ifn alpha beta pathway | -2.0 | -3.4 | 0.04 |
| 3 | resolution of ap sites via the multiple nucleotide patc | +2.3 | +3.1 | 0.01 | p38mapk events | -3.2 | -3.4 | 0.45 |
| 4 | signaling by fgfr3 mutants | +1.9 | +2.9 | 0.00 | tak1 activates nfkb by phosphorylation and activation o | -2.3 | -2.9 | 0.12 |
| 5 | glucose metabolism | +1.8 | +2.7 | 0.02 | traf6 mediated nfkb activation | -3.3 | -2.8 | 0.60 |

**Rituximab** - tier B (first-course outcome), n = 30, responders = 25

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | ligand gated ion channel transport | +2.6 | +2.6 | 0.22 | peroxisomal lipid metabolism* | -3.7 | -3.7 | 0.72 |
| 2 | g2 m dna damage checkpoint | +2.4 | +2.4 | 0.04 | rip mediated nfkb activation via dai* | -3.6 | -3.6 | 0.74 |
| 3 | class c 3 metabotropic glutamate pheromone receptors | +2.3 | +2.3 | 0.06 | traf6 mediated nfkb activation* | -3.6 | -3.6 | 0.72 |
| 4 | repair synthesis for gap filling by dna pol in tc ner | +2.3 | +2.3 | 0.02 | zinc transporters* | -3.5 | -3.5 | 0.58 |
| 5 | base excision repair | +2.3 | +2.3 | 0.00 | metal ion slc transporters* | -3.5 | -3.5 | 0.60 |

**Vincristine** - tier B (first-course outcome), n = 39, responders = 32

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | repair synthesis for gap filling by dna pol in tc ner | +2.8 | +2.8 | 0.25 | traf6 mediated nfkb activation | -3.4 | -3.4 | 0.62 |
| 2 | activation of atr in response to replication stress | +2.7 | +2.7 | 0.32 | p38mapk events | -3.3 | -3.3 | 0.57 |
| 3 | g2 m checkpoints | +2.7 | +2.7 | 0.26 | termination of o glycan biosynthesis | -3.2 | -3.2 | 0.45 |
| 4 | lagging strand synthesis | +2.5 | +2.5 | 0.02 | synthesis of pe | -2.9 | -2.9 | 0.37 |
| 5 | pol switching | +2.5 | +2.5 | 0.02 | rip mediated nfkb activation via dai | -2.9 | -2.9 | 0.35 |

## ESCA - oesophageal

**Cisplatin** - tier A (RECIST response), n = 17, responders = 14

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | regulation of the fanconi anemia pathway | +2.8 | +3.7 | 0.38 | rip mediated nfkb activation via dai | -2.5 | -3.9 | 0.08 |
| 2 | metabolism of rna | +1.7 | +3.5 | 0.00 | synthesis of very long chain fatty acyl coas | -3.0 | -3.2 | 0.65 |
| 3 | cell cycle | +1.3 | +3.2 | 0.00 | oxygen dependent proline hydroxylation of hypoxia induc | -2.3 | -3.2 | 0.10 |
| 4 | regulation of beta cell development | +2.2 | +3.0 | 0.03 | rig i mda5 mediated induction of ifn alpha beta pathway | -1.4 | -3.1 | 0.00 |
| 5 | resolution of ap sites via the multiple nucleotide patc | +2.1 | +2.8 | 0.00 | transferrin endocytosis and recycling | -3.0 | -2.9 | 0.64 |

**Fluorouracil** - tier A (RECIST response), n = 11, responders = 8

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | jnk c jun kinases phosphorylation and activation mediat | +2.0 | +3.0 | 0.03 | metabolism of proteins | -1.6 | -3.4 | 0.00 |
| 2 | inhibition of insulin secretion by adrenaline noradrena | +2.2 | +2.9 | 0.12 | ptm gamma carboxylation hypusine formation and arylsulf | -2.9 | -3.0 | 0.53 |
| 3 | signaling by fgfr1 mutants | +1.1 | +2.8 | 0.00 | metabolism of vitamins and cofactors | -2.5 | -2.9 | 0.26 |
| 4 | downstream signal transduction | +1.2 | +2.7 | 0.02 | reversible hydration of carbon dioxide | -1.6 | -2.9 | 0.00 |
| 5 | regulation of the fanconi anemia pathway | +1.4 | +2.7 | 0.04 | gamma carboxylation transport and amino terminal cleava | -2.2 | -2.7 | 0.05 |

## GBM - glioblastoma

**Bevacizumab** - tier C (PFI among treated (Cox)), n = 31, events = 31

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of rna | +1.9 | +3.3 | 0.00 | +0.6 | rip mediated nfkb activation via dai | -2.4 | -3.3 | 0.03 | -0.3 |
| 2 | conversion from apc c cdc20 to apc c cdh1 in late anaph | +2.3 | +3.1 | 0.10 | +1.4 | p75ntr recruits signalling complexes | -2.6 | -3.3 | 0.20 | -1.7 |
| 3 | mrna decay by 3 to 5 exoribonuclease | +1.5 | +3.0 | 0.00 | +0.2 | semaphorin interactions | -1.8 | -3.2 | 0.03 | -1.1 |
| 4 | base free sugar phosphate removal via the single nucleo | +1.7 | +2.9 | 0.00 | +1.9 | p75ntr signals via nfkb | -2.9 | -3.1 | 0.29 | -1.7 |
| 5 | inhibition of voltage gated ca2 channels via gbeta gamm | +1.5 | +2.8 | 0.01 | +2.1 | o linked glycosylation of mucins | -2.7 | -3.1 | 0.12 | -1.5 |

**Temozolomide** - tier C (PFI among treated (Cox)), n = 100, events = 79

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | mrna decay by 3 to 5 exoribonuclease | +1.3 | +3.3 | 0.01 | -0.5 | acyl chain remodelling of pc | -2.2 | -3.6 | 0.07 | -0.3 |
| 2 | downstream signal transduction | +1.8 | +3.1 | 0.03 | +1.5 | post translational protein modification | -2.7 | -3.3 | 0.08 | +0.2 |
| 3 | signaling by fgfr1 mutants | +1.4 | +2.9 | 0.00 | -0.9 | metabolism of proteins | -1.6 | -3.3 | 0.00 | -0.2 |
| 4 | notch1 intracellular domain regulates transcription | +2.0 | +2.6 | 0.06 | -1.3 | ptm gamma carboxylation hypusine formation and arylsulf | -2.9 | -3.2 | 0.24 | -0.1 |
| 5 | metabolism of rna | +0.8 | +2.5 | 0.00 | -0.8 | acyl chain remodelling of pg | -1.8 | -2.8 | 0.00 | -1.3 |

## HNSC - head & neck

**Cetuximab** - tier A (RECIST response), n = 20, responders = 11

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | g0 and early g1 | +2.7 | +2.7 | 0.39 | dscam interactions | -2.9 | -2.9 | 0.48 |
| 2 | antiviral mechanism by ifn stimulated genes | +2.6 | +2.6 | 0.30 | abc family proteins mediated transport | -2.2 | -2.2 | 0.11 |
| 3 | pecam1 interactions | +2.5 | +2.5 | 0.32 | olfactory signaling pathway | -2.2 | -2.2 | 0.14 |
| 4 | regulation of the fanconi anemia pathway | +2.4 | +2.4 | 0.13 | membrane binding and targetting of gag proteins | -2.1 | -2.1 | 0.05 |
| 5 | host interactions of hiv factors | +2.3 | +2.3 | 0.28 | abca transporters in lipid homeostasis | -1.9 | -1.9 | 0.03 |

## KIRC - kidney clear cell

**Sunitinib** - tier C (PFI among treated (Cox)), n = 27, events = 24

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | peroxisomal lipid metabolism | +2.9 | +2.9 | 0.43 | +0.9 | fgfr1 ligand binding and activation | -3.3 | -3.3 | 0.53 | -3.3 |
| 2 | synthesis of pc | +2.8 | +2.8 | 0.36 | +2.2 | signaling by activated point mutants of fgfr1 | -2.7 | -2.7 | 0.31 | -2.6 |
| 3 | synthesis of bile acids and bile salts via 7alpha hydro | +2.6 | +2.6 | 0.22 | +1.2 | regulatory rna pathways | -2.4 | -2.4 | 0.23 | -0.7 |
| 4 | synthesis of bile acids and bile salts | +2.5 | +2.5 | 0.17 | +0.6 | microrna mirna biogenesis | -2.2 | -2.2 | 0.14 | -0.0 |
| 5 | mitochondrial fatty acid beta oxidation | +2.5 | +2.5 | 0.11 | +0.1 | activated point mutants of fgfr2 | -2.2 | -2.2 | 0.20 | -2.2 |

## LGG - lower-grade glioma

**Bevacizumab** - tier C (PFI among treated (Cox)), n = 46, events = 44

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | jnk c jun kinases phosphorylation and activation mediat | +2.9 | +3.7 | 0.16 | +0.4 | branched chain amino acid catabolism | -3.4 | -3.7 | 0.38 | -3.1 |
| 2 | platelet calcium homeostasis | +2.3 | +3.6 | 0.04 | +1.6 | mhc class ii antigen presentation | -1.6 | -3.1 | 0.01 | +0.5 |
| 3 | downstream signal transduction | +2.1 | +3.6 | 0.00 | +2.1 | tight junction interactions | -2.6 | -3.1 | 0.13 | -0.9 |
| 4 | n glycan antennae elongation | +2.8 | +3.3 | 0.16 | +1.2 | glutathione conjugation | -2.7 | -3.1 | 0.08 | -2.0 |
| 5 | developmental biology | +2.3 | +3.2 | 0.04 | +2.6 | abca transporters in lipid homeostasis | -2.6 | -3.0 | 0.09 | -1.0 |

**Irinotecan** - tier C (PFI among treated (Cox)), n = 21, events = 21

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | map kinase activation in tlr cascade | +1.9 | +3.1 | 0.07 | +3.1 | acyl chain remodelling of pc | -1.0 | -2.7 | 0.00 | +0.3 |
| 2 | creb phosphorylation through the activation of camkii | +1.8 | +3.1 | 0.01 | +1.2 | rip mediated nfkb activation via dai | -1.2 | -2.5 | 0.02 | -0.1 |
| 3 | downstream signal transduction | +1.5 | +3.0 | 0.01 | +1.7 | branched chain amino acid catabolism | -1.7 | -2.5 | 0.05 | -2.1 |
| 4 | platelet calcium homeostasis | +1.3 | +2.9 | 0.01 | +1.1 | rig i mda5 mediated induction of ifn alpha beta pathway | -1.4 | -2.4 | 0.02 | -0.0 |
| 5 | cell death signalling via nrage nrif and nade | +3.1 | +2.9 | 0.50 | +3.1 | integration of provirus | -1.0 | -2.3 | 0.11 | -1.1 |

**Lomustine** - tier A (RECIST response), n = 23, responders = 3

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | myogenesis | +3.1 | +3.1 | 0.56 | transport of vitamins nucleosides and related molecules | -3.1 | -3.1 | 0.49 |
| 2 | signalling to erks | +2.9 | +2.9 | 0.35 | downstream signaling events of b cell receptor bcr | -3.0 | -3.0 | 0.53 |
| 3 | signalling to ras | +2.8 | +2.8 | 0.17 | creation of c4 and c2 activators | -2.9 | -2.9 | 0.46 |
| 4 | amino acid synthesis and interconversion transamination | +2.6 | +2.6 | 0.14 | alpha linolenic acid ala metabolism | -2.9 | -2.9 | 0.51 |
| 5 | signalling to p38 via rit and rin | +2.6 | +2.6 | 0.16 | ikk complex recruitment mediated by rip1 | -2.9 | -2.9 | 0.42 |

**Procarbazine** - tier A (RECIST response), n = 12, responders = 3

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | signalling to p38 via rit and rin | +2.6 | +2.6 | 0.46 | downstream signaling events of b cell receptor bcr | -2.5 | -2.5 | 0.31 |
| 2 | myogenesis | +2.4 | +2.4 | 0.32 | ikk complex recruitment mediated by rip1 | -2.5 | -2.5 | 0.33 |
| 3 | recruitment of mitotic centrosome proteins and complexe | +2.2 | +2.2 | 0.10 | transport of vitamins nucleosides and related molecules | -2.5 | -2.5 | 0.39 |
| 4 | hormone ligand binding receptors | +2.2 | +2.2 | 0.21 | synthesis secretion and deacylation of ghrelin | -2.5 | -2.5 | 0.24 |
| 5 | signalling to erks | +2.1 | +2.1 | 0.05 | regulated proteolysis of p75ntr | -2.5 | -2.5 | 0.34 |

**Vincristine** - tier A (RECIST response), n = 10, responders = 3

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | raf map kinase cascade | +2.8 | +2.8 | 0.52 | synthesis secretion and deacylation of ghrelin | -2.8 | -2.8 | 0.64 |
| 2 | myogenesis | +2.6 | +2.6 | 0.25 | signaling by the b cell receptor bcr | -2.8 | -2.8 | 0.48 |
| 3 | signalling to p38 via rit and rin | +2.6 | +2.6 | 0.21 | trafficking and processing of endosomal tlr | -2.8 | -2.8 | 0.38 |
| 4 | arms mediated activation | +2.4 | +2.4 | 0.08 | synthesis secretion and inactivation of gip | -2.7 | -2.7 | 0.32 |
| 5 | shc1 events in egfr signaling | +2.4 | +2.4 | 0.07 | transport of vitamins nucleosides and related molecules | -2.7 | -2.7 | 0.28 |

## LIHC - liver

**Sorafenib** - tier A (RECIST response), n = 17, responders = 3

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | pre notch transcription and translation | +2.9 | +2.9 | 0.66 | vitamin b5 pantothenate metabolism | -2.4 | -2.4 | 0.24 |
| 2 | signaling by robo receptor | +2.9 | +2.9 | 0.56 | steroid hormones | -2.1 | -2.1 | 0.03 |
| 3 | notch hlh transcription pathway | +2.6 | +2.6 | 0.32 | androgen biosynthesis | -2.1 | -2.1 | 0.00 |
| 4 | activation of bh3 only proteins | +2.5 | +2.5 | 0.32 | metabolism of steroid hormones and vitamins a and d | -2.1 | -2.1 | 0.02 |
| 5 | presynaptic nicotinic acetylcholine receptors | +2.4 | +2.4 | 0.22 | metabolism of vitamins and cofactors | -2.0 | -2.0 | 0.14 |

## LUAD - lung adeno

**Docetaxel** - tier A (RECIST response), n = 12, responders = 4

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | trafficking of glur2 containing ampa receptors | +2.6 | +3.1 | 0.47 | membrane trafficking | -1.6 | -3.0 | 0.00 |
| 2 | sema3a pak dependent axon repulsion | +2.2 | +3.1 | 0.14 | amino acid synthesis and interconversion transamination | -0.9 | -2.9 | 0.00 |
| 3 | p2y receptors | +1.6 | +2.9 | 0.05 | circadian clock | -2.1 | -2.7 | 0.16 |
| 4 | metabolism of carbohydrates | +1.8 | +2.8 | 0.01 | akt phosphorylates targets in the cytosol | -2.1 | -2.7 | 0.12 |
| 5 | regulation of the fanconi anemia pathway | +1.1 | +2.8 | 0.01 | pip3 activates akt signaling | -1.9 | -2.5 | 0.03 |

**Erlotinib** - tier C (PFI among treated (Cox)), n = 21, events = 17

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | downstream tcr signaling | +2.9 | +2.9 | 0.33 | +2.4 | iron uptake and transport* | -3.6 | -3.6 | 0.71 | -3.6 |
| 2 | tcr signaling | +2.5 | +2.5 | 0.10 | +2.0 | transferrin endocytosis and recycling* | -3.6 | -3.6 | 0.78 | -3.6 |
| 3 | translocation of zap 70 to immunological synapse | +2.4 | +2.4 | 0.10 | +1.8 | insulin receptor recycling | -3.3 | -3.3 | 0.58 | -3.3 |
| 4 | regulation of kit signaling | +2.4 | +2.4 | 0.11 | +1.8 | alpha linolenic acid ala metabolism | -2.8 | -2.8 | 0.31 | -2.5 |
| 5 | phosphorylation of cd3 and tcr zeta chains | +2.4 | +2.4 | 0.08 | +1.7 | activation of chaperone genes by xbp1s | -2.7 | -2.7 | 0.34 | -2.9 |

**Etoposide** - tier A (RECIST response), n = 15, responders = 12

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | regulation of the fanconi anemia pathway | +2.0 | +3.8 | 0.23 | gap junction degradation | -0.6 | -2.8 | 0.00 |
| 2 | metabolism of rna | +1.1 | +3.1 | 0.00 | o linked glycosylation of mucins | -1.7 | -2.6 | 0.07 |
| 3 | gluconeogenesis | +1.9 | +2.9 | 0.16 | mhc class ii antigen presentation | -0.7 | -2.5 | 0.00 |
| 4 | signaling by fgfr1 mutants | +0.9 | +2.6 | 0.00 | generic transcription pathway | -1.1 | -2.4 | 0.02 |
| 5 | glucose metabolism | +1.4 | +2.6 | 0.07 | rip mediated nfkb activation via dai | -0.9 | -2.4 | 0.00 |

**Gemcitabine** - tier C (PFI among treated (Cox)), n = 25, events = 20

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | notch1 intracellular domain regulates transcription | +1.2 | +2.4 | 0.02 | +1.4 | regulation of water balance by renal aquaporins | -0.7 | -3.1 | 0.01 | -1.3 |
| 2 | cell cell junction organization | +1.8 | +2.1 | 0.17 | +1.5 | activation of kainate receptors upon glutamate binding | -1.6 | -3.1 | 0.01 | -1.3 |
| 3 | trafficking of glur2 containing ampa receptors | +1.0 | +2.0 | 0.02 | +0.8 | downstream signaling events of b cell receptor bcr | -1.0 | -3.1 | 0.01 | -0.2 |
| 4 | sema3a pak dependent axon repulsion | +1.3 | +2.0 | 0.00 | +1.1 | regulation of ornithine decarboxylase odc | -1.3 | -3.0 | 0.04 | -0.3 |
| 5 | regulation of insulin secretion by acetylcholine | +2.8 | +2.0 | 0.42 | +2.4 | pol switching | -1.0 | -2.7 | 0.00 | +0.4 |

**Vinorelbine** - tier A (RECIST response), n = 16, responders = 13

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | glucose metabolism | +0.8 | +3.1 | 0.00 | gap junction degradation | -1.2 | -2.9 | 0.00 |
| 2 | intrinsic pathway for apoptosis | +2.5 | +3.0 | 0.40 | prostanoid ligand receptors | -0.7 | -2.8 | 0.00 |
| 3 | base free sugar phosphate removal via the single nucleo | +1.5 | +3.0 | 0.02 | membrane trafficking | -1.1 | -2.7 | 0.06 |
| 4 | resolution of ap sites via the multiple nucleotide patc | +1.8 | +3.0 | 0.12 | signaling by bmp | -1.3 | -2.6 | 0.00 |
| 5 | metabolism of rna | +0.5 | +2.7 | 0.00 | downstream signaling events of b cell receptor bcr | -0.7 | -2.5 | 0.00 |

## LUSC - lung squamous

**Docetaxel** - tier A (RECIST response), n = 13, responders = 8

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | glucose metabolism | +2.6 | +4.2 | 0.52 | downregulation of smad2 3 smad4 transcriptional activit | -1.7 | -3.0 | 0.15 |
| 2 | creb phosphorylation through the activation of camkii | +1.0 | +3.0 | 0.00 | phospholipase c mediated cascade | -0.8 | -2.9 | 0.01 |
| 3 | transport of organic anions | +1.7 | +2.8 | 0.04 | amino acid synthesis and interconversion transamination | -0.9 | -2.9 | 0.00 |
| 4 | gluconeogenesis | +1.5 | +2.7 | 0.03 | regulation of water balance by renal aquaporins | -0.1 | -2.8 | 0.00 |
| 5 | metabolism of lipids and lipoproteins | +1.0 | +2.6 | 0.01 | downstream signaling events of b cell receptor bcr | -1.1 | -2.8 | 0.12 |

**Gemcitabine** - tier A (RECIST response), n = 17, responders = 9

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | unfolded protein response | +1.2 | +2.8 | 0.01 | regulation of water balance by renal aquaporins | -1.7 | -3.8 | 0.08 |
| 2 | base free sugar phosphate removal via the single nucleo | +0.9 | +2.6 | 0.00 | activation of kainate receptors upon glutamate binding | -1.5 | -3.1 | 0.01 |
| 3 | platelet calcium homeostasis | +0.8 | +2.6 | 0.01 | regulation of ornithine decarboxylase odc | -1.1 | -2.8 | 0.00 |
| 4 | regulation of the fanconi anemia pathway | +1.2 | +2.6 | 0.00 | downregulation of smad2 3 smad4 transcriptional activit | -1.6 | -2.8 | 0.06 |
| 5 | phospholipid metabolism | +1.5 | +2.5 | 0.00 | transferrin endocytosis and recycling | -1.9 | -2.8 | 0.08 |

**Paclitaxel** - tier A (RECIST response), n = 10, responders = 5

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | creb phosphorylation through the activation of camkii | +2.1 | +3.7 | 0.16 | downstream signaling events of b cell receptor bcr | -2.5 | -3.8 | 0.39 |
| 2 | glucose metabolism | +0.1 | +2.6 | 0.00 | regulation of water balance by renal aquaporins | -1.1 | -3.5 | 0.03 |
| 3 | ras activation uopn ca2 infux through nmda receptor | +2.6 | +2.6 | 0.48 | regulation of ornithine decarboxylase odc | -1.9 | -3.4 | 0.01 |
| 4 | post nmda receptor activation events | +2.5 | +2.6 | 0.47 | metabolism of proteins | -1.2 | -3.1 | 0.00 |
| 5 | developmental biology | +1.5 | +2.5 | 0.00 | metabolism of nucleotides | -1.6 | -3.0 | 0.01 |

**Vinorelbine** - tier A (RECIST response), n = 21, responders = 17

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | glucose metabolism | +2.5 | +4.3 | 0.35 | regulation of water balance by renal aquaporins | -1.0 | -3.4 | 0.00 |
| 2 | glycolysis | +2.6 | +3.3 | 0.50 | gap junction degradation | -1.1 | -2.9 | 0.04 |
| 3 | immunoregulatory interactions between a lymphoid and a  | +1.6 | +2.9 | 0.01 | amino acid synthesis and interconversion transamination | -1.1 | -2.8 | 0.00 |
| 4 | metabolism of rna | +0.7 | +2.9 | 0.00 | chondroitin sulfate dermatan sulfate metabolism | -1.0 | -2.8 | 0.00 |
| 5 | il 2 signaling | +1.5 | +2.6 | 0.11 | transferrin endocytosis and recycling | -2.2 | -2.8 | 0.28 |

## OV - ovarian

**Bevacizumab** - tier C (PFI among treated (Cox)), n = 42, events = 32

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | shc mediated signalling | +2.8 | +3.7 | 0.26 | +0.9 | metabolism of proteins | -1.8 | -3.4 | 0.04 | -0.6 |
| 2 | platelet calcium homeostasis | +1.6 | +3.2 | 0.01 | -0.1 | prostanoid ligand receptors | -1.3 | -3.0 | 0.02 | -1.5 |
| 3 | signalling to erks | +2.9 | +3.1 | 0.30 | +1.8 | pyruvate metabolism | -1.6 | -2.8 | 0.03 | +0.5 |
| 4 | signalling to ras | +3.4 | +3.0 | 0.55 | +2.2 | rna pol iii transcription termination | -1.9 | -2.6 | 0.08 | -1.2 |
| 5 | map kinase activation in tlr cascade | +1.2 | +2.8 | 0.00 | +0.8 | formation of atp by chemiosmotic coupling | -2.6 | -2.6 | 0.05 | -0.6 |

**Carboplatin** - tier B (first-course outcome), n = 282, responders = 239

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | cell cell junction organization | +2.8 | +3.3 | 0.39 | prostanoid ligand receptors | -1.2 | -3.1 | 0.02 |
| 2 | ctla4 inhibitory signaling | +2.0 | +3.2 | 0.07 | chondroitin sulfate dermatan sulfate metabolism | -1.7 | -3.1 | 0.06 |
| 3 | cell cell communication | +2.5 | +3.2 | 0.24 | metabolism of proteins | -1.6 | -2.9 | 0.04 |
| 4 | adherens junctions interactions | +2.6 | +3.1 | 0.31 | activation of kainate receptors upon glutamate binding | -1.3 | -2.7 | 0.06 |
| 5 | cell junction organization | +2.6 | +2.7 | 0.29 | signaling by bmp | -1.5 | -2.7 | 0.07 |

**Cisplatin** - tier B (first-course outcome), n = 78, responders = 68

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | creb phosphorylation through the activation of camkii | +1.3 | +2.8 | 0.02 | mtorc1 mediated signalling | -2.2 | -3.0 | 0.25 |
| 2 | transport of organic anions | +0.7 | +2.5 | 0.00 | metabolism of proteins | -1.3 | -2.7 | 0.10 |
| 3 | jnk c jun kinases phosphorylation and activation mediat | +1.2 | +2.4 | 0.01 | sphingolipid de novo biosynthesis | -2.5 | -2.5 | 0.38 |
| 4 | regulation of beta cell development | +1.3 | +2.4 | 0.01 | insulin receptor recycling | -2.0 | -2.4 | 0.11 |
| 5 | cell cell communication | +1.0 | +2.3 | 0.01 | metabolism of vitamins and cofactors | -1.7 | -2.4 | 0.16 |

**Cyclophosphamide** - tier C (PFI among treated (Cox)), n = 20, events = 19

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | shc mediated signalling | +3.2 | +3.7 | 0.55 | +2.0 | pyruvate metabolism | -3.8 | -4.2 | 0.78 | -2.1 |
| 2 | platelet calcium homeostasis | +1.4 | +2.9 | 0.01 | -0.1 | insulin synthesis and processing | -2.0 | -2.6 | 0.10 | -1.7 |
| 3 | shc1 events in erbb4 signaling | +1.3 | +2.2 | 0.01 | +0.7 | regulation of pyruvate dehydrogenase pdh complex | -3.0 | -2.5 | 0.53 | -1.5 |
| 4 | the role of nef in hiv1 replication and disease pathoge | +0.6 | +2.2 | 0.01 | -0.2 | regulation of water balance by renal aquaporins | -0.4 | -2.5 | 0.01 | -0.6 |
| 5 | shc1 events in egfr signaling | +2.0 | +2.1 | 0.06 | +0.7 | pyruvate metabolism and citric acid tca cycle | -2.2 | -2.4 | 0.14 | -1.1 |

**Docetaxel** - tier B (first-course outcome), n = 35, responders = 30

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | phospholipid metabolism | +2.2 | +3.6 | 0.06 | metabolism of proteins | -2.1 | -3.5 | 0.02 |
| 2 | metabolism of lipids and lipoproteins | +2.1 | +3.4 | 0.03 | amino acid synthesis and interconversion transamination | -1.7 | -3.4 | 0.00 |
| 3 | downstream signal transduction | +2.3 | +3.3 | 0.01 | metabolism of nucleotides | -2.0 | -3.2 | 0.08 |
| 4 | glucose metabolism | +1.2 | +3.2 | 0.00 | formation of transcription coupled ner tc ner repair co | -2.4 | -3.0 | 0.01 |
| 5 | alpha linolenic acid ala metabolism | +2.2 | +3.1 | 0.01 | regulation of ornithine decarboxylase odc | -1.3 | -2.9 | 0.00 |

**Doxorubicin** - tier C (PFI among treated (Cox)), n = 115, events = 112

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | platelet calcium homeostasis | +0.6 | +2.6 | 0.00 | -1.2 | gap junction degradation | -1.7 | -3.8 | 0.04 | -1.5 |
| 2 | signaling by nodal | +2.6 | +2.4 | 0.16 | +1.2 | downstream signaling events of b cell receptor bcr | -2.8 | -3.7 | 0.26 | -1.3 |
| 3 | plc beta mediated events | +2.9 | +2.3 | 0.31 | +0.5 | integration of provirus | -2.9 | -3.6 | 0.27 | -0.0 |
| 4 | downstream signal transduction | +0.2 | +2.3 | 0.00 | -1.7 | akt phosphorylates targets in the cytosol | -3.5 | -3.3 | 0.55 | -2.5 |
| 5 | signaling by fgfr1 mutants | +0.3 | +2.2 | 0.00 | -0.7 | mhc class ii antigen presentation | -1.5 | -3.1 | 0.00 | -1.5 |

**Gemcitabine** - tier B (first-course outcome), n = 22, responders = 16

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | synthesis of pips at the early endosome membrane | +1.7 | +2.7 | 0.01 | metabolism of proteins | -2.5 | -3.7 | 0.25 |
| 2 | shc mediated signalling | +2.2 | +2.6 | 0.14 | pyruvate metabolism | -2.2 | -3.2 | 0.08 |
| 3 | synthesis of pips at the late endosome membrane | +1.7 | +2.5 | 0.06 | bile acid and bile salt metabolism | -1.9 | -2.8 | 0.07 |
| 4 | notch1 intracellular domain regulates transcription | +1.0 | +2.3 | 0.02 | regulation of water balance by renal aquaporins | +0.0 | -2.6 | 0.00 |
| 5 | signaling by fgfr1 mutants | +1.2 | +2.3 | 0.00 | prostanoid ligand receptors | -1.6 | -2.6 | 0.10 |

**Paclitaxel** - tier B (first-course outcome), n = 287, responders = 249

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | glucose metabolism | +0.2 | +2.7 | 0.00 | metabolism of proteins | -1.8 | -3.5 | 0.14 |
| 2 | g1 phase | +2.0 | +2.6 | 0.18 | prostanoid ligand receptors | -1.3 | -3.4 | 0.04 |
| 3 | cell cell junction organization | +2.6 | +2.6 | 0.39 | chondroitin sulfate dermatan sulfate metabolism | -1.8 | -3.2 | 0.06 |
| 4 | adherens junctions interactions | +2.4 | +2.6 | 0.25 | regulation of water balance by renal aquaporins | -0.7 | -3.2 | 0.00 |
| 5 | g0 and early g1 | +2.0 | +2.5 | 0.08 | post translational protein modification | -2.1 | -3.2 | 0.15 |

**Tamoxifen** - tier C (PFI among treated (Cox)), n = 27, events = 26

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of rna | -0.0 | +2.4 | 0.01 | +0.0 | neuronal system | -2.3 | -3.4 | 0.26 | -2.4 |
| 2 | conversion from apc c cdc20 to apc c cdh1 in late anaph | +0.9 | +2.4 | 0.00 | -0.1 | tandem pore domain potassium channels | -2.5 | -2.8 | 0.27 | -2.0 |
| 3 | metabolism of lipids and lipoproteins | +0.8 | +2.3 | 0.01 | -0.2 | phospholipase c mediated cascade | -0.5 | -2.5 | 0.00 | +0.3 |
| 4 | steroid hormones | +1.9 | +2.3 | 0.11 | +1.6 | rip mediated nfkb activation via dai | -1.2 | -2.4 | 0.01 | -0.6 |
| 5 | metabolism of steroid hormones and vitamins a and d | +2.3 | +2.3 | 0.23 | +1.8 | platelet aggregation plug formation | -1.2 | -2.4 | 0.00 | -1.6 |

**Topotecan** - tier C (PFI among treated (Cox)), n = 97, events = 97

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | amino acid synthesis and interconversion transamination | +2.7 | +2.7 | 0.37 | +2.0 | tandem pore domain potassium channels | -2.3 | -2.3 | 0.14 | -0.9 |
| 2 | passive transport by aquaporins | +2.6 | +2.6 | 0.38 | +2.0 | il1 signaling | -1.9 | -1.9 | 0.05 | +0.9 |
| 3 | interaction between l1 and ankyrins | +2.5 | +2.5 | 0.28 | -0.2 | prostanoid ligand receptors | -1.9 | -1.9 | 0.10 | -1.5 |
| 4 | ion transport by p type atpases | +2.4 | +2.4 | 0.16 | +1.2 | rna pol iii chain elongation | -1.8 | -1.8 | 0.02 | -0.7 |
| 5 | basigin interactions | +2.2 | +2.2 | 0.19 | +1.7 | traf6 mediated induction of tak1 complex | -1.8 | -1.8 | 0.05 | -1.8 |

## PAAD - pancreatic

**Oxaliplatin** - tier C (PFI among treated (Cox)), n = 21, events = 19

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of lipids and lipoproteins | +2.1 | +3.6 | 0.09 | +2.1 | pol switching | -2.3 | -3.5 | 0.00 | -0.9 |
| 2 | phospholipid metabolism | +2.6 | +3.3 | 0.10 | +2.6 | packaging of telomere ends | -2.3 | -3.4 | 0.01 | -1.8 |
| 3 | signaling by fgfr1 mutants | +1.5 | +3.1 | 0.00 | +0.9 | amino acid synthesis and interconversion transamination | -2.3 | -3.3 | 0.06 | -2.0 |
| 4 | heparan sulfate heparin hs gag metabolism | +2.3 | +3.0 | 0.00 | +2.1 | mitotic g2 g2 m phases | -2.4 | -3.3 | 0.09 | -2.4 |
| 5 | downstream signal transduction | +1.7 | +3.0 | 0.00 | +2.1 | telomere maintenance | -2.5 | -3.3 | 0.03 | -1.7 |

## PRAD - prostate

**Bicalutamide** - tier A (RECIST response), n = 18, responders = 15

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | nuclear events kinase and transcription factor activati | +2.7 | +2.7 | 0.50 | formation of the hiv1 early elongation complex | -2.5 | -2.5 | 0.20 |
| 2 | signaling by erbb2 | +2.6 | +2.6 | 0.31 | viral messenger rna synthesis | -2.5 | -2.5 | 0.17 |
| 3 | erk mapk targets | +2.6 | +2.6 | 0.42 | abortive elongation of hiv1 transcript in the absence o | -2.5 | -2.5 | 0.09 |
| 4 | mapk targets nuclear events mediated by map kinases | +2.6 | +2.6 | 0.38 | formation of transcription coupled ner tc ner repair co | -2.3 | -2.3 | 0.03 |
| 5 | g protein activation | +2.5 | +2.5 | 0.25 | formation of rna pol ii elongation complex  | -2.3 | -2.3 | 0.03 |

**Leuprolide** - tier A (RECIST response), n = 21, responders = 16

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | na cl dependent neurotransmitter transporters | +3.3 | +3.3 | 0.64 | tryptophan catabolism | -3.6 | -3.6 | 0.86 |
| 2 | map kinase activation in tlr cascade | +2.6 | +2.6 | 0.34 | traf6 mediated irf7 activation | -2.9 | -2.9 | 0.37 |
| 3 | prefoldin mediated transfer of substrate to cct tric | +2.6 | +2.6 | 0.30 | synthesis of pc | -2.8 | -2.8 | 0.36 |
| 4 | amino acid transport across the plasma membrane | +2.4 | +2.4 | 0.05 | nef mediates down modulation of cell surface receptors  | -2.6 | -2.6 | 0.18 |
| 5 | formation of tubulin folding intermediates by cct tric | +2.4 | +2.4 | 0.14 | er phagosome pathway | -2.5 | -2.5 | 0.23 |

## READ - rectal

**Oxaliplatin** - tier A (RECIST response), n = 19, responders = 14

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | alpha linolenic acid ala metabolism | +2.2 | +3.2 | 0.17 | gap junction degradation | -1.4 | -3.1 | 0.08 |
| 2 | synthesis of bile acids and bile salts via 24 hydroxych | +2.5 | +3.0 | 0.33 | packaging of telomere ends | -1.5 | -2.9 | 0.03 |
| 3 | signaling by fgfr3 mutants | +0.8 | +2.4 | 0.00 | reversible hydration of carbon dioxide | -1.8 | -2.9 | 0.15 |
| 4 | trafficking of glur2 containing ampa receptors | +1.3 | +2.4 | 0.01 | nef mediated downregulation of mhc class i complex cell | -1.2 | -2.7 | 0.05 |
| 5 | metabolism of lipids and lipoproteins | +0.3 | +2.4 | 0.00 | signaling by bmp | -1.4 | -2.6 | 0.00 |

## SARC - sarcoma

**Ifosfamide** - tier A (RECIST response), n = 16, responders = 9

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | transport of organic anions | +2.4 | +3.2 | 0.09 | regulation of water balance by renal aquaporins | -2.1 | -3.8 | 0.04 |
| 2 | heparan sulfate heparin hs gag metabolism | +2.0 | +2.9 | 0.07 | pol switching | -1.9 | -3.2 | 0.00 |
| 3 | regulated proteolysis of p75ntr | +3.0 | +2.9 | 0.69 | cgmp effects | -2.1 | -3.0 | 0.14 |
| 4 | signaling by fgfr1 mutants | +1.4 | +2.9 | 0.00 | neuronal system | -1.8 | -2.8 | 0.01 |
| 5 | phospholipid metabolism | +1.8 | +2.9 | 0.15 | reversible hydration of carbon dioxide | -2.4 | -2.7 | 0.10 |

## SKCM - melanoma

**Dacarbazine** - tier A (RECIST response), n = 27, responders = 12

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | amine derived hormones | +3.0 | +3.3 | 0.11 | dag and ip3 signaling | -2.9 | -3.6 | 0.09 |
| 2 | common pathway | +2.2 | +2.9 | 0.03 | downregulation of smad2 3 smad4 transcriptional activit | -2.5 | -3.5 | 0.07 |
| 3 | intrinsic pathway for apoptosis | +1.7 | +2.8 | 0.00 | generic transcription pathway | -2.6 | -3.4 | 0.04 |
| 4 | mrna decay by 3 to 5 exoribonuclease | +0.9 | +2.8 | 0.00 | circadian clock* | -3.2 | -3.3 | 0.10 |
| 5 | immunoregulatory interactions between a lymphoid and a  | +1.1 | +2.7 | 0.00 | bmal1 clock npas2 activates circadian expression | -3.0 | -3.2 | 0.10 |

**Interferon Alfa-2B, Recombinant** - tier A (RECIST response), n = 18, responders = 10

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | synthesis of bile acids and bile salts via 7alpha hydro | +2.6 | +2.6 | 0.40 | signaling by hippo | -3.1 | -3.1 | 0.64 |
| 2 | synthesis of bile acids and bile salts | +2.6 | +2.6 | 0.26 | signaling by fgfr1 fusion mutants | -3.0 | -3.0 | 0.48 |
| 3 | synthesis of substrates in n glycan biosythesis | +2.5 | +2.5 | 0.14 | irak2 mediated activation of tak1 complex upon tlr7 8 o | -3.0 | -3.0 | 0.42 |
| 4 | glycolysis | +2.4 | +2.4 | 0.14 | nuclear signaling by erbb4 | -2.7 | -2.7 | 0.35 |
| 5 | respiratory electron transport atp synthesis by chemios | +2.3 | +2.3 | 0.10 | jnk c jun kinases phosphorylation and activation mediat | -2.7 | -2.7 | 0.22 |

**Ipilimumab** - tier A (RECIST response), n = 15, responders = 5

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | yap1 and wwtr1 taz stimulated gene expression | +2.0 | +2.0 | 0.04 | cyclin a b1 associated events during g2 m transition | -2.5 | -2.5 | 0.40 |
| 2 | notch hlh transcription pathway | +2.0 | +2.0 | 0.06 | gap junction trafficking | -2.2 | -2.2 | 0.31 |
| 3 | signaling by notch | +1.9 | +1.9 | 0.18 | cyclin e associated events during g1 s transition  | -2.1 | -2.1 | 0.08 |
| 4 | inwardly rectifying k channels | +1.9 | +1.9 | 0.11 | gap junction assembly | -2.1 | -2.1 | 0.15 |
| 5 | regulation of insulin secretion by acetylcholine | +1.9 | +1.9 | 0.05 | erks are inactivated | -2.1 | -2.1 | 0.11 |

## STAD - stomach

**Etoposide** - tier A (RECIST response), n = 17, responders = 12

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of lipids and lipoproteins | +2.4 | +3.8 | 0.04 | regulation of water balance by renal aquaporins | -1.4 | -3.4 | 0.00 |
| 2 | regulation of the fanconi anemia pathway | +1.3 | +3.3 | 0.01 | activation of kainate receptors upon glutamate binding | -1.8 | -3.1 | 0.07 |
| 3 | adherens junctions interactions | +2.8 | +3.3 | 0.50 | generic transcription pathway | -2.0 | -3.0 | 0.14 |
| 4 | glycosphingolipid metabolism | +2.7 | +3.2 | 0.17 | oxygen dependent proline hydroxylation of hypoxia induc | -2.4 | -3.0 | 0.14 |
| 5 | keratan sulfate keratin metabolism | +2.0 | +3.2 | 0.09 | neuronal system | -1.6 | -2.9 | 0.00 |

## TGCT - testicular germ cell

**Bleomycin** - tier C (PFI among treated (Cox)), n = 53, events = 22

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | platelet calcium homeostasis | +1.1 | +3.1 | 0.01 | +0.9 | metabolism of proteins | -1.2 | -2.8 | 0.01 | -1.5 |
| 2 | ctla4 inhibitory signaling | +1.2 | +2.9 | 0.02 | +0.7 | pol switching | -1.2 | -2.7 | 0.00 | -1.3 |
| 3 | p2y receptors | +1.5 | +2.6 | 0.03 | +1.0 | nef mediated downregulation of mhc class i complex cell | -1.5 | -2.7 | 0.06 | -0.8 |
| 4 | downstream signal transduction | +0.7 | +2.6 | 0.02 | +1.4 | microrna mirna biogenesis | -1.3 | -2.7 | 0.07 | -0.8 |
| 5 | il 3 5 and gm csf signaling | +1.3 | +2.5 | 0.00 | +0.6 | gap junction degradation | -0.4 | -2.5 | 0.00 | +1.0 |

**Cisplatin** - tier C (PFI among treated (Cox)), n = 53, events = 23

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | ctla4 inhibitory signaling | +1.3 | +2.9 | 0.02 | +0.7 | platelet aggregation plug formation | -1.4 | -2.8 | 0.03 | -0.5 |
| 2 | platelet calcium homeostasis | +1.3 | +2.8 | 0.05 | +1.1 | microrna mirna biogenesis | -1.1 | -2.7 | 0.07 | -0.7 |
| 3 | inflammasomes | +1.3 | +2.5 | 0.00 | +0.8 | signaling by bmp | -1.6 | -2.6 | 0.04 | -0.5 |
| 4 | il 3 5 and gm csf signaling | +1.3 | +2.3 | 0.00 | +0.5 | nef mediated downregulation of mhc class i complex cell | -1.2 | -2.6 | 0.07 | -0.7 |
| 5 | il 2 signaling | +1.2 | +2.2 | 0.00 | +0.5 | metabolism of proteins | -1.0 | -2.5 | 0.02 | -1.4 |

**Etoposide** - tier C (PFI among treated (Cox)), n = 54, events = 23

| # | Better-outcome pathway | local z | comb z | stab | pred | Worse-outcome pathway | local z | comb z | stab | pred |
|---|---|---|---|---|---|---|---|---|---|
| 1 | platelet calcium homeostasis | +1.1 | +3.0 | 0.01 | +0.9 | gap junction degradation | -0.5 | -2.8 | 0.00 | +0.9 |
| 2 | downstream signal transduction | +0.4 | +2.5 | 0.00 | +1.2 | microrna mirna biogenesis | -1.0 | -2.6 | 0.04 | -0.6 |
| 3 | ctla4 inhibitory signaling | +0.7 | +2.3 | 0.00 | +0.4 | platelet aggregation plug formation | -1.1 | -2.6 | 0.06 | -0.3 |
| 4 | inflammasomes | +1.2 | +2.3 | 0.00 | +0.7 | pol switching | -0.6 | -2.5 | 0.00 | -1.0 |
| 5 | acyl chain remodelling of ps | +1.8 | +2.2 | 0.16 | +0.3 | metabolism of proteins | -0.8 | -2.4 | 0.00 | -1.3 |

## UCEC - endometrial

**Doxorubicin** - tier A (RECIST response), n = 11, responders = 6

| # | Better-outcome pathway | local z | comb z | stab | Worse-outcome pathway | local z | comb z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | immunoregulatory interactions between a lymphoid and a  | +2.6 | +3.6 | 0.19 | neuronal system | -1.5 | -2.8 | 0.00 |
| 2 | downstream signal transduction | +1.9 | +3.4 | 0.01 | pol switching | -0.9 | -2.7 | 0.00 |
| 3 | the role of nef in hiv1 replication and disease pathoge | +2.5 | +3.4 | 0.08 | microrna mirna biogenesis | -0.8 | -2.4 | 0.00 |
| 4 | il 3 5 and gm csf signaling | +2.4 | +3.3 | 0.07 | generic transcription pathway | -1.1 | -2.4 | 0.00 |
| 5 | synthesis of pips at the early endosome membrane | +2.0 | +3.2 | 0.06 | glutamate neurotransmitter release cycle | -2.9 | -2.4 | 0.64 |

