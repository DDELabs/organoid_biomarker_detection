# Top-10 pathway biomarkers per cancer type and drug (TCGA RECIST response)

Each list = 5 response + 5 resistance pathways ranked by the combined z (local evidence in that cancer, stabilised by the drug's pan-cancer ATLAS effect). `local z` = evidence inside the cancer only; `stab` = fraction of 200 bootstraps in which the pathway is in that cancer's local top 10. Local q < 0.1 is marked *. These are exploratory (20-130 patients per pair).

## BLCA - bladder

**Cisplatin** (n = 64, responders = 41)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | ctla4 inhibitory signaling | +1.9 | +3.3 | 0.03 | microrna mirna biogenesis | -2.9 | -3.9 | 0.27 |
| 2 | signaling by fgfr3 mutants | +1.4 | +3.1 | 0.04 | asparagine n linked glycosylation | -2.7 | -3.2 | 0.23 |
| 3 | immunoregulatory interactions between a lymphoid and a  | +2.2 | +3.1 | 0.04 | metabolism of proteins | -1.9 | -3.2 | 0.07 |
| 4 | il 2 signaling | +2.3 | +3.0 | 0.06 | regulatory rna pathways | -2.8 | -3.1 | 0.24 |
| 5 | platelet calcium homeostasis | +1.6 | +3.0 | 0.01 | downregulation of smad2 3 smad4 transcriptional activit | -1.9 | -3.1 | 0.08 |

**Gemcitabine** (n = 72, responders = 44)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | base free sugar phosphate removal via the single nucleo | +1.4 | +3.0 | 0.02 | gap junction degradation | -2.0 | -3.5 | 0.08 |
| 2 | lipid digestion mobilization and transport | +2.3 | +2.7 | 0.19 | semaphorin interactions | -1.2 | -3.4 | 0.01 |
| 3 | resolution of ap sites via the multiple nucleotide patc | +0.9 | +2.5 | 0.01 | metabolism of proteins | -2.1 | -3.4 | 0.10 |
| 4 | unfolded protein response | +0.7 | +2.4 | 0.00 | asparagine n linked glycosylation | -3.1 | -3.2 | 0.53 |
| 5 | perk regulated gene expression | +1.9 | +2.4 | 0.06 | regulation of water balance by renal aquaporins | -0.8 | -3.2 | 0.00 |

## BRCA - breast

**Cyclophosphamide** (n = 161, responders = 152)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | regulation of the fanconi anemia pathway | +2.6 | +3.9 | 0.23 | metabolism of proteins | -1.6 | -3.3 | 0.01 |
| 2 | metabolism of rna | +2.2 | +3.6 | 0.01 | ionotropic activity of kainate receptors | -2.1 | -3.2 | 0.07 |
| 3 | shc mediated signalling | +2.2 | +3.1 | 0.15 | prostanoid ligand receptors | -2.0 | -3.0 | 0.12 |
| 4 | platelet calcium homeostasis | +1.5 | +3.0 | 0.00 | activation of kainate receptors upon glutamate binding | -1.0 | -3.0 | 0.01 |
| 5 | ctla4 inhibitory signaling | +1.5 | +2.9 | 0.00 | gap junction degradation | -1.2 | -2.9 | 0.00 |

**Docetaxel** (n = 73, responders = 67)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | synthesis of bile acids and bile salts via 24 hydroxych | +2.1 | +3.3 | 0.20 | regulation of water balance by renal aquaporins | -1.5 | -3.8 | 0.01 |
| 2 | amine derived hormones | +1.8 | +2.9 | 0.05 | activation of kainate receptors upon glutamate binding | -2.1 | -3.4 | 0.14 |
| 3 | alpha linolenic acid ala metabolism | +1.8 | +2.9 | 0.19 | phospholipase c mediated cascade | -1.3 | -3.2 | 0.00 |
| 4 | glucose metabolism | +0.4 | +2.7 | 0.00 | membrane trafficking | -1.4 | -2.9 | 0.09 |
| 5 | creb phosphorylation through the activation of camkii | +0.4 | +2.6 | 0.00 | gap junction degradation | -1.2 | -2.8 | 0.06 |

**Doxorubicin** (n = 96, responders = 87)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | regulation of the fanconi anemia pathway | +2.2 | +3.8 | 0.14 | gap junction degradation | -1.3 | -3.5 | 0.01 |
| 2 | metabolism of rna | +1.9 | +3.5 | 0.00 | prostanoid ligand receptors | -1.7 | -3.0 | 0.06 |
| 3 | shc mediated signalling | +2.3 | +3.1 | 0.34 | regulation of water balance by renal aquaporins | -0.8 | -2.8 | 0.00 |
| 4 | platelet calcium homeostasis | +1.2 | +3.0 | 0.00 | activation of kainate receptors upon glutamate binding | -1.7 | -2.8 | 0.02 |
| 5 | metabolism of non coding rna | +2.2 | +3.0 | 0.04 | acyl chain remodelling of pc | -1.1 | -2.8 | 0.00 |

**Paclitaxel** (n = 64, responders = 55)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of rna | +1.8 | +3.4 | 0.02 | metabolism of proteins | -2.0 | -3.7 | 0.17 |
| 2 | regulation of the fanconi anemia pathway | +1.7 | +3.1 | 0.01 | prostanoid ligand receptors | -1.7 | -3.7 | 0.01 |
| 3 | g1 phase | +2.5 | +3.0 | 0.45 | gap junction degradation | -2.1 | -3.5 | 0.05 |
| 4 | glucose metabolism | +0.5 | +2.9 | 0.00 | regulation of water balance by renal aquaporins | -0.9 | -3.4 | 0.01 |
| 5 | ctla4 inhibitory signaling | +1.4 | +2.9 | 0.01 | chondroitin sulfate dermatan sulfate metabolism | -1.6 | -3.1 | 0.02 |

**Tamoxifen** (n = 25, responders = 18)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of rna | +1.2 | +3.2 | 0.00 | neuronal system | -1.3 | -2.6 | 0.01 |
| 2 | common pathway | +1.9 | +2.9 | 0.23 | phospholipase c mediated cascade | -0.7 | -2.6 | 0.03 |
| 3 | glucose metabolism | +0.9 | +2.8 | 0.01 | acyl chain remodelling of pc | -0.6 | -2.5 | 0.03 |
| 4 | base free sugar phosphate removal via the single nucleo | +1.1 | +2.7 | 0.07 | effects of pip2 hydrolysis | -0.9 | -2.5 | 0.02 |
| 5 | conversion from apc c cdc20 to apc c cdh1 in late anaph | +1.2 | +2.6 | 0.01 | amino acid synthesis and interconversion transamination | -0.9 | -2.4 | 0.03 |

## CESC - cervical

**Cisplatin** (n = 65, responders = 56)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | cell cycle | +2.4 | +4.0 | 0.07 | regulation of water balance by renal aquaporins | -1.8 | -3.5 | 0.03 |
| 2 | metabolism of rna | +1.9 | +3.6 | 0.01 | rig i mda5 mediated induction of ifn alpha beta pathway | -1.1 | -2.9 | 0.00 |
| 3 | mrna decay by 3 to 5 exoribonuclease | +2.6 | +3.5 | 0.27 | rip mediated nfkb activation via dai | -1.2 | -2.9 | 0.00 |
| 4 | metabolism of amino acids and derivatives | +2.8 | +3.2 | 0.38 | semaphorin interactions | -1.6 | -2.8 | 0.03 |
| 5 | signaling by fgfr3 mutants | +1.3 | +3.0 | 0.02 | circadian clock | -1.8 | -2.8 | 0.01 |

## COAD - colon

**Bevacizumab** (n = 20, responders = 7)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | downstream signal transduction | +1.5 | +3.1 | 0.02 | metabolism of nucleotides | -1.8 | -2.9 | 0.04 |
| 2 | jnk c jun kinases phosphorylation and activation mediat | +1.9 | +3.0 | 0.09 | pol switching | -1.2 | -2.9 | 0.01 |
| 3 | glucose metabolism | +1.8 | +3.0 | 0.10 | regulation of ornithine decarboxylase odc | -1.0 | -2.8 | 0.01 |
| 4 | platelet calcium homeostasis | +1.2 | +2.9 | 0.01 | microrna mirna biogenesis | -1.5 | -2.8 | 0.07 |
| 5 | map kinase activation in tlr cascade | +1.2 | +2.8 | 0.02 | gap junction degradation | -0.7 | -2.6 | 0.00 |

**Fluorouracil** (n = 64, responders = 43)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | glucose metabolism | +2.4 | +3.9 | 0.26 | gap junction degradation | -1.7 | -3.5 | 0.12 |
| 2 | alpha linolenic acid ala metabolism | +1.9 | +3.3 | 0.17 | integration of provirus | -2.0 | -3.2 | 0.15 |
| 3 | gluconeogenesis | +2.2 | +3.1 | 0.18 | activation of kainate receptors upon glutamate binding | -1.8 | -3.2 | 0.07 |
| 4 | metabolism of polyamines | +2.2 | +3.0 | 0.29 | other semaphorin interactions | -1.2 | -2.9 | 0.01 |
| 5 | metabolism of lipids and lipoproteins | +1.6 | +2.9 | 0.04 | regulation of water balance by renal aquaporins | -1.1 | -2.7 | 0.01 |

**Oxaliplatin** (n = 51, responders = 37)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | glucose metabolism | +2.1 | +3.4 | 0.16 | gap junction degradation | -2.4 | -3.7 | 0.29 |
| 2 | jnk c jun kinases phosphorylation and activation mediat | +2.2 | +3.3 | 0.14 | insulin receptor recycling | -2.2 | -2.9 | 0.23 |
| 3 | transport of organic anions | +2.0 | +3.3 | 0.09 | prostanoid ligand receptors | -0.9 | -2.9 | 0.01 |
| 4 | metabolism of lipids and lipoproteins | +1.3 | +3.0 | 0.02 | il 7 signaling | -2.2 | -2.8 | 0.17 |
| 5 | metabolism of polyamines | +2.5 | +2.9 | 0.47 | semaphorin interactions | -1.2 | -2.8 | 0.00 |

## ESCA - oesophageal

**Capecitabine** (n = 22, responders = 17)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | alpha linolenic acid ala metabolism | +2.4 | +3.8 | 0.29 | metabolism of proteins | -1.5 | -3.3 | 0.02 |
| 2 | heparan sulfate heparin hs gag metabolism | +2.0 | +2.8 | 0.09 | integration of provirus | -1.9 | -3.2 | 0.09 |
| 3 | muscle contraction | +2.0 | +2.6 | 0.21 | mtorc1 mediated signalling | -2.1 | -3.2 | 0.18 |
| 4 | signaling by fgfr1 mutants | +0.7 | +2.5 | 0.00 | packaging of telomere ends | -1.7 | -3.0 | 0.01 |
| 5 | metabolism of carbohydrates | +1.8 | +2.5 | 0.04 | telomere maintenance | -1.5 | -2.9 | 0.01 |

## HNSC - head & neck

**Carboplatin** (n = 28, responders = 22)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | creb phosphorylation through the activation of camkii | +0.9 | +3.1 | 0.00 | metabolism of nucleotides | -2.3 | -3.4 | 0.14 |
| 2 | downstream signal transduction | +1.9 | +3.0 | 0.03 | metabolism of proteins | -2.0 | -3.2 | 0.06 |
| 3 | platelet calcium homeostasis | +2.0 | +3.0 | 0.03 | packaging of telomere ends | -2.0 | -3.1 | 0.01 |
| 4 | alpha linolenic acid ala metabolism | +2.1 | +2.9 | 0.09 | telomere maintenance | -2.0 | -3.0 | 0.03 |
| 5 | developmental biology | +1.8 | +2.8 | 0.02 | regulation of ornithine decarboxylase odc | -1.7 | -3.0 | 0.03 |

**Cisplatin** (n = 50, responders = 45)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | signaling by fgfr3 mutants | +1.5 | +3.1 | 0.01 | rig i mda5 mediated induction of ifn alpha beta pathway | -1.2 | -3.0 | 0.00 |
| 2 | keratan sulfate biosynthesis | +1.6 | +2.7 | 0.04 | downstream signaling events of b cell receptor bcr | -1.5 | -3.0 | 0.01 |
| 3 | platelet homeostasis | +1.7 | +2.6 | 0.13 | microrna mirna biogenesis | -1.6 | -3.0 | 0.04 |
| 4 | keratan sulfate keratin metabolism | +1.6 | +2.6 | 0.07 | oxygen dependent proline hydroxylation of hypoxia induc | -1.8 | -2.8 | 0.12 |
| 5 | transport of organic anions | +0.8 | +2.6 | 0.12 | circadian clock | -1.6 | -2.7 | 0.07 |

**Paclitaxel** (n = 22, responders = 16)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | signaling by fgfr1 mutants | +1.8 | +3.4 | 0.01 | regulation of ornithine decarboxylase odc | -1.9 | -3.4 | 0.01 |
| 2 | shc mediated signalling | +2.2 | +3.2 | 0.19 | metabolism of proteins | -1.6 | -3.4 | 0.01 |
| 3 | downstream signal transduction | +2.0 | +3.1 | 0.09 | metabolism of nucleotides | -1.7 | -3.1 | 0.03 |
| 4 | phospholipid metabolism | +1.5 | +3.0 | 0.01 | pyruvate metabolism | -1.9 | -3.0 | 0.04 |
| 5 | alpha linolenic acid ala metabolism | +2.0 | +2.9 | 0.06 | downstream signaling events of b cell receptor bcr | -1.3 | -2.9 | 0.00 |

## LGG - lower-grade glioma

**Temozolomide** (n = 131, responders = 21)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | mrna decay by 3 to 5 exoribonuclease | +1.7 | +3.5 | 0.13 | activation of kainate receptors upon glutamate binding | -0.9 | -3.0 | 0.01 |
| 2 | developmental biology | +1.4 | +2.9 | 0.03 | acyl chain remodelling of pc | -1.2 | -2.8 | 0.03 |
| 3 | platelet calcium homeostasis | +1.3 | +2.8 | 0.03 | downregulation of smad2 3 smad4 transcriptional activit | -1.4 | -2.7 | 0.09 |
| 4 | map kinase activation in tlr cascade | +1.4 | +2.8 | 0.04 | packaging of telomere ends | -1.0 | -2.6 | 0.04 |
| 5 | shc mediated signalling | +1.4 | +2.5 | 0.05 | branched chain amino acid catabolism | -1.5 | -2.5 | 0.14 |

## LUAD - lung adeno

**Carboplatin** (n = 43, responders = 24)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of lipids and lipoproteins | +1.8 | +3.4 | 0.02 | fgfr4 ligand binding and activation | -2.0 | -3.1 | 0.13 |
| 2 | creb phosphorylation through the activation of camkii | +1.2 | +3.3 | 0.01 | downregulation of smad2 3 smad4 transcriptional activit | -1.4 | -2.9 | 0.04 |
| 3 | alpha linolenic acid ala metabolism | +2.6 | +3.2 | 0.41 | signaling by bmp | -1.8 | -2.9 | 0.03 |
| 4 | glucose metabolism | +1.9 | +3.1 | 0.07 | fgfr ligand binding and activation | -2.0 | -2.8 | 0.07 |
| 5 | platelet calcium homeostasis | +2.0 | +3.0 | 0.07 | packaging of telomere ends | -1.4 | -2.7 | 0.01 |

**Cisplatin** (n = 51, responders = 40)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | platelet homeostasis | +2.4 | +3.1 | 0.21 | generic transcription pathway | -1.9 | -3.1 | 0.12 |
| 2 | regulation of the fanconi anemia pathway | +1.8 | +3.0 | 0.20 | rig i mda5 mediated induction of ifn alpha beta pathway | -1.1 | -2.9 | 0.00 |
| 3 | ctla4 inhibitory signaling | +1.4 | +3.0 | 0.02 | nef mediated downregulation of mhc class i complex cell | -1.4 | -2.8 | 0.01 |
| 4 | keratan sulfate biosynthesis | +1.8 | +2.9 | 0.07 | traf6 mediated irf7 activation | -1.8 | -2.7 | 0.10 |
| 5 | transport of organic anions | +1.2 | +2.9 | 0.03 | rip mediated nfkb activation via dai | -0.6 | -2.5 | 0.00 |

**Paclitaxel** (n = 27, responders = 14)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | glucose metabolism | +1.2 | +3.4 | 0.03 | metabolism of proteins | -1.7 | -3.5 | 0.10 |
| 2 | creb phosphorylation through the activation of camkii | +1.6 | +3.4 | 0.01 | prostanoid ligand receptors | -1.2 | -3.3 | 0.02 |
| 3 | signaling by fgfr1 mutants | +1.6 | +3.3 | 0.07 | regulation of water balance by renal aquaporins | -0.9 | -3.3 | 0.01 |
| 4 | sema3a pak dependent axon repulsion | +2.2 | +3.2 | 0.33 | cgmp effects | -2.3 | -3.0 | 0.34 |
| 5 | platelet calcium homeostasis | +1.9 | +3.2 | 0.03 | signaling by bmp | -1.5 | -2.9 | 0.01 |

**Pemetrexed** (n = 28, responders = 16)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | phospholipid metabolism | +1.4 | +2.8 | 0.01 | downstream signaling events of b cell receptor bcr | -0.7 | -2.8 | 0.01 |
| 2 | sphingolipid metabolism | +1.7 | +2.6 | 0.06 | mrna 3 end processing | -1.3 | -2.7 | 0.01 |
| 3 | metabolism of lipids and lipoproteins | +1.4 | +2.6 | 0.02 | regulation of ornithine decarboxylase odc | -0.6 | -2.6 | 0.00 |
| 4 | unfolded protein response | +1.2 | +2.5 | 0.00 | fgfr4 ligand binding and activation | -1.7 | -2.6 | 0.07 |
| 5 | downstream signal transduction | +1.3 | +2.4 | 0.01 | effects of pip2 hydrolysis | -0.8 | -2.5 | 0.03 |

## LUSC - lung squamous

**Carboplatin** (n = 28, responders = 20)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of rna | +1.0 | +2.9 | 0.01 | phospholipase c mediated cascade | -2.2 | -3.7 | 0.18 |
| 2 | creb phosphorylation through the activation of camkii | +0.3 | +2.7 | 0.00 | regulation of water balance by renal aquaporins | -1.9 | -3.6 | 0.01 |
| 3 | mrna decay by 3 to 5 exoribonuclease | +1.4 | +2.6 | 0.04 | prostanoid ligand receptors | -1.4 | -3.3 | 0.01 |
| 4 | base free sugar phosphate removal via the single nucleo | +1.3 | +2.5 | 0.01 | downregulation of smad2 3 smad4 transcriptional activit | -1.8 | -3.2 | 0.00 |
| 5 | glucuronidation | +1.7 | +2.4 | 0.09 | downstream signaling events of b cell receptor bcr | -1.4 | -3.2 | 0.06 |

**Cisplatin** (n = 32, responders = 27)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of rna | +1.9 | +3.6 | 0.05 | regulation of water balance by renal aquaporins | -1.6 | -3.4 | 0.01 |
| 2 | cell cycle | +1.6 | +3.4 | 0.04 | prostanoid ligand receptors | -1.5 | -3.1 | 0.01 |
| 3 | base free sugar phosphate removal via the single nucleo | +1.4 | +2.7 | 0.00 | acyl chain remodelling of pc | -1.4 | -3.0 | 0.04 |
| 4 | metabolism of mrna | +1.8 | +2.6 | 0.06 | nef mediated downregulation of mhc class i complex cell | -1.6 | -2.9 | 0.01 |
| 5 | regulation of the fanconi anemia pathway | +1.2 | +2.6 | 0.02 | gap junction degradation | -1.3 | -2.8 | 0.01 |

## MESO - mesothelioma

**Cisplatin** (n = 21, responders = 8)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | base free sugar phosphate removal via the single nucleo | +1.6 | +2.9 | 0.00 | prostanoid ligand receptors | -1.9 | -3.4 | 0.07 |
| 2 | cell cycle | +0.7 | +2.8 | 0.00 | acyl chain remodelling of pc | -1.9 | -3.3 | 0.07 |
| 3 | resolution of ap sites via the multiple nucleotide patc | +2.0 | +2.8 | 0.14 | phospholipase c mediated cascade | -1.9 | -3.2 | 0.07 |
| 4 | mitochondrial protein import | +1.8 | +2.5 | 0.01 | rig i mda5 mediated induction of ifn alpha beta pathway | -1.4 | -3.2 | 0.02 |
| 5 | transport of glucose and other sugars bile salts and or | +1.1 | +2.4 | 0.02 | platelet aggregation plug formation | -1.7 | -3.0 | 0.04 |

**Pemetrexed** (n = 24, responders = 8)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | base free sugar phosphate removal via the single nucleo | +1.8 | +3.3 | 0.01 | regulation of water balance by renal aquaporins | -2.0 | -3.9 | 0.06 |
| 2 | resolution of ap sites via the multiple nucleotide patc | +1.5 | +3.0 | 0.01 | phospholipase c mediated cascade | -1.7 | -3.4 | 0.04 |
| 3 | glucose metabolism | +1.6 | +2.8 | 0.06 | activation of kainate receptors upon glutamate binding | -1.7 | -3.3 | 0.08 |
| 4 | resolution of ap sites via the single nucleotide replac | +1.7 | +2.7 | 0.00 | mhc class ii antigen presentation | -1.7 | -3.3 | 0.03 |
| 5 | base excision repair | +1.7 | +2.7 | 0.01 | semaphorin interactions | -1.2 | -3.2 | 0.00 |

## PAAD - pancreatic

**Fluorouracil** (n = 20, responders = 9)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | ctla4 inhibitory signaling | +1.5 | +3.2 | 0.01 | metabolism of proteins | -1.6 | -3.3 | 0.12 |
| 2 | downstream signal transduction | +1.7 | +3.1 | 0.01 | amino acid synthesis and interconversion transamination | -1.8 | -3.0 | 0.13 |
| 3 | signaling by fgfr1 mutants | +1.4 | +3.0 | 0.02 | reversible hydration of carbon dioxide | -1.7 | -3.0 | 0.07 |
| 4 | notch1 intracellular domain regulates transcription | +2.1 | +2.8 | 0.23 | telomere maintenance | -1.9 | -3.0 | 0.01 |
| 5 | alpha linolenic acid ala metabolism | +1.0 | +2.6 | 0.01 | pol switching | -1.2 | -2.9 | 0.00 |

**Gemcitabine** (n = 66, responders = 29)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | inflammasomes | +1.1 | +2.5 | 0.01 | regulation of water balance by renal aquaporins | -1.3 | -3.6 | 0.04 |
| 2 | mrna decay by 3 to 5 exoribonuclease | +1.0 | +2.4 | 0.01 | downstream signaling events of b cell receptor bcr | -1.5 | -3.4 | 0.07 |
| 3 | signaling by fgfr3 mutants | +1.2 | +2.4 | 0.01 | semaphorin interactions | -1.0 | -3.2 | 0.01 |
| 4 | metabolism of rna | +0.1 | +2.3 | 0.00 | erks are inactivated | -2.3 | -3.1 | 0.30 |
| 5 | base free sugar phosphate removal via the single nucleo | +0.5 | +2.3 | 0.00 | pyruvate metabolism | -2.0 | -3.1 | 0.19 |

## READ - rectal

**Fluorouracil** (n = 33, responders = 26)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | immunoregulatory interactions between a lymphoid and a  | +2.2 | +3.4 | 0.11 | pol switching | -1.8 | -3.3 | 0.01 |
| 2 | regulation of ifng signaling | +2.3 | +3.1 | 0.23 | telomere maintenance | -2.0 | -3.1 | 0.03 |
| 3 | il 3 5 and gm csf signaling | +2.1 | +2.9 | 0.10 | metabolism of proteins | -1.1 | -3.0 | 0.00 |
| 4 | costimulation by the cd28 family | +2.0 | +2.9 | 0.06 | amino acid synthesis and interconversion transamination | -1.7 | -2.9 | 0.01 |
| 5 | il 2 signaling | +1.9 | +2.8 | 0.10 | mtorc1 mediated signalling | -1.9 | -2.9 | 0.06 |

## SARC - sarcoma

**Docetaxel** (n = 31, responders = 14)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of rna | +1.4 | +3.3 | 0.09 | regulation of water balance by renal aquaporins | -1.3 | -3.7 | 0.01 |
| 2 | regulation of the fanconi anemia pathway | +1.8 | +3.3 | 0.10 | neuronal system | -2.0 | -3.1 | 0.14 |
| 3 | map kinase activation in tlr cascade | +1.6 | +3.0 | 0.01 | activation of kainate receptors upon glutamate binding | -1.2 | -2.8 | 0.00 |
| 4 | keratan sulfate biosynthesis | +1.6 | +3.0 | 0.05 | phospholipase c mediated cascade | -0.6 | -2.7 | 0.00 |
| 5 | keratan sulfate keratin metabolism | +1.6 | +2.9 | 0.10 | pyruvate metabolism | -1.2 | -2.6 | 0.07 |

**Doxorubicin** (n = 31, responders = 11)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | base free sugar phosphate removal via the single nucleo | +1.9 | +3.3 | 0.03 | regulation of water balance by renal aquaporins | -2.4 | -3.9 | 0.14 |
| 2 | metabolism of lipids and lipoproteins | +1.4 | +3.0 | 0.01 | phospholipase c mediated cascade | -2.5 | -3.5 | 0.34 |
| 3 | sphingolipid metabolism | +2.1 | +3.0 | 0.11 | effects of pip2 hydrolysis | -1.8 | -3.3 | 0.03 |
| 4 | glycosphingolipid metabolism | +2.3 | +3.0 | 0.12 | platelet aggregation plug formation | -1.7 | -3.0 | 0.00 |
| 5 | heparan sulfate heparin hs gag metabolism | +1.7 | +2.9 | 0.01 | neuronal system | -1.9 | -3.0 | 0.02 |

**Gemcitabine** (n = 34, responders = 13)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | metabolism of rna | +1.4 | +3.3 | 0.09 | regulation of water balance by renal aquaporins | -1.2 | -3.5 | 0.00 |
| 2 | regulation of the fanconi anemia pathway | +1.8 | +3.0 | 0.07 | semaphorin interactions | -0.7 | -3.1 | 0.00 |
| 3 | unfolded protein response | +1.3 | +2.9 | 0.01 | neuronal system | -2.0 | -2.9 | 0.14 |
| 4 | g1 phase | +2.0 | +2.9 | 0.15 | activation of kainate receptors upon glutamate binding | -1.2 | -2.9 | 0.01 |
| 5 | base free sugar phosphate removal via the single nucleo | +0.8 | +2.6 | 0.01 | recycling of bile acids and salts | -1.7 | -2.6 | 0.04 |

## STAD - stomach

**Capecitabine** (n = 32, responders = 22)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | glucose metabolism | +0.7 | +2.8 | 0.01 | activation of kainate receptors upon glutamate binding | -1.8 | -3.2 | 0.12 |
| 2 | adherens junctions interactions | +1.5 | +2.6 | 0.08 | platelet aggregation plug formation | -1.7 | -2.9 | 0.07 |
| 3 | alpha linolenic acid ala metabolism | +0.5 | +2.5 | 0.00 | nef mediated downregulation of mhc class i complex cell | -1.8 | -2.8 | 0.10 |
| 4 | signaling by fgfr3 mutants | +1.3 | +2.4 | 0.01 | regulation of water balance by renal aquaporins | -1.2 | -2.7 | 0.03 |
| 5 | base free sugar phosphate removal via the single nucleo | +1.1 | +2.4 | 0.00 | downstream signaling events of b cell receptor bcr | -0.8 | -2.7 | 0.01 |

**Cisplatin** (n = 36, responders = 21)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | cell cycle | +2.1 | +3.8 | 0.01 | regulation of water balance by renal aquaporins | -1.6 | -3.4 | 0.01 |
| 2 | regulation of the fanconi anemia pathway | +2.1 | +3.2 | 0.18 | prostanoid ligand receptors | -1.2 | -2.9 | 0.00 |
| 3 | transport of organic anions | +1.6 | +3.1 | 0.04 | acyl chain remodelling of pc | -1.1 | -2.8 | 0.00 |
| 4 | metabolism of rna | +1.2 | +3.1 | 0.00 | rig i mda5 mediated induction of ifn alpha beta pathway | -0.8 | -2.8 | 0.02 |
| 5 | g0 and early g1 | +2.4 | +2.9 | 0.16 | activation of kainate receptors upon glutamate binding | -1.2 | -2.5 | 0.01 |

**Epirubicin** (n = 25, responders = 18)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | regulation of the fanconi anemia pathway | +2.1 | +3.8 | 0.24 | gap junction degradation | -1.6 | -3.7 | 0.04 |
| 2 | ctla4 inhibitory signaling | +1.5 | +3.0 | 0.11 | regulation of water balance by renal aquaporins | -0.7 | -2.8 | 0.00 |
| 3 | base free sugar phosphate removal via the single nucleo | +1.3 | +2.9 | 0.04 | acyl chain remodelling of pc | -0.9 | -2.7 | 0.00 |
| 4 | signaling by fgfr1 mutants | +1.2 | +2.9 | 0.05 | platelet aggregation plug formation | -1.2 | -2.7 | 0.03 |
| 5 | resolution of ap sites via the multiple nucleotide patc | +1.7 | +2.8 | 0.09 | nef mediated downregulation of mhc class i complex cell | -1.5 | -2.5 | 0.03 |

**Fluorouracil** (n = 88, responders = 58)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | mrna decay by 3 to 5 exoribonuclease | +1.9 | +3.4 | 0.04 | gap junction degradation | -2.1 | -3.8 | 0.18 |
| 2 | metabolism of polyamines | +2.6 | +3.4 | 0.40 | circadian clock | -2.3 | -3.3 | 0.28 |
| 3 | ctla4 inhibitory signaling | +1.7 | +3.3 | 0.08 | other semaphorin interactions | -1.4 | -2.9 | 0.00 |
| 4 | metabolism of rna | +1.3 | +3.0 | 0.00 | chondroitin sulfate dermatan sulfate metabolism | -1.0 | -2.9 | 0.00 |
| 5 | presynaptic nicotinic acetylcholine receptors | +2.3 | +2.9 | 0.30 | signaling by bmp | -1.7 | -2.8 | 0.04 |

**Oxaliplatin** (n = 21, responders = 13)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | creb phosphorylation through the activation of camkii | +1.1 | +3.0 | 0.12 | mitotic g2 g2 m phases | -2.2 | -3.2 | 0.32 |
| 2 | jnk c jun kinases phosphorylation and activation mediat | +1.4 | +2.8 | 0.01 | microrna mirna biogenesis | -1.6 | -3.0 | 0.07 |
| 3 | signaling by fgfr1 mutants | +0.9 | +2.7 | 0.01 | mhc class ii antigen presentation | -1.4 | -2.9 | 0.08 |
| 4 | signaling by fgfr3 mutants | +1.1 | +2.6 | 0.01 | pol switching | -1.1 | -2.8 | 0.01 |
| 5 | ctla4 inhibitory signaling | +0.7 | +2.4 | 0.00 | platelet aggregation plug formation | -1.5 | -2.8 | 0.01 |

## UCEC - endometrial

**Carboplatin** (n = 54, responders = 46)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | creb phosphorylation through the activation of camkii | +0.3 | +2.7 | 0.00 | prostanoid ligand receptors | -2.1 | -3.8 | 0.29 |
| 2 | base free sugar phosphate removal via the single nucleo | +1.5 | +2.6 | 0.09 | amino acid synthesis and interconversion transamination | -2.4 | -3.5 | 0.38 |
| 3 | signaling by fgfr1 mutants | +1.1 | +2.5 | 0.00 | chondroitin sulfate dermatan sulfate metabolism | -1.9 | -3.2 | 0.17 |
| 4 | conversion from apc c cdc20 to apc c cdh1 in late anaph | +1.4 | +2.5 | 0.03 | regulation of water balance by renal aquaporins | -1.2 | -3.1 | 0.04 |
| 5 | sema3a pak dependent axon repulsion | +1.4 | +2.5 | 0.06 | activated notch1 transmits signal to the nucleus | -2.1 | -2.7 | 0.23 |

**Paclitaxel** (n = 52, responders = 44)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | signaling by fgfr1 mutants | +1.4 | +3.1 | 0.06 | prostanoid ligand receptors | -2.5 | -4.2 | 0.39 |
| 2 | base free sugar phosphate removal via the single nucleo | +1.6 | +2.9 | 0.14 | amino acid synthesis and interconversion transamination | -2.3 | -3.8 | 0.32 |
| 3 | conversion from apc c cdc20 to apc c cdh1 in late anaph | +1.7 | +2.7 | 0.01 | regulation of water balance by renal aquaporins | -1.2 | -3.5 | 0.02 |
| 4 | sema3a pak dependent axon repulsion | +1.3 | +2.6 | 0.06 | chondroitin sulfate dermatan sulfate metabolism | -1.9 | -3.3 | 0.15 |
| 5 | signaling by fgfr1 fusion mutants | +1.5 | +2.5 | 0.06 | activated notch1 transmits signal to the nucleus | -2.2 | -2.9 | 0.26 |

## UCS - uterine carcinosarcoma

**Carboplatin** (n = 29, responders = 19)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | trafficking of glur2 containing ampa receptors | +2.0 | +2.9 | 0.23 | pol switching | -1.3 | -3.0 | 0.06 |
| 2 | role of dcc in regulating apoptosis | +2.4 | +2.9 | 0.33 | pyruvate metabolism | -2.0 | -2.9 | 0.17 |
| 3 | metabolism of lipids and lipoproteins | +0.8 | +2.8 | 0.00 | microrna mirna biogenesis | -1.6 | -2.8 | 0.04 |
| 4 | creb phosphorylation through the activation of camkii | +0.4 | +2.7 | 0.00 | dscam interactions | -2.1 | -2.8 | 0.26 |
| 5 | p2y receptors | +1.7 | +2.7 | 0.01 | downstream signaling events of b cell receptor bcr | -0.9 | -2.8 | 0.01 |

**Paclitaxel** (n = 30, responders = 22)

| # | Response pathway | local z | combined z | stab | Resistance pathway | local z | combined z | stab |
|---|---|---|---|---|---|---|---|---|
| 1 | il 3 5 and gm csf signaling | +1.7 | +3.0 | 0.01 | prostanoid ligand receptors | -1.1 | -3.2 | 0.03 |
| 2 | il 2 signaling | +1.7 | +2.8 | 0.01 | pyruvate metabolism | -2.2 | -3.2 | 0.25 |
| 3 | ctla4 inhibitory signaling | +1.1 | +2.7 | 0.01 | metabolism of proteins | -1.1 | -3.0 | 0.00 |
| 4 | immunoregulatory interactions between a lymphoid and a  | +1.2 | +2.6 | 0.00 | amino acid synthesis and interconversion transamination | -1.2 | -3.0 | 0.03 |
| 5 | shc mediated signalling | +1.2 | +2.5 | 0.02 | formation of rna pol ii elongation complex  | -2.0 | -2.9 | 0.10 |

