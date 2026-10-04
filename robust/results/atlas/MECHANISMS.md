# ATLAS drug-response pathways: literature annotation

**Scope.** Annotation of pathway–response associations from the ATLAS hierarchical model (TCGA pan-cancer, 2,769 patient–drug RECIST labels, 1,388 patients, 29 cancers, 21 drugs in 10 classes). Effect = log-odds of response per SD of pathway activity (within-cancer centred, adjusted for treatment setting and proliferation). **Response** means higher activity goes with response; **resistance** means higher activity goes with non-response. z is bootstrap-based. q is BH-FDR within each layer (Tier 1) or across all drug × pathway totals (Tier 2).

**Classification.**
- **KNOWN**: there is direct prior evidence for this pathway or its genes with this drug or class, in the same direction.
- **PLAUSIBLE**: the mechanism fits, but the evidence is indirect or the direction is only partly supported.
- **NOVEL**: we found no prior evidence.
- **ARTEFACT?**: likely reflects tumour or tissue composition (neuronal, renal and similar pathways).
- **CONTRA**: prior evidence points the other way.

References were checked against PubMed (esearch/esummary) unless marked *(unverified)*.

**Caveats.** (i) No drug-total effect reaches q < 0.1. All Tier 2 items have q ≥ 0.32 and are hypotheses only. (ii) Bevacizumab is the only drug in the antiangiogenic class, so its class and drug layers are the same signal, not two confirmations. Its labels come largely from COAD and LGG, so neuronal pathways in LGG tumours are expected. (iii) Drug-layer effects are deviations from the shared and class layers. (iv) Pathways are adjusted for proliferation, so cell-cycle effects are residual effects.

## Tier 1: q < 0.1 (21 effects, 15 distinct signals)

| Pathway | Layer / drug | Dir. | z | q | Class | Evidence summary | Refs |
|---|---|---|---|---|---|---|---|
| Presynaptic nicotinic AChR | drug:FLUOROURACIL | response | 4.52 | 0.004 | ARTEFACT? / CONTRA | Neuronal CHRN genes. Experimental work goes the other way: nicotine/α7-nAChR signalling protects GI cancer cells from 5-FU. Most likely innervation or tissue composition. | Chen 2015 Tumour Biol PMID 26136123; Dinicola 2013 Toxicol In Vitro PMID 24095863 |
| Signaling by FGFR3 mutants | drug:CISPLATIN | response | 3.78 | 0.035 | PLAUSIBLE (mixed) | Reflects FGFR3-high, luminal-papillary urothelial biology. Subtype data on cisplatin NAC conflict: basal tumours benefit most (Seiler), while FGFR/KRT subtyping is predictive in Ecke. | Seiler 2017 Eur Urol PMID 28390739; Ecke 2022 Int J Mol Sci PMID 35887247 |
| Cell cycle (residual after proliferation adjustment) | drug:CISPLATIN | response | 3.83 | 0.035 | KNOWN | High proliferation predicts chemosensitivity and pCR. The effect survives the proliferation covariate, which points to cycle/checkpoint content beyond Ki67-like signal. | Alba 2016 Oncologist PMID 26786263 |
| RIG-I/MDA5 induction of IFN-α/β | drug:CISPLATIN | resistance | −3.82 | 0.035 | KNOWN | Chronic IFN/IRDS signalling predicts resistance to DNA-damaging therapy. RIG-I-like receptors (DDX58, IFIH1) drive the therapy-induced IFN response. Acute type-I IFN can instead aid efficacy (Sistigu), so context matters. | Weichselbaum 2008 PNAS PMID 19001271; Ranoa 2016 Oncotarget PMID 27034163 |
| Antigen activates BCR → second messengers | drug:PACLITAXEL | response | 4.03 | 0.038 | KNOWN | B-cell/immune signatures and TILs predict pCR to taxane/anthracycline neoadjuvant chemotherapy in breast cancer. | Denkert 2010 JCO PMID 19917869; Iglesia 2014 Clin Cancer Res PMID 24916698 |
| Prostanoid ligand receptors | drug:PACLITAXEL | resistance | −3.65 | 0.090 | PLAUSIBLE | Signalling through COX-2/PGE2–EP receptors is linked to taxane resistance, and COX-2 inhibition restores paclitaxel sensitivity in vitro. | Hasegawa 2013 Oncol Rep PMID 24100466 |
| Fanconi anemia pathway | class:topoisomerase_II | response | 3.84 | 0.071 | PLAUSIBLE | FA/BRCA status governs sensitivity to DNA-damaging agents. Higher FA expression likely marks replication stress or HRD-associated transcription rather than repair capacity, so the direction is unexpected for a canonical repair-resistance model. Supported by Tier 2 regulation-of-FA (z = 3.32). | van der Heijden 2005 Clin Cancer Res PMID 16243825; Murai 2017 Int J Clin Oncol PMID 28643177 |
| Asparagine N-linked glycosylation | class:topoisomerase_II | response | 3.70 | 0.071 | NOVEL | Altered N-glycosylation appears in doxorubicin-resistant lines, and P-gp is a glycoprotein, but those findings run toward resistance. No evidence found for a response association. | Ji 2017 Oncotarget PMID 28077793 (indirect) |
| GABA, GABA-B receptor activation | class:antiangiogenic + drug:BEVACIZUMAB | response | 3.44, 3.55 | 0.089 | ARTEFACT? | Neuronal pathways in a cohort that includes LGG. Most likely brain-parenchyma or neural-subtype composition. | Verhaak 2010 Cancer Cell PMID 20129251 (neural subtype) |
| Inhibition of insulin secretion by adrenaline; Gβγ inhibition of Ca²⁺ channels | class:antiangiogenic + drug:BEVACIZUMAB | response | 3.56, 3.39 | 0.089 | ARTEFACT? | Neuronal/GPCR gene content (GNB/GNG, CACNA). Same composition caveat as above. | — |
| Netrin-1 signalling | class:antiangiogenic + drug:BEVACIZUMAB | response | 3.36 | 0.089 | PLAUSIBLE | Netrin-1/UNC5B is an axon-guidance cue that regulates sprouting angiogenesis and is proposed as a non-VEGF antiangiogenic axis. Also neuronal, so artefact risk remains. | Larrivée 2007 Genes Dev PMID 17908930; Pircher 2014 Oncology PMID 24401553 |
| Purine ribonucleoside monophosphate biosynthesis | class:antiangiogenic + drug:BEVACIZUMAB | resistance | −3.45 | 0.089 | NOVEL | No direct evidence linking de novo purine synthesis to bevacizumab outcome. May reflect metabolic adaptation to hypoxia. | — |
| Semaphorin interactions | class:antimetabolite | resistance | −3.81 | 0.095 | PLAUSIBLE | Semaphorin signalling promotes stemness and dissemination in pancreatic cancer (SEMA3C), and SEMA6A is altered in drug-resistant ovarian lines. No direct gemcitabine data. | Tomizawa 2023 Cancer Cell Int PMID 37537633; Prislei 2008 Mol Cancer Ther PMID 18187809 |

## Tier 2: relaxed class layers and drug totals (selected from top-5 lists; q ≥ 0.11)

Drug totals are dominated by the shared layer. The same pathways recur across platinum, taxane, anthracycline and alkylator totals, so they are summarised once in Tier 3. The items below are class- or drug-specific.

| Pathway | Layer / drug | Dir. | z | q | Class | Evidence summary | Refs |
|---|---|---|---|---|---|---|---|
| Regulation of FA pathway | class:topoII; DOXORUBICIN, ETOPOSIDE totals | response | 3.32; 3.23, 3.33 | 0.20; 0.32 | PLAUSIBLE | Concordant with the Tier 1 FA signal. | as FA above |
| Fanconi anemia pathway | class:other_dna (bleomycin) | resistance | −2.38 | 0.25 | PLAUSIBLE | Repair capacity opposing DNA-strand-break agents. Direction is opposite to topoII. | van der Heijden 2005 |
| PI3K cascade | class:microtubule | resistance | −3.12 | 0.20 | KNOWN | PIK3CA mutations reduce pCR to taxane-containing neoadjuvant chemotherapy. | Loibl 2016 Ann Oncol PMID 27177864; Guo 2020 Cancer Res Treat PMID 32019278 |
| Netrin-1 signalling | class:microtubule | resistance | −3.77 | 0.11 | NOVEL | Opposite sign to bevacizumab. No taxane data found. | — |
| Renal aquaporin water balance | class:microtubule; all drug totals | resistance | −3.59; −3.1 to −3.9 | 0.11–0.32 | ARTEFACT? | Kidney/vasopressin genes (AQP2, AVPR2). AQP1 polymorphisms have been tied to cisplatin outcome only in mesothelioma. | Senk 2019 Radiol Oncol PMID 30840592 |
| Glycolysis / glucose metabolism | class:microtubule; PACLITAXEL, DOCETAXEL totals | response | 2.80; 3.59, 3.31 | 0.29; 0.32 | CONTRA | Most preclinical data link glycolysis (LDHA, HIF-1α) to paclitaxel resistance. | Liu 2023 Cancer Gene Ther PMID 36241702 |
| Antigen cross-presentation | class:alkylating | response | 3.11 | 0.53 | KNOWN (mechanistic) | Cyclophosphamide efficacy depends on T-cell/DC priming and immunogenic cell death. | Viaud 2013 Science PMID 24264990 |
| PD-1 signalling | class:fluoropyrimidine | response | 2.81 | 0.78 | PLAUSIBLE | 5-FU depletes MDSCs and enhances T-cell–dependent immunity. A T-cell-rich context should therefore favour response. | Vincent 2010 Cancer Res PMID 20388795 |
| SMAD2/3/4 transcription; RIP-mediated NF-κB via DAI (ZBP1) | class:platinum | resistance | −2.81, −2.81 | 0.59 | PLAUSIBLE | TGF-β/SMAD activity promotes cisplatin resistance (ovarian). DAI-NF-κB is a cytosolic-DNA sensing arm that parallels the RIG-I/IRDS finding. | Mo 2020 J Interferon Cytokine Res PMID 32701410; Weichselbaum 2008 |
| NO stimulates guanylate cyclase | class:alkylating | resistance | −3.38 | 0.48 | NOVEL | No direct evidence for alkylator response. Possibly vascular/stromal composition. | — |
| CTLA4 inhibitory signalling | CISPLATIN, FLUOROURACIL, OXALIPLATIN totals | response | 2.83, 3.02, 2.69 | 0.32 | PLAUSIBLE | Marks activated or regulatory T-cell infiltrate, consistent with the TIL–chemo-response literature. | Denkert 2018 Lancet Oncol PMID 29233559 |

**Anchor checks (drug totals, none significant).** Most anchors point the expected way but are weak:
- NER is resistance-leaning for every drug (cisplatin z −1.30, carboplatin −1.16), consistent with ERCC1 (Olaussen 2006 NEJM PMID 16957145).
- ABC-family transporters are weakly negative for taxanes (−0.20 to −0.27), in line with Gottesman 2002 Nat Rev Cancer PMID 11902585 but uninformative.
- Pyrimidine metabolism with 5-FU is −0.74, a weak TYMS-direction signal (Johnston 1995 JNCI PMID 7563193; Peters 1995 Eur J Cancer PMID 7577040).
- IFN-α/β signalling with 5-FU is −1.31, IRDS-concordant.
- MGMT has no dedicated Reactome set. The "DNA damage reversal" pathway was not in the top lists for temozolomide (Hegi 2005 NEJM PMID 15758010 not testable at pathway level).

## Tier 3: shared layer (top 10 each direction; all q = 0.26–0.28)

| Pathway | Dir. | z | Class | Evidence summary | Refs |
|---|---|---|---|---|---|
| Metabolism of RNA | response | 3.23 | PLAUSIBLE | Biosynthetic and proliferative state beyond the proliferation covariate. | Alba 2016 |
| Glucose metabolism | response | 2.88 | CONTRA / NOVEL | See the taxane row in Tier 2. | — |
| Regulation of Fanconi anemia pathway | response | 2.80 | PLAUSIBLE | Recurs across topoII, anthracycline, alkylator and taxane totals. Most coherent DNA-repair theme. | van der Heijden 2005 |
| CTLA4 inhibitory signalling | response | 2.68 | PLAUSIBLE | T-cell infiltrate marker. | Denkert 2018 |
| BER single-nucleotide replacement | response | 2.64 | NOVEL | Base-excision repair activity is usually linked to resistance. | — |
| Signalling by FGFR1 mutants | response | 2.80 | NOVEL | — | — |
| Renal aquaporins | resistance | −3.31 | ARTEFACT? | Tissue composition. | — |
| Gap junction degradation | resistance | −3.14 | PLAUSIBLE | Loss of gap-junction communication reduces bystander cisplatin toxicity; degradation activity implies less coupling. | Arora 2018 Cancers PMID 30279363 |
| Downstream BCR signalling | resistance | −2.87 | CONTRA | Opposite to the paclitaxel-specific BCR response effect and to B-cell/pCR data. Possibly reflects lymphoid-aggregate composition. | Iglesia 2014 |
| MHC class II antigen presentation | resistance | −2.72 | CONTRA | Immune/APC content usually predicts chemo-response. The sign is unexpected and may reflect macrophage-rich (M2) stroma. | Denkert 2010 |
| Kainate receptor activation | resistance | −2.84 | ARTEFACT? | Neuronal genes. | — |
| Prostanoid ligand receptors | resistance | −2.77 | PLAUSIBLE | Also seen with paclitaxel (Tier 1) and carboplatin. | Hasegawa 2013 |
| POL switching | resistance | −2.77 | NOVEL | Replication polymerase switch (POLD/RFC/PCNA). Unexpected after proliferation adjustment. | — |

## Synthesis

**Recurring themes.**

1. **Interferon and innate nucleic-acid sensing → resistance.** This theme includes RIG-I/MDA5 with cisplatin (Tier 1), DAI-NF-κB with platinum, and IFN-α/β with 5-FU and doxorubicin (weak). It is the best literature-anchored theme, matching the IRDS model (Weichselbaum 2008; Minn 2015 Trends Immunol PMID 26604042).
2. **Adaptive immune context → response.** BCR with paclitaxel (Tier 1), CTLA4 and PD-1 signalling with platinum and 5-FU, and cross-presentation with alkylators all point to response, consistent with TIL/B-cell data (Denkert; Iglesia). This theme is internally inconsistent: shared MHC-II and shared BCR point toward resistance, which needs cell-type deconvolution.
3. **DNA repair / FA → response to topoisomerase II poisons.** This appears in Tier 1 and is reinforced by the regulation-of-FA pathway across drug totals. The direction most likely marks HRD/replication-stress biology rather than repair proficiency. NER→platinum resistance is directionally correct but weak.
4. **Proliferation residual → response.** Cell cycle with cisplatin and metabolism of RNA in the shared layer are consistent with classical chemosensitivity.
5. **Metabolism.** Glucose metabolism with taxanes and purine synthesis with bevacizumab are novel or contrary to the literature.

Neuronal, synaptic and renal pathways (nicotinic AChR, GABA, kainate, aquaporins) recur at the top of several lists. They should be treated as composition artefacts until reproduced with tumour-purity and tissue-signature adjustment.

**Five most credible candidate biomarkers for follow-up.**

1. **RIG-I/MDA5 → IFN-α/β (IRDS-like) activity → cisplatin resistance.** Tier 1, sign consistency 1.0, strong prior.
2. **B-cell receptor / B-cell signature → paclitaxel response.** Tier 1, strong breast neoadjuvant prior.
3. **Fanconi anemia pathway (and its regulation) → topoisomerase II poison response.** Tier 1 plus Tier 2 replication. Direction needs mechanistic clarification against HRD scores.
4. **Residual cell-cycle activity → cisplatin response.** Tier 1, classical prior. Test whether it adds to Ki67/proliferation indices.
5. **Prostanoid (PGE2/EP) receptor signalling → taxane/platinum resistance.** Tier 1 for paclitaxel, also in the shared layer and carboplatin totals. Plausible and druggable (COX-2 inhibition).

The FGFR3-cisplatin finding is a sixth, bladder-specific candidate, but it needs a subtype-aware re-analysis.
