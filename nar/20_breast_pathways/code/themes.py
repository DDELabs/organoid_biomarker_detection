"""Plain-language theme for each pathway, by keywords in the Reactome name. First matching theme wins. Crude on purpose; full list saved for review."""
from common import *
IMMUNE = ["IMMUN", "INTERLEUKIN", "INTERFERON", "CYTOKINE", "T_CELL", "B_CELL", "TCR", "BCR", "ANTIGEN", "MHC", "COMPLEMENT", "TOLL_LIKE", "INFLAMMASOME",
          "FCGAMMA", "FCERI", "CD28", "PD_1", "CHEMOKINE", "NEUTROPHIL", "ADAPTIVE", "INNATE", "LEUKOCYTE", "PHAGOCYT", "DEFENSIN", "TLR", "NLR", "RIG_I", "CTLA4",
          "COSTIMULATION", "ZAP70", "LCK", "PD1", "TNF", "NF_KB", "TRAF6", "MYD88", "LYMPH", "MAST_CELL", "ANTIMICROBIAL", "PEPTIDE_LIGAND_BINDING_RECEPTORS", "CD22", "GRB2_EVENTS_IN_ERBB"]
IMMUNE = [k for k in IMMUNE if k not in ("PEPTIDE_LIGAND_BINDING_RECEPTORS", "GRB2_EVENTS_IN_ERBB")]
SIGNAL = ["SIGNAL", "RECEPTOR", "GPCR", "G_PROTEIN", "KINASE", "MAPK", "PI3K", "AKT", "MTOR", "WNT", "NOTCH", "HEDGEHOG", "TGF", "EGFR", "ERBB", "FGFR", "VEGF", "RAS_", "RHO_",
          "RAC", "CDC42", "GTPASE", "ESTROGEN", "NUCLEAR_RECEPTOR", "PDGF", "INSULIN", "IGF", "NTRK", "MET_", "SMAD", "BMP", "ERK", "PKA", "CAMP", "CALCIUM", "PLC", "GAB1", "SHC", "SOS", "GRB2", "FRS", "TRK", "EPH", "SEMAPHORIN", "INTEGRIN", "ADHESION", "ECM"]
METAB = ["METABOLISM", "BIOSYNTHESIS", "SYNTHESIS", "OXIDATION", "LIPID", "CHOLESTEROL", "GLYCOLYSIS", "GLUCONEOGENESIS", "TCA", "RESPIRATORY", "ELECTRON_TRANSPORT", "AMINO_ACID",
         "FATTY_ACID", "STEROID", "BILE", "VITAMIN", "HEME", "NUCLEOTIDE", "PURINE", "PYRIMIDINE", "GLYCOSAMINOGLYCAN", "SPHINGOLIPID", "PHOSPHOLIPID", "GLUTATHIONE", "COENZYME",
         "GLYCOGEN", "UREA_CYCLE", "KETONE", "FOLATE", "ONE_CARBON", "PENTOSE", "CARNITINE", "SULFUR", "ASPARTATE", "GLUTAMATE", "SERINE", "TRYPTOPHAN", "XENOBIOTIC", "CYTOCHROME_P450", "ACETYLCHOLINE_SYNTHESIS"]
def theme(name):
    n = name.replace("REACTOME_", "")
    if any(k in n for k in PROLIF_KW): return "proliferation"
    if any(k in n for k in IMMUNE): return "immune"
    if any(k in n for k in METAB): return "metabolism"
    if any(k in n for k in SIGNAL): return "signalling"
    return "other"

BROAD = ["CDC6", "CDC20", "CDC25", "CDC14", "CDK", "ORC", "MEIOSIS", "G0_AND", "CYCLIN", "KINETOCHORE", "CENP", "SPINDLE", "CHROMOSOME", "TELOMERE", "S_PHASE", "MITOSIS", "APC_C", "DNA_REPAIR",
         "MITOTIC", "CONDENSATION", "SISTER_CHROMATID", "MINICHROMOSOME", "MCM", "FANCONI", "HOMOLOGY_DIRECTED"]
def broad_prolif(name):
    n = name.replace("REACTOME_", "")
    return theme(name) == "proliferation" or any(k in n for k in BROAD)
