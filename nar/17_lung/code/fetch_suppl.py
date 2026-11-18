"""Download the small supplementary files (expression tables) the lung builders need (W1 step 1/3)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "robust"))
import fetch_geo_matrices as F
FILES = {"GSE161537": ["GSE161537_nivobio_log2cpm.csv.gz"], "GSE162520": ["GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz"],
         "GSE135222": ["GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"], "GSE126044": ["GSE126044_counts.txt.gz"],
         "GSE218989": ["GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz"],
         "GSE207422": ["GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx"],
         "GSE81089": ["GSE81089_FPKM_cufflinks.tsv.gz", "GSE81089_readcounts_featurecounts.tsv.gz"],
         "GSE232569": ["filelist.txt", "GSE232569_RAW.tar"], "GSE14814": ["filelist.txt"]}
for g, fs in FILES.items():
    for f in fs:
        dest = F.OUT / g / f
        if not dest.exists():
            print(F.get(F.base(g) + "suppl/" + f, dest), flush=True)
