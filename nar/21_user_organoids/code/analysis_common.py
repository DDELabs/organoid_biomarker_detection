"""Small helpers shared by the two analysis runners (no heavy imports here)."""
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
FOLDER = HERE.parent
REPO = FOLDER.parents[1]
NAR = REPO / "nar"
PRECLIN_REAL = REPO / "data" / "curated" / "preclinical"
CACHE = FOLDER / "_cache"                      # shared downloads (TCGA expression etc.), safe to delete

# what the user may type for --cancer  ->  TCGA projects used for the patient side
CANCER_TCGA = {
    "colorectal": ["COAD", "READ"], "colon": ["COAD"], "rectal": ["READ"], "crc": ["COAD", "READ"], "coad": ["COAD"], "read": ["READ"],
    "bladder": ["BLCA"], "blca": ["BLCA"], "pancreas": ["PAAD"], "pancreatic": ["PAAD"], "paad": ["PAAD"],
    "liver": ["LIHC"], "hepatocellular": ["LIHC"], "lihc": ["LIHC"], "ovarian": ["OV"], "ovary": ["OV"], "ov": ["OV"],
    "lung": ["LUAD", "LUSC"], "luad": ["LUAD"], "lusc": ["LUSC"], "breast": ["BRCA"], "brca": ["BRCA"],
    "sarcoma": ["SARC"], "sarc": ["SARC"], "stomach": ["STAD"], "gastric": ["STAD"], "stad": ["STAD"],
    "kidney": ["KIRC"], "kirc": ["KIRC"], "head and neck": ["HNSC"], "hnsc": ["HNSC"], "esophageal": ["ESCA"], "esca": ["ESCA"],
    "glioma": ["LGG"], "glioblastoma": ["GBM"], "gbm": ["GBM"], "melanoma": ["SKCM"], "skcm": ["SKCM"],
    "prostate": ["PRAD"], "prad": ["PRAD"], "endometrial": ["UCEC"], "uterine": ["UCEC"], "cervical": ["CESC"],
}
# the group name used by nar/03 units (UNITS[i]["cancer"])
CANCER_UNIT_GROUP = {"COAD": "colorectal", "READ": "colorectal", "BLCA": "bladder", "PAAD": "pancreas", "LIHC": "liver", "OV": "ovarian",
                     "LUAD": "lung", "LUSC": "lung", "SARC": "sarcoma"}


def tcga_projects(cancer):
    c = str(cancer or "").strip().lower()
    if c in CANCER_TCGA:
        return CANCER_TCGA[c]
    for k, v in CANCER_TCGA.items():             # e.g. 'colorectal cancer', 'bladder urothelial carcinoma'
        if len(k) > 3 and k in c:
            return v
    return []


def unit_group(projects):
    return next((CANCER_UNIT_GROUP[p] for p in projects if p in CANCER_UNIT_GROUP), None)


def link(src, dst):
    if dst.is_symlink() or dst.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(src, dst)


def make_raw_cache():
    """Folder with symlinks to every raw file the repo already downloaded (TCGA expression, BCR drug tables, STRING,
    survival table). New downloads (a cancer not downloaded before) land here, never in nar/02 or nar/09."""
    raw = CACHE / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    for src in (NAR / "02_organoid_only" / "_data" / "raw", NAR / "09_registered_test" / "_data" / "tcga_response" / "raw"):
        if src.exists():
            for p in src.iterdir():
                link(p, raw / p.name)
    return raw


def make_registry(set_dir, run_dir, name):
    """A preclinical folder that holds every curated set plus the user's set (symlinks), so the repo loaders find all of them."""
    reg = Path(run_dir) / "_registry"
    reg.mkdir(parents=True, exist_ok=True)
    for p in PRECLIN_REAL.iterdir():
        if p.is_dir():
            link(p, reg / p.name)
    if (PRECLIN_REAL / name).exists():
        raise SystemExit(f"The set name '{name}' is already used by a curated set in data/curated/preclinical. Pick another name.")
    link(Path(set_dir).resolve(), reg / name)
    return reg
