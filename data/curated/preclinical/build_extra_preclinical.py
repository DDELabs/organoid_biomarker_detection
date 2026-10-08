#!/usr/bin/env python3
"""Build additional curated pre-clinical sets (patient-derived xenografts / organoids / primary cultures).

Re-uses the helpers of build_preclinical.py (download cache, HGNC mapping, drug-name harmonisation,
byte-reproducible writer) and writes, per set, data/curated/preclinical/<SET>/
  expression.tsv.gz      genes (HGNC symbol, protein-coding) x models, log2 scale
  response.tsv           sample, drug, response, metric, direction, n_screens, screens (+ extra columns)
  patient_response.tsv   (only when the matched patient's clinical response is public)
  mrecist.tsv            (PDXE only)
  SOURCE.md
then regenerates CATALOG_PRECLINICAL.tsv for all sets.

    python3 data/curated/preclinical/build_extra_preclinical.py                 # all extra sets
    python3 data/curated/preclinical/build_extra_preclinical.py pdx_novartis_gao2015 --raw /tmp/cache
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_preclinical as bp  # noqa: E402

SPRINGER = "https://static-content.springer.com/esm/{}"
bp.SOURCES.update({
    # Gao et al. 2015 Nat Med 21:1318 (Novartis PDX Encyclopedia) - Supplementary Table S1 (nm.3954-S2.xlsx)
    "pdxe_xlsx": (SPRINGER.format("art%3A10.1038%2Fnm.3954/MediaObjects/41591_2015_BFnm3954_MOESM10_ESM.xlsx"),
                  "Gao2015_nm3954_TableS1.xlsx"),
})

LOWER = "lower = more sensitive"
HIGHER = "higher = more sensitive"


def write(name, expr, resp, md, catalog, direction=LOWER, extra_files=None):
    """bp.write_set with a per-set response direction; ensures n_screens / screens columns."""
    resp = resp.copy()
    if "n_screens" not in resp:
        resp["n_screens"] = 1
    if "screens" not in resp:
        resp["screens"] = resp["sample"]
    bp.DIRECTION = direction
    try:
        return bp.write_set(name, expr, resp, md, catalog, extra_files=extra_files)
    finally:
        bp.DIRECTION = LOWER


# ============================================================================ Novartis PDXE (Gao 2015)
# Novartis compound codes -> generic names (INN where one exists; otherwise the code is kept).
PDXE_DRUGS = {
    "5FU": "FLUOROURACIL", "BGJ398": "INFIGRATINIB", "BKM120": "BUPARLISIB", "BYL719": "ALPELISIB",
    "CGM097": "CGM097", "CKX620": "CKX620", "CLR457": "CLR457", "HDM201": "SIREMADLIN", "HSP990": "HSP990",
    "INC280": "CAPMATINIB", "INC424": "RUXOLITINIB", "LCL161": "LCL161", "LDE225": "SONIDEGIB",
    "LDK378": "CERITINIB", "LEE011": "RIBOCICLIB", "LFA102": "LFA102", "LFW527": "LFW527", "LGH447": "PIM447",
    "LGW813": "LGW813", "LJC049": "LJC049", "LJM716": "ELGEMTUMAB", "LKA136": "LKA136", "LLM871": "LLM871",
    "TAS266": "TAS266", "WNT974": "WNT-974", "abraxane": "NAB-PACLITAXEL", "binimetinib": "BINIMETINIB",
    "binimetinib-3.5mpk": "BINIMETINIB", "cetuximab": "CETUXIMAB", "dacarbazine": "DACARBAZINE",
    "encorafenib": "ENCORAFENIB", "erlotinib": "ERLOTINIB", "everolimus": "EVEROLIMUS",
    'figitumumab"': "FIGITUMUMAB", "gemcitabine-50mpk": "GEMCITABINE", "gemcitabine": "GEMCITABINE", "paclitaxel": "PACLITAXEL",
    "tamoxifen": "TAMOXIFEN", "trametinib": "TRAMETINIB", "trastuzumab": "TRASTUZUMAB",
}
PDXE_TUMOUR = {"BRCA": "breast", "CRC": "colorectal", "NSCLC": "lung (NSCLC)", "GC": "gastric",
               "PDAC": "pancreas", "CM": "melanoma"}


def mrecist(best, best_avg):
    """mRECIST call of Gao et al. 2015 (Online Methods), on % tumour-volume change from day 0."""
    if best < -95 and best_avg < -40:
        return "CR"
    if best < -50 and best_avg < -20:
        return "PR"
    if best < 35 and best_avg < 30:
        return "SD"
    return "PD"


def build_pdxe():
    path = bp.fetch("pdxe_xlsx")
    fpkm = pd.read_excel(path, sheet_name="RNAseq_fpkm", index_col=0)
    met = pd.read_excel(path, sheet_name="PCT curve metrics")
    raw = pd.read_excel(path, sheet_name="PCT raw data", usecols=["Model", "Tumor Type"]).drop_duplicates()
    ttype = raw.drop_duplicates("Model").set_index("Model")["Tumor Type"]

    met = met[met["Treatment"] != "untreated"].copy()
    met["drug"] = [" + ".join(PDXE_DRUGS[p.strip()] for p in t.split(" + ")) for t in met["Treatment"]]
    met["mrecist"] = [mrecist(b, a) for b, a in zip(met["BestResponse"], met["BestAvgResponse"])]
    met["responder"] = met["mrecist"].isin(["CR", "PR"]).astype(int)
    met["tumour_type"] = met["Model"].map(ttype)
    met["paper_category"] = met["ResponseCategory"]
    agree = (met["mrecist"] == met["paper_category"].str.split("-").str[0]).mean()

    mr = met.rename(columns={"Model": "sample", "BestAvgResponse": "best_avg_response",
                             "BestResponse": "best_response", "Treatment": "treatment"})
    mr = mr[["sample", "drug", "best_avg_response", "mrecist", "responder", "best_response", "treatment",
             "tumour_type", "paper_category", "TimeToDouble", "Day_Last"]].rename(
        columns={"TimeToDouble": "time_to_double", "Day_Last": "day_last"})
    mr = mr.sort_values(["drug", "sample"])

    # response.tsv: BestAvgResponse. The low-dose arms (binimetinib-3.5mpk) are kept only in mrecist.tsv;
    # gemcitabine-50mpk is the only gemcitabine arm and is kept as GEMCITABINE.
    r = met[met["Treatment"] != "binimetinib-3.5mpk"]
    resp = pd.DataFrame({"sample": r["Model"], "drug": r["drug"], "response": r["BestAvgResponse"],
                         "metric": "BestAvgResponse", "n_screens": 1, "screens": r["Model"] + "|" + r["Treatment"],
                         "mrecist": r["mrecist"], "treatment_type": r["Treatment type"],
                         "tumour_type": r["tumour_type"]})
    assert not resp.duplicated(["sample", "drug"]).any()

    expr = bp.map_symbols(np.log2(fpkm.astype(float) + 1))
    samples = pd.DataFrame({"sample": sorted(set(expr.columns) | set(ttype.index))})
    samples["tumour_type"] = samples["sample"].map(ttype)
    samples["cancer"] = samples["tumour_type"].map(PDXE_TUMOUR)
    samples["has_expression"] = samples["sample"].isin(expr.columns).astype(int)
    samples["has_pct"] = samples["sample"].isin(resp["sample"]).astype(int)

    both = set(expr.columns) & set(resp["sample"])
    per_type = resp[resp["sample"].isin(both)].drop_duplicates("sample")["tumour_type"].value_counts()
    cnt = mr["mrecist"].value_counts()
    md = f"""# Novartis PDX Encyclopedia (Gao et al. 2015)

**Paper**: Gao H, Korn JM, Ferretti S, et al. *High-throughput screening using patient-derived tumor xenografts to
predict clinical trial drug response.* Nat Med 2015;21(11):1318-1325. doi:10.1038/nm.3954 (PMID 26479923).

**Accession**: Supplementary Table S1 (`nm.3954-S2.xlsx`, Springer static content), sheets `RNAseq_fpkm`,
`PCT curve metrics`, `PCT raw data`. Raw sequencing is not used.

{bp.provenance(["pdxe_xlsx", "gene_info"])}

## Design
PDX clinical trial (PCT), 1 x 1 x 1 design: one mouse per PDX model x treatment, tumour volume followed for
~3 weeks or more. {{n_resp_models}} models are treated; tumour types of the models with both data types:
{", ".join(f"{k} {v}" for k, v in per_type.items())} (BRCA breast, CRC colorectal, GC gastric, NSCLC lung,
PDAC pancreas, CM cutaneous melanoma).

## Derivation
* Response (`response.tsv`): **BestAvgResponse** from the paper (minimum over t >= 10 d of the running mean of the
  % tumour-volume change from day 0). **Lower = more sensitive** (negative = shrinkage). One row per model x
  treatment arm; the `untreated` arm and the low-dose `binimetinib-3.5mpk` arm are excluded (the latter is kept in
  `mrecist.tsv`). Combination arms are kept with the components joined by ` + ` (e.g. `ALPELISIB + ELGEMTUMAB`),
  `treatment_type` = single / combo. Novartis codes were mapped to INNs where one exists
  (BYL719 alpelisib, BKM120 buparlisib, LEE011 ribociclib, INC280 capmatinib, INC424 ruxolitinib, LDK378 ceritinib,
  LDE225 sonidegib, BGJ398 infigratinib, HDM201 siremadlin, LGH447 PIM447, LJM716 elgemtumab, WNT974,
  abraxane nab-paclitaxel); research codes without an INN are kept (CGM097, CLR457, CKX620, HSP990, LCL161,
  LFA102, LFW527, LGW813, LJC049, LKA136, LLM871, TAS266).
* `mrecist.tsv`: mRECIST class recomputed from BestResponse / BestAvgResponse with the thresholds of the paper's
  Online Methods: CR if BestResponse < -95% and BestAvgResponse < -40%; PR if BestResponse < -50% and
  BestAvgResponse < -20%; SD if BestResponse < 35% and BestAvgResponse < 30%; otherwise PD. `responder` = 1 for
  CR/PR. Counts: {", ".join(f"{k} {cnt.get(k, 0)}" for k in ["CR", "PR", "SD", "PD"])}. The paper's own
  `ResponseCategory` (which also tracks progression after the best response, e.g. `SD-->PD`) is kept as
  `paper_category`; its first class agrees with the recomputed call in {agree:.1%} of arms.
* Expression: `RNAseq_fpkm` (FPKM, Novartis pipeline) -> log2(FPKM + 1); symbols re-mapped to current HGNC
  protein-coding symbols. Units: **log2(FPKM + 1)**. Mouse reads were removed by the authors.

## Caveats
* No matched patient clinical response is public for PDXE (there is no `patient_response.tsv`).
* Single animal per arm: noisy. Treat the binary responder call as the most robust endpoint.
* Response rates are low (CR+PR ~ 6% of arms), and most arms are targeted agents.

## Sample-ID matching
Model IDs (`X-1004` ...) shared by both sheets. Expression models: {{n_expr_models}}; treated models:
{{n_resp_models}}; **overlap n = {{n_both}}**. Genes: {{n_genes}}; drugs/regimens: {{n_drugs}}.
"""
    return write("pdx_novartis_gao2015", expr, resp, md,
                 {"type": "PDX", "tissue": "pan-cancer", "source": "Nat Med supplement"},
                 extra_files={"mrecist.tsv": mr, "samples.tsv": samples})


BUILDERS = {
    "pdx_novartis_gao2015": build_pdxe,
}

bp.META.update({
    "pdx_novartis_gao2015": dict(model="PDX", tissue="pan-cancer (breast, CRC, gastric, NSCLC, pancreas, melanoma)",
                                 expression_units="log2(FPKM+1)", reference="Gao 2015 Nat Med 21:1318",
                                 accession="nm.3954 Supplementary Table S1"),
})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sets", nargs="*", help=f"subset of {list(BUILDERS)}")
    ap.add_argument("--raw", default=str(bp.RAW), help="download cache directory")
    a = ap.parse_args()
    bp.RAW = Path(a.raw)
    for s in a.sets or BUILDERS:
        BUILDERS[s]()
    bp.write_catalog()


if __name__ == "__main__":
    sys.exit(main())
