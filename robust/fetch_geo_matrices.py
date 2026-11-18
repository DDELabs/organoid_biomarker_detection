"""Download GEO series-matrix files and supplementary-file listings for the training-table cohorts."""
import re, sys, time, urllib.request
from pathlib import Path
OUT = Path(__file__).resolve().parent.parent / "data" / "external" / "geo_train"
COHORTS = {
 "bladder": ["GSE279390","GSE87304","GSE235067","GSE48276","GSE70691","GSE169455","GSE31684"],
 "breast": ["GSE194040","GSE164458","GSE319641","GSE96058","GSE25066","GSE186102","GSE163882","GSE123845","GSE240671","GSE116335","GSE234519"],
 "colorectal": ["GSE156915","GSE39582","GSE14333","GSE103479","GSE104645","GSE87211","GSE72970","GSE28702"],
 # lung: priority 1-3 rows of nar/16_plan/COHORTS_TO_ADD.tsv, plus GSE3141, GSE81089, GSE126044 named in PLAN.md W1
 "lung": ["GSE14814","GSE68465","GSE37745","GSE29013","GSE31210","GSE50081","GSE8894","GSE3141","GSE30219","GSE72094",
          "GSE41271","GSE157011","GSE33072","GSE161537","GSE162520","GSE135222","GSE126044","GSE218989","GSE207422",
          "GSE232569","GSE81089"],
}
def get(url, dest=None, tries=4):
    for i in range(tries):
        try:
            if dest is None:
                return urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "replace")
            urllib.request.urlretrieve(url, dest); return dest
        except Exception as e:
            err = e; time.sleep(2 ** (i + 1))
    print("FAILED", url, err); return None
def base(gse): return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/"
if __name__ == "__main__":
    want = set(sys.argv[1:])   # optional: only these cancer groups, e.g. "lung"
    for cancer, ids in COHORTS.items():
        if want and cancer not in want:
            continue
        for gse in ids:
            d = OUT / gse; d.mkdir(parents=True, exist_ok=True)
            for sub in ("matrix", "suppl"):
                html = get(base(gse) + sub + "/") or ""
                rows = re.findall(r'href="([^"?/][^"]*)"[^\n]*?(\d[\d.]*[KMG]?)\s*$', html, flags=re.M)
                (d / f"{sub}_listing.txt").write_text("\n".join(f"{n}\t{s}" for n, s in rows))
                if sub == "matrix":
                    for n, s in rows:
                        if n.endswith("series_matrix.txt.gz") and not (d / n).exists():
                            get(base(gse) + "matrix/" + n, d / n)
            print(cancer, gse, [p.name for p in d.glob("*matrix.txt.gz")], flush=True)
