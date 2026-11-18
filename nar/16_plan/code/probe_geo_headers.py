"""Stream only the header of each GEO series matrix (no full download) and summarise sample characteristics.
Usage: python probe_geo_headers.py OUTDIR GSE1 GSE2 ...   -> OUTDIR/<GSE>.txt  (title, platforms, n, characteristic keys + values)"""
import gzip, io, re, sys, urllib.request
from collections import defaultdict, Counter
out = sys.argv[1]
def base(g): return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:-3]}nnn/{g}/matrix/"
def probe(g):
    lines = []
    try:
        html = urllib.request.urlopen(base(g), timeout=60).read().decode()
        files = sorted(set(re.findall(r'href="([^"]*series_matrix\.txt\.gz)"', html)))
        rep = [f"# {g} matrix files: {files}"]
        for f in files:
            r = urllib.request.urlopen(base(g) + f, timeout=120)
            gz = gzip.GzipFile(fileobj=r)
            chars = defaultdict(Counter); n = 0; title = plat = ""
            for raw in io.TextIOWrapper(gz, encoding="utf-8", errors="replace"):
                if raw.startswith("!series_matrix_table_begin"): break
                p = raw.rstrip("\n").split("\t")
                k = p[0]
                v = [x.strip('"') for x in p[1:]]
                if k == "!Series_title": title = v[0]
                elif k == "!Series_platform_id": plat = ",".join(v)
                elif k == "!Sample_geo_accession": n = len(v)
                elif k in ("!Sample_characteristics_ch1", "!Sample_source_name_ch1", "!Sample_description"):
                    for x in v:
                        key, _, val = x.partition(":") if ":" in x and k == "!Sample_characteristics_ch1" else (k, "", x)
                        chars[key.strip()][val.strip()[:40]] += 1
            r.close()
            rep.append(f"## {f}\ntitle: {title}\nplatform: {plat}\nn_samples: {n}")
            for key, c in chars.items():
                rep.append(f"  [{key}] distinct={len(c)} :: " + " | ".join(f"{a}({b})" for a, b in c.most_common(6)))
        lines = rep
    except Exception as e:
        lines = [f"# {g} ERROR {e}"]
    return lines
for g in sys.argv[2:]:
    import time
    for t in range(4):
        lines = probe(g)
        txt = "\n".join(lines)
        if "n_samples: 0" not in txt and "matrix files: []" not in txt and "ERROR" not in txt: break
        time.sleep(3 * (t + 1))
    open(f"{out}/{g}.txt", "w").write("\n".join(lines) + "\n")
    print(g, "done", flush=True)
