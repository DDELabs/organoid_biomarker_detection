"""Record size + sha256 of the 02 GI raw inputs (02 REPORT_GI lists sizes only). Read-only on 02."""
import hashlib, os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "02_organoid_only" / "_data" / "raw"
out = ROOT / "12_audit" / "tables" / "data_versions_02GI_sha256.tsv"
rows = []
for p in sorted(RAW.rglob("*")):
    if p.is_file():
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for c in iter(lambda: f.read(1 << 20), b""):
                h.update(c)
        rows.append(f"{p.relative_to(ROOT.parent)}\t{p.stat().st_size}\t{h.hexdigest()}")
out.write_text("file\tbytes\tsha256\n" + "\n".join(rows) + "\n")
print(len(rows), "files")
