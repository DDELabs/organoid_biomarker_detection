"""Extract the sheets of the paper's source data (MOESM3) that this analysis needs into TSV files."""
from pathlib import Path
import pandas as pd
B = Path(__file__).resolve().parents[1]
SRC = B.parents[1] / "_cache" / "41467_2020_19313_MOESM3_ESM.xlsx"
OUT = B / "_data" / "moesm3"; OUT.mkdir(parents=True, exist_ok=True)
spec = {"Figure 2a-f": 2, "Figure 2a,d": 2, "Figure 2e": 3, "Figure 2h": 2, "Figure 4c": None,
        "Figure 4f": 2, "FigureS5": 2, "Figure 2b": 3, "Figure 2g": 2, "Figure 4b": None, "Figure 3e": 3}
for sheet, hdr in spec.items():
    df = pd.read_excel(SRC, sheet_name=sheet, header=hdr)
    df.to_csv(OUT / (sheet.replace(" ", "_").replace(",", "").replace("-", "_") + ".tsv"), sep="\t", index=False)
    print(sheet, df.shape)
