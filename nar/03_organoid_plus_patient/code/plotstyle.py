"""Figure style: colour-blind-safe Okabe-Ito colours, fonts >= 8 pt, panel letters, PNG (300 dpi) + PDF + caption."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CB = ["#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2", "#F0E442", "#000000"]
GREY = "#9a9a9a"
MODEL_COL = {"organoid-only": "#E69F00", "patient-only": "#56B4E9", "joint": "#009E73", "joint (w=0.5)": "#CC79A7"}
SETNAME = {"colorectal_vandewetering2015": "CRC organoids (van de Wetering 2015)",
           "colorectal_pdx_isella2017": "CRC PDX (Isella 2017)",
           "bladder_lee2018": "Bladder organoids (Lee 2018)",
           "pancreas_tiriac2018": "Pancreas organoids (Tiriac 2018)",
           "pancreas_shi2022": "Pancreas organoids (Shi 2022)",
           "liver_ji2023": "Liver organoids (Ji 2023)",
           "liver_broutier2017": "Liver organoids (Broutier 2017)",
           "sarcoma_alshihabi2024": "Sarcoma organoids (Al Shihabi 2024)",
           "ovarian_vias2023": "Ovarian organoids (Vias 2023)",
           "ctrpv2_ccle": "Cell lines CTRPv2", "prism_repurposing": "Cell lines PRISM"}


def setup():
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 8,
                         "ytick.labelsize": 8, "legend.fontsize": 8, "pdf.fonttype": 42, "font.family": "DejaVu Sans",
                         "axes.spines.top": False, "axes.spines.right": False})


def letters(axes, x=-0.12, y=1.06):
    for i, a in enumerate(axes):
        a.text(x, y, "abcdefghijklmnop"[i], transform=a.transAxes, fontsize=11, fontweight="bold", va="bottom")


def save(fig, path, caption):
    path = str(path)
    fig.savefig(path + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(path + ".pdf", bbox_inches="tight")
    plt.close(fig)
    with open(path + ".caption.md", "w") as f:
        f.write(caption.strip() + "\n")
