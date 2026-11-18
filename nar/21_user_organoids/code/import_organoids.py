"""One command for the user's own organoids.

  python3 nar/21_user_organoids/code/import_organoids.py <input_folder> <set_name> [--cancer colorectal] [options]

Step 1  read + check + score the files  -> nar/21_user_organoids/imports/<set_name>/   (IMPORT_REPORT.md says what was found)
Step 2  paper method (organoid-only)    -> nar/21_user_organoids/runs/<set_name>/organoid_only/     (reuses nar/02_organoid_only/code)
Step 3  organoid + patient method       -> nar/21_user_organoids/runs/<set_name>/organoid_plus_patient/ (reuses nar/03_organoid_plus_patient/code)
Step 4  plain-language summary          -> nar/21_user_organoids/runs/<set_name>/ANALYSIS_SUMMARY.md
Steps 2-3 run one after the other, each as ONE process on ONE CPU.  Use --import-only to stop after step 1.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import analysis_common as AC  # noqa: E402
import load_user_organoids as L  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("set_name", help="short name without spaces, e.g. mylab_colorectal2025")
    ap.add_argument("--cancer", help="colorectal, bladder, pancreas, liver, ovarian, lung, breast, sarcoma ... or a TCGA code (COAD). "
                                     "If missing, taken from the cancer_type column of samples.tsv")
    ap.add_argument("--metric", help="what the drug numbers are (AUC, IC50, viability ...) when the file does not say")
    ap.add_argument("--direction", choices=["lower", "higher"], help="'lower' = a lower number means more sensitive (AUC, IC50). Default: guessed from the metric name")
    ap.add_argument("--value-type", choices=["counts", "tpm", "log"], help="force the expression value type instead of auto-detection")
    ap.add_argument("--drugs", default="", help="comma list: analyse only these drugs, e.g. FLUOROURACIL,OXALIPLATIN")
    ap.add_argument("--features", default="shortlist", choices=["shortlist", "all", "both"], help="nar/03 feature set (default: network shortlist)")
    ap.add_argument("--fast", action="store_true", help="fewer bootstraps/permutations in the paper method (1000/100 instead of 10000/1000)")
    ap.add_argument("--import-only", action="store_true", help="only import and score; do not run the two methods")
    ap.add_argument("--skip-organoid-only", action="store_true")
    ap.add_argument("--skip-joint", action="store_true")
    ap.add_argument("--out-root", default=str(AC.FOLDER), help="where imports/ and runs/ are created (default: nar/21_user_organoids)")
    a = ap.parse_args()

    root = Path(a.out_root)
    imp = root / "imports" / a.set_name
    run = root / "runs" / a.set_name
    print(f"[1/4] Importing {a.input_folder} ...", flush=True)
    res = L.import_folder(a.input_folder, a.set_name, out_dir=imp, cancer=a.cancer, metric=a.metric, direction=a.direction, value_type=a.value_type)
    print(res["report"])
    print(f"Import report: {imp / 'IMPORT_REPORT.md'}", flush=True)
    if res["status"] != "ok":
        sys.exit(1)
    cancer = a.cancer or res.get("cancer")
    if a.import_only:
        return
    if not cancer or not AC.tcga_projects(cancer):
        print("The patient analysis needs a cancer type that matches a TCGA project. Re-run with --cancer (for example --cancer colorectal). "
              "The import itself is finished.")
        return
    run.mkdir(parents=True, exist_ok=True)
    py = sys.executable
    base = [str(imp), a.set_name, cancer]
    codes = {}
    if not a.skip_organoid_only:
        print("[2/4] Paper method (organoid-only) ...", flush=True)
        cmd = [py, str(HERE / "run_organoid_only_user.py"), *base, str(run / "organoid_only"), "--drugs", a.drugs] + (["--fast"] if a.fast else [])
        codes["organoid_only"] = subprocess.call(cmd)
    if not a.skip_joint:
        print("[3/4] Organoid + patient method ...", flush=True)
        cmd = [py, str(HERE / "run_joint_user.py"), *base, str(run / "organoid_plus_patient"), "--drugs", a.drugs, "--features", a.features]
        codes["joint"] = subprocess.call(cmd)
    print("[4/4] Writing summary ...", flush=True)
    write_summary(a.set_name, cancer, imp, run, codes)
    print(f"Summary: {run / 'ANALYSIS_SUMMARY.md'}")


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return "n/a"
    if isinstance(x, (int,)) or (isinstance(x, float) and float(x).is_integer() and abs(x) > 50):
        return f"{int(x)}"
    return f"{x:.{nd}g}"


def write_summary(set_name, cancer, imp, run, codes):
    lines = [f"# Analysis summary: {set_name} ({cancer})", "",
             f"Import report: `{imp / 'IMPORT_REPORT.md'}`. Numbers below are copied from the tables in this folder; nothing is typed by hand.", ""]
    s = run / "organoid_only" / "organoid_only_summary.tsv"
    lines += ["## 1. Paper method: organoids only, then TCGA patients", "",
              "How it works: pick Reactome pathways that sit close to the drug's targets in the protein network (z <= -1.28); score them in your "
              "organoids; train Ridge regression on your organoid drug response; apply it to TCGA patients who got the drug; split them at the "
              "median predicted response; compare survival (log-rank test, hazard ratio = HR). HR below 1 means the group predicted to be "
              "sensitive lived longer.", ""]
    if s.exists() and len(pd.read_csv(s, sep="\t")):
        t = pd.read_csv(s, sep="\t")
        lines += ["| Drug | TCGA cohort | Organoids | Treated patients (events) | Top pathway (Ridge) | HR (95% CI) | log-rank p | Untreated p |", "|---|---|---|---|---|---|---|---|"]
        for _, r in t.iterrows():
            lines.append(f"| {r['drug']} | {r['cohort']} | {_fmt(r['n_organoids'])} | {_fmt(r['n_treated_patients'])} ({_fmt(r['treated_events'])}) | "
                         f"{r['ridge_top_pathway']} | {_fmt(r['HR'])} ({_fmt(r['HR_low'])} to {_fmt(r['HR_high'])}) | {_fmt(r['logrank_p'])} | {_fmt(r['untreated_logrank_p'])} |")
        lines += ["", "What I think: these are the numbers as they came out. With few organoids and few treated patients the intervals are wide; "
                  "look at the number of events before reading anything into a p-value. Full tables: `organoid_only/organoid_only_summary.tsv`.", ""]
    else:
        lines += ["No analysis was possible. Reasons:", ""]
        for f in (run / "organoid_only" / "run_log.tsv", run / "organoid_only" / "drugs_not_analysed.tsv"):
            if f.exists():
                d = pd.read_csv(f, sep="\t")
                col = "skip" if "skip" in d else "reason"
                if col in d:
                    for _, r in d.iterrows():
                        lines.append(f"- {r.get('model_drug', r.get('drug', ''))} {r.get('cohort', '')}: {r[col]}")
        if codes.get("organoid_only", 0):
            lines.append(f"- The runner stopped with exit code {codes['organoid_only']} (see the console output).")
        lines.append("")
    j = run / "organoid_plus_patient" / "joint_summary.tsv"
    lines += ["## 2. Organoids + patients together", "",
              "How it works: for each cancer-drug pair the repo already studies (nar/03), your organoids replace the lab models. Three models are "
              "compared on patient cohorts that were held out from training: organoid-only, patient-only, and joint (organoids + other patient "
              "cohorts). C = concordance (0.5 = no better than chance, 1 = perfect).", ""]
    if j.exists():
        t = pd.read_csv(j, sep="\t")
        t = t[t["cohort"] == "POOLED"] if (t["cohort"] == "POOLED").any() else t
        lines += ["| Unit | Features | Model | Outcome | n | Events | C | HR (95% CI) | p |", "|---|---|---|---|---|---|---|---|---|"]
        for _, r in t.iterrows():
            hr = f"{_fmt(r.get('HR'))} ({_fmt(r.get('HR_low'))} to {_fmt(r.get('HR_high'))})" if "HR" in t and not pd.isna(r.get("HR")) else "n/a"
            lines.append(f"| {r['unit']} | {r['features']} | {r['model']} | {r['outcome']} | {_fmt(r['n'])} | {_fmt(r['events'])} | {_fmt(r['C'])} | {hr} | {_fmt(r.get('p'))} |")
        lines += ["", "What I think: compare the three models within one unit; the joint model is only useful if it beats organoid-only and patient-only on held-out cohorts. "
                  "Per-cohort tables: `organoid_plus_patient/tables/`.", ""]
    else:
        lines += ["No unit could be run."]
        f = run / "organoid_plus_patient" / "units_not_run.tsv"
        if f.exists():
            for _, r in pd.read_csv(f, sep="\t").iterrows():
                lines.append(f"- {r.iloc[0]}: {r.iloc[1]}")
        lines.append("")
    (run / "ANALYSIS_SUMMARY.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
