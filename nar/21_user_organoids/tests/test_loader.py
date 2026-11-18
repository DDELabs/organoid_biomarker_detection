"""Tests for the user-organoid loader.  Run:  python3 nar/21_user_organoids/tests/test_loader.py

Main test: data/curated/preclinical/colorectal_vandewetering2015 is written out in three different "user" formats
(synthetic folders under tests/_work/), imported, and the pathway scores must equal the scores computed directly from
the original curated set. Then small unit tests for orientation, Ensembl ids, counts vs log vs TPM, wide vs long, etc.
"""
import shutil
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code"))
import load_user_organoids as L  # noqa: E402
from obd.curated import PRECLIN  # noqa: E402
from obd.reference import ensembl_to_symbol  # noqa: E402

WORK = HERE / "_work"
SRC = PRECLIN / "colorectal_vandewetering2015"
ORIG_EXPR = pd.read_csv(SRC / "expression.tsv.gz", sep="\t", index_col=0)
ORIG_RESP = pd.read_csv(SRC / "response.tsv", sep="\t")
INV = {}
for _k, _v in ensembl_to_symbol().items():
    INV.setdefault(_v, _k)
RESULTS = {}


def fresh(name):
    d = WORK / name
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    return d


def reference_scores(expr):
    return L.score_pathways(expr)


def make_format_a():
    """Symbols, genes in rows, log values, csv; long drug table with explicit metric column."""
    d = fresh("format_A_symbols_rows_log_long")
    ORIG_EXPR.to_csv(d / "expression.csv")
    r = ORIG_RESP.rename(columns={"response": "value"})[["sample", "drug", "metric", "value"]]
    r.to_csv(d / "drug_response.csv", index=False)
    return d


def make_format_b():
    """Ensembl ids (with version suffix), SAMPLES in rows and GENES in columns, gzipped tsv; wide drug table (AUC in the name)."""
    d = fresh("format_B_ensembl_samplesrows_wide")
    keep = [g for g in ORIG_EXPR.index if g in INV]
    e = ORIG_EXPR.loc[keep]
    e.index = [INV[g] + ".7" for g in keep]
    e.T.to_csv(d / "expression.tsv.gz", sep="\t", index_label="sample")
    wide = ORIG_RESP.groupby(["sample", "drug"])["response"].median().unstack()
    wide.to_csv(d / "drug_response_AUC.tsv", sep="\t", index_label="sample")
    pd.DataFrame({"sample": ORIG_EXPR.columns, "cancer_type": "colorectal", "patient": ORIG_EXPR.columns}).to_csv(d / "samples.tsv", sep="\t", index=False)
    return d, ORIG_EXPR.loc[keep]


def make_format_c():
    """TPM (not log), symbols in rows, xlsx; drug table wide with DRUGS in rows and samples in columns (csv)."""
    d = fresh("format_C_tpm_xlsx_drugsrows")
    (2.0 ** ORIG_EXPR - 1).to_excel(d / "expression.xlsx", index_label="gene")
    wide = ORIG_RESP.groupby(["sample", "drug"])["response"].median().unstack()
    wide.T.to_csv(d / "drug_response_AUC.csv", index_label="drug")
    return d


def run(d, name, **kw):
    return L.import_folder(d, name, out_dir=WORK / "out" / name, **kw)


class SyntheticReproduction(unittest.TestCase):
    def check(self, res, ref_expr, tag):
        self.assertEqual(res["status"], "ok", res["report"])
        ref = reference_scores(ref_expr)
        got = pd.read_csv(res["out"] / "pathway_scores.tsv.gz", sep="\t", index_col=0)
        self.assertEqual(list(got.columns), list(ref.columns))
        self.assertEqual(sorted(got.index), sorted(ref.index))
        diff = float((got.loc[ref.index] - ref).abs().to_numpy().max())
        RESULTS[tag] = (ref.shape, diff)
        self.assertLess(diff, 1e-9, f"{tag}: max abs difference {diff}")
        # drug response must equal the original median AUC per sample and drug
        resp = pd.read_csv(res["out"] / "response.tsv", sep="\t").groupby(["sample", "drug"])["response"].median()
        orig = ORIG_RESP.assign(drug=ORIG_RESP["drug"].map(L.common_drug_name)).groupby(["sample", "drug"])["response"].median()
        both = pd.concat([resp, orig], axis=1, keys=["got", "orig"])
        self.assertEqual(int(both["got"].isna().sum() + both["orig"].isna().sum()), 0, tag)
        self.assertLess(float((both["got"] - both["orig"]).abs().max()), 1e-9)

    def test_a_symbols_rows_log_long(self):
        res = run(make_format_a(), "synthetic_A")
        self.assertIn("genes in rows", res["einfo"]["orientation"])
        self.check(res, ORIG_EXPR, "A symbols / genes in rows / log / long drug table")
        # this curated set has 10,912 genes, so 533 of the 673 Reactome pathways have >= 5 genes in it
        patients = pd.read_csv(L.REPO / "robust/results/training_table/pathway_scores_geo.tsv.gz", sep="\t", index_col=0, nrows=2)
        self.assertEqual(res["n_pathways"], 533)
        self.assertEqual(patients.shape[1], 673)
        self.assertTrue(set(res["scores"].columns) <= set(patients.columns))

    def test_b_ensembl_samples_rows_wide(self):
        d, ref = make_format_b()
        res = run(d, "synthetic_B")
        self.assertEqual(res["einfo"]["gene_id_type"], "Ensembl gene ids")
        self.assertIn("flipped", res["einfo"]["orientation"])
        self.assertEqual(res["rinfo"]["layout"], "wide")
        self.assertEqual(res["rinfo"]["metric"], "AUC")
        self.check(res, ref, f"B Ensembl ids / samples in rows / log / wide ({len(ref)} of {len(ORIG_EXPR)} genes have an Ensembl id)")

    def test_c_tpm_xlsx_drugs_rows(self):
        res = run(make_format_c(), "synthetic_C")
        self.assertEqual(res["einfo"]["value_type"], "tpm")
        self.assertEqual(res["rinfo"]["layout"], "wide")
        self.check(res, ORIG_EXPR, "C TPM / xlsx / drugs in rows")


class Units(unittest.TestCase):
    def table(self, genes, samples, scale="log", seed=0):
        rng = np.random.default_rng(seed)
        x = rng.normal(7, 2, (len(genes), len(samples)))
        if scale == "counts":
            x = np.round(2 ** x * 20)
        elif scale == "tpm":
            x = 2 ** x
        return pd.DataFrame(x, index=genes, columns=samples)

    def folder(self, name, expr, resp, expr_name="expression.tsv", resp_name="drug_response.tsv", index_label="gene"):
        d = fresh(name)
        expr.to_csv(d / expr_name, sep="\t" if expr_name.endswith("tsv") else ",", index_label=index_label)
        resp.to_csv(d / resp_name, sep="\t" if resp_name.endswith("tsv") else ",", index=resp.index.name is not None)
        return d

    def base_genes(self):
        return list(ORIG_EXPR.index[:1500])

    def samples(self):
        return [f"S{i}" for i in range(10)]

    def long_resp(self, metric="AUC"):
        rng = np.random.default_rng(1)
        rows = [(s, dn, metric, rng.uniform(0.2, 0.9)) for s in self.samples() for dn in ("5-FU", "Oxaliplatin", "Irinotecan")]
        return pd.DataFrame(rows, columns=["sample", "drug", "metric", "value"])

    def test_value_type_detection(self):
        g, s = self.base_genes(), self.samples()
        for scale, want in (("log", "log"), ("counts", "counts"), ("tpm", "tpm")):
            rep = L.Report()
            _, info = L.parse_expression(self.table(g, s, scale).reset_index().rename(columns={"index": "gene"}), rep)
            self.assertEqual(info["value_type"], want, scale)

    def test_counts_give_same_scores_as_log(self):
        g, s = self.base_genes(), self.samples()
        log = self.table(g, s, "log")
        counts = np.round(2 ** log * 50)         # counts whose log2 ranks are ~ those of the log table
        a = L.score_pathways(L.parse_expression(counts.reset_index().rename(columns={"index": "gene"}), L.Report())[0])
        b = L.score_pathways(np.log2(counts / counts.sum() * 1e6 + 1))
        self.assertLess(float((a - b).abs().to_numpy().max()), 1e-9)

    def test_orientation(self):
        g, s = self.base_genes(), self.samples()
        e = self.table(g, s)
        for flip, text in ((False, "genes in rows"), (True, "samples in rows")):
            t = e.T if flip else e
            df = t.reset_index().rename(columns={"index": "x"})
            _, info = L.parse_expression(df, L.Report())
            self.assertIn(text, info["orientation"])

    def test_ensembl_ids_with_and_without_version(self):
        g = [x for x in self.base_genes() if x in INV][:800]
        e = self.table(g, self.samples())
        for suffix in ("", ".12"):
            e2 = e.copy()
            e2.index = [INV[x] + suffix for x in g]
            out, info = L.parse_expression(e2.reset_index().rename(columns={"index": "gene"}), L.Report())
            self.assertEqual(info["gene_id_type"], "Ensembl gene ids")
            self.assertEqual(sorted(out.index), sorted(set(g)))

    def test_entrez_ids_stop_with_clear_message(self):
        e = self.table([str(i) for i in range(1000, 1500)], self.samples())
        with self.assertRaises(L.ImportProblem):
            L.parse_expression(e.reset_index().rename(columns={"index": "gene"}), L.Report())

    def test_long_and_wide_response_agree(self):
        long = self.long_resp()
        wide = long.pivot(index="sample", columns="drug", values="value")
        s = self.samples()
        a, _ = L.parse_response(long, L.Report(), s)
        w = wide.reset_index()
        b, _ = L.parse_response(w, L.Report(), s, metric_hint="AUC")
        a = a.set_index(["sample", "drug"])["value"].sort_index()
        b = b.set_index(["sample", "drug"])["value"].sort_index()
        pd.testing.assert_series_equal(a, b)
        self.assertEqual(sorted(set(a.index.get_level_values("drug"))), ["FLUOROURACIL", "IRINOTECAN", "OXALIPLATIN"])

    def test_ic50_is_logged_and_direction(self):
        long = self.long_resp("IC50")
        res, info = L.parse_response(long, L.Report(), self.samples())
        self.assertEqual(info["metric"], "ln_IC50")
        self.assertTrue(np.allclose(np.exp(res["value"]), res["value_raw"]))
        self.assertEqual(info["direction"], "lower")
        _, info2 = L.parse_response(self.long_resp("inhibition"), L.Report(), self.samples())
        self.assertEqual(info2["direction"], "higher")

    def test_unknown_metric_warns_and_assumes_lower(self):
        rep = L.Report()
        _, info = L.parse_response(self.long_resp("mystery"), rep, self.samples())
        self.assertEqual(info["direction"], "lower")
        self.assertTrue(any("LOWER or a HIGHER" in w for w in rep.warn))

    def test_sample_name_mismatch_is_reported_and_fuzzy_matched(self):
        rep = L.Report()
        m, miss = L.reconcile(["P19TA", "P5", "P6"], ["p-19 ta", "P5", "Z9"], "drug-response", rep)
        self.assertEqual(m["p-19 ta"], "P19TA")
        self.assertEqual(miss, ["Z9"])
        self.assertTrue(rep.warn and rep.did)

    def test_no_overlap_stops_with_plain_message(self):
        g = self.base_genes()
        e = self.table(g, self.samples())
        r = self.long_resp()
        r["sample"] = "other_" + r["sample"]
        d = self.folder("no_overlap", e, r.set_index("sample"))
        res = L.import_folder(d, "x", out_dir=WORK / "out" / "no_overlap")
        self.assertEqual(res["status"], "stopped")
        self.assertIn("do not match", (WORK / "out" / "no_overlap" / "IMPORT_REPORT.md").read_text())


if __name__ == "__main__":
    shutil.rmtree(WORK / "out", ignore_errors=True)
    prog = unittest.main(exit=False, verbosity=2)
    print("\nPathway-score reproduction (max absolute difference to scores computed from the original curated set):")
    for k, (shape, diff) in RESULTS.items():
        print(f"  {k}: {shape[0]} samples x {shape[1]} pathways, max |difference| = {diff:.2e}")
    sys.exit(0 if prog.result.wasSuccessful() else 1)
