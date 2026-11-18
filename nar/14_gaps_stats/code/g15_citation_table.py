"""G15: literature table from PubMed records saved in _data/pubmed/ (citation APIs not reachable from this machine:
api.openalex.org, api.crossref.org and api.semanticscholar.org were refused by the network proxy on 2026-11-18).
PubMed has no 'cited by' link here, so this is a related-articles + targeted-keyword substitute, not a citation list.
The 'relation' column is my reading of each title/abstract."""
import json

import pandas as pd

from g_common import DATA, TAB

NOTES = {
    "35764641": "same group (Kong, Kim); extends the network-proximity + ML idea to immunotherapy response; no re-test of the organoid biomarkers",
    "41999762": "same group (Kim S lab, POSTECH); network-propagation deep learning on mutations, 10 TCGA cancer-drug pairs; different method, does not re-test the organoid biomarkers",
    "34154646": "systems pharmacogenomics with network methods (medulloblastoma); related approach, no re-test",
    "33882314": "organoid targeted RNA-seq screening (colorectal); transcriptomic readout in organoids, no patient re-test of Kong biomarkers",
    "36675023": "new transcriptomic biomarkers of 5-FU resistance; does not test the BH3-only pathway marker",
    "35805961": "biomarkers for drug-treatment monitoring in bladder organoids; does not test the amino-acid pathway marker",
    "33709991": "blinded PDO chemotherapy-response study in stage IV colorectal cancer (functional avatar evidence)",
    "35696669": "PDO prediction of chemotherapy response in colorectal cancer (functional avatar evidence)",
    "37345888": "drug repurposing screen in colorectal organoids",
    "38767747": "ML prognosis/drug-susceptibility signature in bladder cancer (TCGA); unrelated method",
}


def main():
    rows = []
    for f, src in (("related_top20.json", "PubMed similar-articles to PMID 33127883 (top 20)"),
                   ("search_organoid_ml.json", "PubMed keyword searches (organoid / network-based / drug response, 2021-2026)")):
        d = json.load(open(DATA / "pubmed" / f))
        for a in d["articles"]:
            ids = a["identifiers"]
            pm = ids.get("pmid")
            rows.append({"pmid": pm, "doi": ids.get("doi") or a.get("doi"), "year": a["publication_date"].get("year"),
                         "journal": a["journal"]["title"], "title": a["title"], "found_by": src,
                         "relation_to_Kong2020": NOTES.get(pm, "not a re-test or extension (topic overlap only)")})
    extra = [{"pmid": "41999762", "doi": "10.1016/j.crmeth.2026.101411", "year": "2026", "journal": "Cell reports methods",
              "title": "A network-based deep learning model integrating subclonal architecture for therapy response prediction in cancer.",
              "found_by": "PubMed author/affiliation search (Kim S, POSTECH, network-based)", "relation_to_Kong2020": NOTES["41999762"]}]
    T = pd.DataFrame(rows + extra).drop_duplicates("pmid")
    T["retests_or_extends"] = T.relation_to_Kong2020.str.contains("same group|extends")
    T.sort_values(["retests_or_extends", "year"], ascending=[False, False]).to_csv(TAB / "G15_literature_related_to_Kong2020.tsv", sep="\t", index=False)
    print(T[T.relation_to_Kong2020 != "not a re-test or extension (topic overlap only)"][["pmid", "doi", "year", "relation_to_Kong2020"]].to_string())
    print(len(T))


if __name__ == "__main__":
    main()
