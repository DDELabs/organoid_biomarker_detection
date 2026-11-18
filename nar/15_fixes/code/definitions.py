"""Shared patient definitions (audit fixes I01, I02). Imported by nar/03 units.py and nar/05 b1_external.py.

I01  GSE39582 "5-FU-treated" (adjuvant):
     treated  = arm in {5FU, FUFOL, FOLFOX}  (training-table drugs FLUOROURACIL, FLUOROURACIL;LEUCOVORIN,
                FLUOROURACIL;LEUCOVORIN;OXALIPLATIN), stage II or III only.
     excluded = 'chemo (type not given)' / 'other' / 'unknown' arms (drug not recorded), FOLFIRI (irinotecan,
                stated separately as its own subset), stage 0, I and IV, missing stage.
     comparator = 'no chemotherapy', stage II or III only.
I02  FOLFIRINOX (oxaliplatin AND irinotecan in the same record) is removed from both arms of every
     oxaliplatin-vs-irinotecan comparison and from the oxaliplatin and irinotecan analysis units
     (as nar/09_registered_test does). In 5-FU-backbone analyses they stay (they did receive 5-FU).
"""
GSE39582_STAGES = (2, 3)
# regexes on the training-table `drugs` column (exact match, so CHEMOTHERAPY and FOLFIRI never match)
GSE39582_DRUGS = {
    "5FU_based": r"^FLUOROURACIL(?:;LEUCOVORIN(?:;OXALIPLATIN)?)?$",     # primary: 5FU + FUFOL + FOLFOX
    "5FU_alone": r"^FLUOROURACIL$",
    "FUFOL": r"^FLUOROURACIL;LEUCOVORIN$",
    "FU_no_oxaliplatin": r"^FLUOROURACIL(?:;LEUCOVORIN)?$",            # 5FU + FUFOL
    "FOLFOX": r"^FLUOROURACIL;LEUCOVORIN;OXALIPLATIN$",
    "FOLFIRI": r"^FLUOROURACIL;LEUCOVORIN;IRINOTECAN$",               # optional subset, never in the primary
    "untreated": r"^NONE$",
}
# GSE39582 `arm` labels (used by nar/05 which matches on 'arm|drugs')
GSE39582_ARMS = {"5FU_based": ("5FU", "FUFOL", "FOLFOX"), "5FU_alone": ("5FU",), "FUFOL": ("FUFOL",),
                 "FU_no_oxaliplatin": ("5FU", "FUFOL"), "FOLFOX": ("FOLFOX",), "FOLFIRI": ("FOLFIRI",),
                 "untreated": ("no chemotherapy",)}
BOTH_OX_IRI = r"OXALIPLATIN.*IRINOTECAN|IRINOTECAN.*OXALIPLATIN"       # FOLFIRINOX-type records (I02)
