# B3 progress
Code complete in code/. Running (relaunch if dead; all steps resume/skip finished outputs):
1. `cd code; python3 bench.py 0 1` (organoid benchmark, one pickle per pair in _data/bench/)
2. chain: `python3 cells.py cv; python3 cells.py gcn; python3 cells.py lc; python3 bench.py 0 1 --reverse; python3 learning.py` (writes _data/chain2.done)
Then: `python3 tables.py; python3 figures.py all; python3 suppfig7.py` (suppfig7 already done) and write REPORT.md.
