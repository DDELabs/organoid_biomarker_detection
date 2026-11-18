# Fig. 1b (redrawn): how "network proximity" is measured

a, b. A toy protein network (grey dots = proteins, lines = known interactions). Red dots = the drug's target
proteins. Blue dots = the genes of one pathway. For each pathway gene we count the steps to the nearest target and
average them (closest distance d). In (a) the pathway sits next to the targets; in (b) it sits far away.
c. To know whether d is small, the same distance is computed for 1,000 random gene sets whose genes have the same
number of network links (degree-matched). The z-score compares the real d with these random ones. A pathway is
called "proximal" when z <= -1.2816 (the lowest 10% of the random sets), the paper's cut-off.
Illustration only: the network and the histogram are simulated drawings, not data.
