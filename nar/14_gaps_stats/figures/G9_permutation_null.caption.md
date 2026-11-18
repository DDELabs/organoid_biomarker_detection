**Figure G9. How many small p-values would the 57-analysis survey give by chance, taking shared patients into account?**
We shuffled survival (time and death status together) between patients inside each TCGA cohort, 5,000 times. Treated
and untreated patients were shuffled separately. Colon, rectum and colon+rectum analyses share patients, so they were
shuffled together, with one shuffle per cohort used for every analysis. The organoid scores were not changed. This
keeps the overlap between analyses but breaks any real link between score and survival.
(a) Treated patients: blue bars = how often each count of p < 0.05 appears under shuffling; black dots = the simple
binomial count that assumes independent tests; red line = observed 7. Shuffle p = 0.030.
(b) The same for untreated patients: observed 12, shuffle p = 0.0006.
(c) Smallest of the 57 treated p-values under shuffling; observed 0.0011, shuffle p = 0.050
(this is the family-wise adjusted p of the best single analysis).
