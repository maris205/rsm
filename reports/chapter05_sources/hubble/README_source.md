# Fixed-clock calibrated time models

This completed stage uses the native 1657-row Pantheon+SH0ES likelihood and
the 14-entry FS14 distance product. The calibration comes from the archived
2022 release, with full calibrator/flow covariance. No summary-H0 prior is
added. The three time-law families share fixed q0=0.9512450389523021 and
explicitly declared shape domains; q/q0=0.9 and 1.1 are sensitivity cases.

- `protocol.md`: pre-calculation analysis specification.
- `logarithmic.json`, `exponential.json`, `power_law.json`: four fits per
  family, competing optimizer solutions, baseline H0 profiles, convergence
  checks and source hashes.
- `comparison.csv`: all twelve conditional fits.
- `null_profile.json`: independent calculation of the common null H0 profile
  using the production likelihood.
- `summary_provenance.json`: fit-input, null-profile and table hashes.
- `validation.json`: independent background/QR likelihood and optimization
  checks; `passed` is true for the completed files.

The common null gives H0=73.7318091 km/s/Mpc and chi2=1467.2859634.
At q0, finite-future-positive logarithmic, exponential and power histories
improve chi2 by 0.841858, 4.944324 and 0.842242, respectively. The exponential
nondecreasing alternative improves it by 1.171197. H0 remains near 73;
the ruler near 136 Mpc is a free late-time fit, not an early-universe prediction.

See the [Chinese results report](../../papers/calibrated_models_2026-09-19/results.md),
[likelihood derivation](../../papers/calibrated_models_2026-09-19/likelihood_notes.md),
and [calibrated scalar realization](../calibrated_scalar_forward_2026-09-19/README.md).
The completed manuscript stage was preserved before the wider clock study in
[the revision snapshot](../../papers/revision_2026-09-19/calibrated_stage/README.md).

The follow-up clock-profile stage is stored separately in
`../clock_profile_2026-09-19/`; its output does not replace these conditional
results. Reproduction commands are in the project and code README files.
