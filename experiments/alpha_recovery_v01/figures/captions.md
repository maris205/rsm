# Figure captions and provenance

These figures display synthetic recovery and deterministic response geometry at the existing observing design. They are not new fits to measured astronomical or clock means. Reproduce both PNG and vector PDF files with:

```bash
python riemann_model/experiments/alpha_recovery_v01/code/plot_recovery.py
```

## Recovery capability

![Synthetic amplitude recovery and drift precision](recovery_capability.png)

**Figure 1. Amplitude recovery under the specified current observing design.** (a) Rejection rates for a fixed inverse-logarithmic-square response, using a two-sided 5% test of zero amplitude. Positive and negative injected amplitudes are plotted separately, with a small horizontal displacement for visibility. The null simulation has a 4.795% false-positive rate. Error bars are approximate 95% Monte Carlo uncertainties from 20,000 evaluation repetitions; the smooth line is the exact Gaussian power function under the assumed design. The horizontal coordinate is the injected present-day drift in units of the adopted clock standard error, 2.5 × 10⁻¹⁹ yr⁻¹. (b) Expected drift errors after fitting the same two instrument offsets where applicable. The astronomical-only error is 1.8426 × 10⁻¹⁶ yr⁻¹, whereas the joint error is 2.4999977 × 10⁻¹⁹ yr⁻¹; the clock accounts for 99.9998159% of the joint amplitude information. These are conditional statistical errors within the response and noise assumptions. Detecting a nonzero amplitude does not identify its time dependence.

Sources: `../results/recovery.csv` (main design `king_clock_offsets`, injected `s=2`, plus the null); `../results/summary.json` (design precision and thresholds). Original data means are not plotted. [Vector PDF](recovery_capability.pdf).

## Shape identifiability

![AIC selection frequencies and ideal precision budget](shape_identifiability.png)

**Figure 2. Nonzero-amplitude recovery does not separate the three logarithmic histories.** (a) AIC selection frequencies for injected exponents 1, 2 and 3, each with positive drift 1.25 × 10⁻¹⁸ yr⁻¹ (five adopted clock standard errors). Each row uses 20,000 repetitions of the King-plus-clock design with two free instrument offsets. The same noise realizations are reused across injections, as specified in the protocol. Near-collinear histories give essentially identical selection frequencies for all three true exponents; the small selection rate for the middle exponent does not establish that it is physically disfavored. Winner frequencies are not evidence probabilities, and an AIC selection is not an amplitude-detection test. The five-standard-error injection is a synthetic stress case, not a claim of compatibility with current measured clock means. (b) The minimum noise-free Δχ² against either adjacent exponent after refitting amplitude and offsets, for an injected exponent of 2. All astronomical standard errors, including the assumed extra scatter, are ideally scaled by the horizontal factor; clock error, sampling, background and reference scale remain fixed. Dots show the crossings of reference residual levels 1 and 9, which are not discovery significances. The reference scale is the Planck time; profiling that scale can further weaken identifiability. The curve is a conditional information budget, not a telescope performance forecast or a claim about a ten-year schedule.

Sources: `../results/recovery.csv` (main design, positive `d=12.5`); `../results/summary.json` (ideal precision curve and crossings); `../results/shape_geometry.csv` (independent check of each curve's present-precision starting value). [Vector PDF](shape_identifiability.pdf).
