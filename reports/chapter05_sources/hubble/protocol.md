# Calibrated time-model comparison

Continuation authorized by the user. Written before fitting this extension.
This is exploratory development with existing 2022 compressed Cepheid/SN
calibration and archived DESI distance products, not a new 2026 H0 measurement
and not a refit of Cepheid photometry or the H0DN network.

Use the pinned author combined Pantheon+SH0ES selection: zHD>.01 OR calibrator,
1657 rows =77 calibrators+1580 cosmological rows. Calibrators use CEPH_DIST even
when zHD>.01. Retain the entire selected STAT+SYS covariance including cross
blocks. Predict a common absolute magnitude M for all rows. No extra Gaussian
SH0ES/H0DN H0 prior and no extra uncalibrated Pantheon likelihood is added.

Primary distance combination: FS14 plus calibrated Pantheon+SH0ES. Baseline
dimensionless clock q=67.4*(GYR_S/MPC_KM)*13.8. H0 is now fit and t0=q/H0 in
consistent units, hence t0 is a model-imposed late-time boundary, not a measured
Big-Bang age. Compare all three families with finite-future positive g, and
separately with the nondecreasing hypothesis at baseline q. Future-positive
sensitivity cases q/q0=.9 and1.1 are also fit. Fixed q is not an added datum.

Bounds: Omega_ref [.05,.60], s[-.60,.60], H0[50,90] km/s/Mpc, intersected with
the declared future-positive or nondecreasing domain. The logarithmic lower
domain endpoint depends on t0 through L0 and must be recomputed when H0 varies.
All fits explicitly search parameter-domain edges. Preserve boundary flags and
multiple local solutions; no automatic expansion or discovery interpretation.

For exponential/power law at fixed (Omega_ref,s,q), solve the two-column full
covariance GLS with columns common1 and cosmological-row indicator:
calibrator base=CEPH_DIST; cosmological base=5log10[(1+zHEL)I];
coefficients=(M,alpha), alpha=25+5log10(c/H0). This profiles H0 exactly, subject
to its broad bounds. The ruler S=c/(H0 rd) is separately profiled from BAO.
For logarithmic models, I also depends weakly on H0 through L0: explicitly
optimize H0, recompute the background, and profile only M and S at each trial.
The scalar H0 optimization includes both BAO and calibrated-SN contributions.

Use a normalized coordinate in each declared slope domain to handle the small
H0-dependence of the logarithmic future boundary without changing the domain.
This is an optimizer coordinate, not an extra parameter. Null model s=0 has
4 fit parameters (Omega_ref,H0,rd,M); a fixed-q evolving model has5.
Report Δchi2 relative to the null for the identical observations and ΔAIC=2−Δchi2.
Neither directional boundaries nor exploratory model/clock choices are
translated to Gaussian discovery claims.

For baseline-q future-positive cases, profile H0 with shape parameters
reoptimized at every trial. At fixed H0 profile only M, not a free second
intercept that would remove the calibration. Also obtain the null H0 profile.
Ranges are conditional likelihood level sets, not marginalized posteriors.

Independent benchmark: calibrated SN-only flat LCDM, plus calibrated FS14 null.
Reproduce the native likelihood with independent background quadrature and
linear nuisance algebra. Validate parameter counts, crosscovariance, all-edge
optimization, shape/H0 profile nesting, source hashes and tighter integration.
The inferred free rd is a late-time ruler value, not an early-universe sound
horizon calculation. No Planck likelihood, age likelihood or growth likelihood
is added, so a claimed resolution of Hubble tension remains outside scope.
