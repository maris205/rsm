# Persistent restoration and correlated matter feedback

This bounded stage tests a declared four-dimensional modulus current that can preserve positive transverse clock curvature when the prescribed clock potential tends to zero. It also computes the correlated matter soft mass and the resulting selected quantum feedback. This is a conditional theoretical calculation, not a new observational fit or a complete UV model.

The worked example passes the stated finite-interval local conditions and remains locally stable after the selected charged loop correction. It requires an unexplained input coupling of approximately **2.93 × 10⁻⁶¹**. A persistent soft term leaves a quantum slope that eventually exceeds the prescribed exponentially decreasing clock force. Thus persistent transverse restoration alone does not establish permanent slow evolution.

- [Chinese results and interpretation](reports/results_cn.md)
- [Declared action, protocol and Material Passport](protocol.md)
- [Independent calculation and retained numerical failures](reports/independent_results_cn.md)
- [Conditional parameter window](figures/persistent_window.pdf)
- [Persistent mass and future feedback limit](figures/persistent_future_limit.pdf)
- [Complete local geometry](../../reports/persistent_current_geometry_20260926.md)
- [Selected charged Gaussian determinant](../../reports/persistent_charged_loop_20260926.md)
- [Primary sources and symmetry limits](../../reports/persistent_current_sources_20260926.md)
- [Implementation review and cancellation limits](../../reports/persistent_current_geometry_20260926_review.md)
- [Machine-readable summary](results/summary.json)
- [Material manifest](materials_manifest.json) and [archive audit](results/archive_audit.json)

From the repository root:

```bash
python experiments/persistent_restoration_v01/code/run_restoration.py
python experiments/persistent_restoration_v01/code/independent_check.py --derive --compare
python reports/persistent_current_geometry_20260926_checks.py
python reports/persistent_charged_loop_20260926_checks.py
python reports/persistent_current_sources_20260926_checks.py
python experiments/persistent_restoration_v01/code/plot_results.py
python experiments/persistent_restoration_v01/code/audit_archive.py
```

Re-execution may change recorded timestamps; preserve the historical manifest and use `--write-manifest` only when deliberately creating a new material snapshot. The implementation review JSON records a bounded internal review, not an executable independent experiment.

Only declared sectors and fixed finite matching are included. Tiny-coupling protection, complete heavy-spectrum loops, global rolling stability, real cosmic time, and a Riemann derivation of the inverse-log law remain open. The exact-cancellation branch has uncontrolled second-order terms and omits a known vacuum soft contribution in its formal far-tail kernel; its signed future root is not a full-model physical prediction. Historical experiments and the 86-page manuscript are unchanged.
