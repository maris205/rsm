# Threshold matching and low-energy running for the mediator chain

This bounded stage matches the selected heavy Gaussian sector onto the light modulus-scalar sector at leading external-source order S². The finite high-scale boundary is fixed. Changing renormalization scale moves the large logarithm into the local coefficient; it does not remove the physical force.

For the 2 TeV, 127-site, mG = 10⁻¹⁵ eV control, the new tree plus selected one-loop force is about **2.00 × 10⁻⁵** of the prescribed clock force. The one-loop/new-tree ratio remains 1.23447, while the actual light-sector repeated-loop parameters are of order 10⁻²³. Thus that ratio alone does not establish perturbative failure. This is not a complete UV or supergravity matching result.

The companion locality audit quantifies sensitivity to explicitly prescribed bypasses. A vectorlike Higgs construction realizes the mass matrix but leaves an unstabilized chiral direction in its minimal form and does not protect the existing bulk modulus's endpoint isolation.

- [Chinese results](reports/results_cn.md) and [figure PDF](figures/threshold_running_and_force.pdf)
- [Protocol and Material Passport](protocol.md)
- [Independent output comparison](reports/independent_results_cn.md)
- [Main summary](results/summary.json)
- [Threshold/RG derivation](../../reports/threshold_rg_derivation_20260926.md)
- [Low-energy power counting and scoped two-loop graph group](../../reports/threshold_power_counting_20260926.md)
- [Prescribed bypasses and concrete mass-sector origin](../../reports/chain_locality_breaking_20260926.md)
- [Material manifest](materials_manifest.json) and [archive audit](results/archive_audit.json)

From the repository root:

```bash
python experiments/threshold_rg_v01/code/run_thresholds.py
python experiments/threshold_rg_v01/code/independent_check.py
python experiments/threshold_rg_v01/code/plot_results.py
python experiments/threshold_rg_v01/code/audit_archive.py
```

The specialized reports link their own scripts and retained failures. Re-execution may change timestamps; use `--write-manifest` only for an intentional new material snapshot. Numerical and archive consistency are separate from physical validation.

No new observations or cosmic integration are included. Historical experiments and the 86-page manuscript are unchanged. Protected locality, complete threshold and supergravity matching, common vacuum rematching, and a Riemann origin for interactions or the inverse-log-squared law remain open.
