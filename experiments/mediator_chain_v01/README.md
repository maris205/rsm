# Gapped mediator chain: endpoint suppression and finite auxiliary response

A declared nearest-neighbor chain of massive Stückelberg multiplets produces an exponentially small endpoint cross coupling from order-one local parameters. Its complete tree matching also produces an unsuppressed modulus self interaction. This changes the finite-auxiliary-field response and invalidates the previous worked scale as a slow-clock example of this particular chain.

The lower-scale controls pass selected static-response and mass-gap conditions. They do not establish a complete protected theory: theory-space locality, isolation of the existing bulk modulus at one endpoint, full vacuum matching, and threshold/RG treatment remain open. The fixed-scale rigid one-loop response can exceed the new static tree response because of a large logarithm, although both selected forces are small in one low-scale control.

- [Chinese results and figure](reports/results_cn.md)
- [Protocol and Material Passport](protocol.md)
- [Independent output comparison](reports/independent_comparison_cn.md)
- [Main summary](results/summary.json) and [conditional response thresholds](results/response_limits.csv)
- [Matching and response figure, PDF](figures/chain_matching_and_response.pdf)
- [Primary sources and protection limits](../../reports/mediator_chain_sources_20260926.md)
- [Independent endpoint algebra and unbroken spectrum](../../reports/mediator_chain_algebra_20260926.md)
- [Finite auxiliary response and selected rigid Gaussian spectrum](../../reports/mediator_auxiliary_response_20260926.md)
- [Full local supergravity response and modulus Hessian](../../reports/mediator_sugra_response_20260926.md)
- [Material manifest](materials_manifest.json) and [archive consistency audit](results/archive_audit.json)

From the repository root:

```bash
python experiments/mediator_chain_v01/code/run_chain.py
python experiments/mediator_chain_v01/code/independent_compare.py
python experiments/mediator_chain_v01/code/plot_results.py
python experiments/mediator_chain_v01/code/audit_archive.py
```

The source-specific scripts linked from the four reports reproduce their own scoped checks. The auxiliary-response script deliberately returns a failure status for the retained old-scale second-order approximation accuracy test: its 0.2478% slope error exceeds the declared 0.2% goal. The old benchmark uses the exact spectrum; the smaller-splitting convergence test passes the same accuracy goal. Do not erase this failure or replace it by a looser tolerance.

Re-execution can change timestamps and therefore manifest hashes. Use `--write-manifest` only when deliberately creating a new snapshot. Material consistency is separate from scientific acceptance: an archive can correctly preserve an unsuccessful approximation or candidate.

No observations or cosmic integration are added. Historical experiments and the 86-page manuscript are unchanged. A Riemann origin of real interactions and the inverse-log-squared law remains an unestablished part of the research idea.
