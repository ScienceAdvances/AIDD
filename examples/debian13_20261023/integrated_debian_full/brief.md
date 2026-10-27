# Case study brief

## Research question

How can structural and sequence models prioritize variants for further validation?

## Evidence

| Case | Metric | Value | Unit |
|---|---|---|---|
| HADDOCK Debian full | top_score_dockq | 0.569 | dimensionless |
| HADDOCK Debian full | best_output_dockq | 0.569 | dimensionless |
| Sequence optimization | WT_nll | 1.9176931381225586 | mean negative natural log probability |
| Sequence optimization | single_D5E_nll | 1.7470042705535889 | mean negative natural log probability |
| Sequence optimization | single_S4E_nll | 1.7673014402389526 | mean negative natural log probability |
| Sequence optimization | single_L9K_nll | 1.7776105403900146 | mean negative natural log probability |
| CAR-T | 4H_[B]-L3 | 30.3719875 | relative IFNg |
| CAR-T | 4H_[B]-H2L3wt | 12.65709325 | relative IFNg |
| CAR-T | S7B_[B]-L3 | 10695.333333333334 | MFI |
| CAR-T | S7B_[B]-H2L3wt | 2577.3333333333335 | MFI |
| CAR-T | S7I_5:1_[B]-L3 | 9.288652 | relative viability percent |
| CAR-T | S7I_5:1_[B]-H2L3wt | 80.13090333333334 | relative viability percent |
| Factor Xa MD | official_1ns_rmsd_to_em | 0.10297157412935323 | nm |
| Factor Xa MD | local_run_rmsd_to_em | 0.10307621940298507 | nm |
| GNINA | workflow_state | engine_run | state |

## Decision and next work

The sequence run contains 197 candidates and 4 nominated sequences, including WT.
Use candidate_validation_plan.csv to connect each nomination to a testable hypothesis.
Use evidence.csv to recover the measurement definition and source of each claim.
Describe Factor Xa MD and CAR-T reanalysis as separate case studies with their own biological contexts.
GNINA execution state follows its run manifest.
