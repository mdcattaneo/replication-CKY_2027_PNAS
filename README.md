# Replication: Cattaneo, Klusowski and Yu (2027, PNAS)

Replication files for simulation and empirical results.

## Replication

Requires R, make, a C++ compiler, Cairo support, Liberation Sans fonts, and `sha256sum`.
Run from the repository folder:

```sh
Rscript CKY_2027_PNAS--simuls.R
Rscript CKY_2027_PNAS--empapp.R
Rscript CKY_2027_PNAS--plots.R
Rscript src/checks.R test
Rscript src/checks.R check
```

Each design defaults to 1,000,000 replications. Override individual counts with
`CKY_M_REPS_P1`, `CKY_M_REPS_P2`, `CKY_M_REPS_JTPA_P1`, and
`CKY_M_REPS_JTPA_P2` (integers of at least 2). Set `CKY_OUTPUT_DIR` and
`CKY_PLOT_DIR` for separate short or custom runs; use the same output directory
for both runners, plotting, and output validation. Plotting and output validation
require results for all four designs.

Figure axes label the covariates `X` for one-covariate designs and `X_1`,
`X_2` for two-covariate designs; the SI defines their correspondence to
observations `x_i` and `x_{i,1}`, `x_{i,2}`. In the JTPA tables, `y` and `d`
denote outcomes and treatment indicators; `X1` and `X2` denote preprogram
earnings and age. Earnings axes and RMSE are displayed in thousands of
dollars, while age is displayed in years.

## Reference

Cattaneo, Klusowski and Yu (2027): [Accuracy Limits of Causal Trees for Individualized Treatment Effects](https://mdcattaneo.github.io/papers/Cattaneo-Klusowski-Yu_2027_PNAS.pdf).<br>
PNAS manuscript.<br>
[Supplemental Appendix](https://mdcattaneo.github.io/papers/Cattaneo-Klusowski-Yu_2027_PNAS--Supplement.pdf)<br>
DOI: 
