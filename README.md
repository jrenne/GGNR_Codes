# The Shadow-Rate Model: Let's Make it Real

Replication package for the U.S. results in Golinski, Guilloux-Nefussi, and Renne, *The Shadow-Rate Model: Let's Make it Real*.

## Contents

- `matlab/`: Matlab code that evaluates the estimated model and reproduces the main U.S. tables and figures.
- `data/data_JPR.mat`: frozen Matlab input used for the paper results.
- `data_pipeline/`: R scripts that rebuild the U.S. public-data database and export a Matlab-ready input file.
- `Figures/`: output folder for replicated figures.
- `Tables/`: output folder for replicated LaTeX tables.

## Quick Replication

Open Matlab from the repository root and run:

```matlab
run('matlab/run_replication.m')
```

The script loads `data/data_JPR.mat`, evaluates the model at the parameter values used in the paper, and writes outputs to `Figures/` and `Tables/`.

## Rebuilding The Data

The R pipeline downloads public data and reconstructs the U.S. estimation database:

```bash
Rscript data_pipeline/build_us_database.R
```

This creates:

- `data/Data_US/data.Rda`
- `data/Data_US/Estimation_data_US.xlsx`
- `data/data_from_r.mat`

After building `data/Data_US/data.Rda`, the supplementary normality table can be regenerated with:

```bash
Rscript data_pipeline/03_make_real_yield_normality_table.R
```

The pipeline uses data from FRED, the Federal Reserve Board, the Philadelphia Fed Survey of Professional Forecasters, the New York Fed, and the FRB/US data package. FRED access requires an API key. Set it before running R:

```bash
export FRED_API_KEY="your_fred_key"
```

The regenerated `data/data_from_r.mat` contains the core model inputs constructed from the public database. The frozen file `data/data_JPR.mat` remains the exact input used for the paper replication; it also contains comparison term-premium series used in the term-premium figure.

### Replacing The Frozen Matlab Input

To run the Matlab replication with a newly rebuilt dataset rather than the frozen paper input:

1. Rebuild the data:

   ```bash
   Rscript data_pipeline/build_us_database.R
   ```

2. Keep a backup of the frozen paper input:

   ```bash
   cp data/data_JPR.mat data/data_JPR_frozen_paper.mat
   ```

3. Replace the Matlab input with the newly generated file:

   ```bash
   cp data/data_from_r.mat data/data_JPR.mat
   ```

4. Re-run the Matlab replication:

   ```matlab
   run('matlab/run_replication.m')
   ```

This replacement should be interpreted as an updated-data run, not as the exact paper replication. The frozen `data/data_JPR.mat` contains 670 observations and includes comparison term-premium series used in the term-premium figure. The regenerated file fills those comparison term-premium series with missing values unless they are added separately. Public macro series are also revised over time, so a newly rebuilt file need not reproduce the frozen paper input exactly.

### Data Consistency Check

The frozen paper input uses the monthly sample from October 1968 through July 2024, for 670 observations. The R pipeline is configured to rebuild that same sample. To compare the rebuilt Matlab input with the frozen paper input, run:

```bash
Rscript data_pipeline/05_compare_matlab_inputs.R
```

After aligning the dates and variable ordering, the nominal yields, real yields, and survey-expectation series reproduce the frozen file to numerical precision. The frozen `ffr` series is exactly the FRED `FEDFUNDS` monthly series over the same sample; if the exported `ffr` is missing, re-run the full R download/preparation pipeline so that `FEDFUNDS` is included in `data/Data_US/Estimation_data_US.xlsx`.

The remaining differences are in the macro block: output growth (`dy`), inflation (`Pi`), and the output gap (`z`). These series depend on public macro data that are revised over time, including the Brave-Butters-Kelley monthly GDP series, CPI, and GDP/GDPPOT. The Taylor-rule policy rate (`PTR`) reproduces the frozen input to numerical precision. For exact replication of the submitted paper, use the frozen `data/data_JPR.mat`; use `data/data_from_r.mat` when intentionally updating the data vintage or checking the data-construction pipeline.

## R Package Dependencies

The data pipeline uses:

- `fredr`
- `Hmisc`
- `readxl`
- `openxlsx`
- `dplyr`
- `zoo`
- `R.matlab`

Install missing packages with, for example:

```r
install.packages(c("fredr", "Hmisc", "readxl", "openxlsx", "dplyr", "zoo", "R.matlab"))
```

## Matlab Notes

The main entry point is `matlab/run_replication.m`. It is intentionally an evaluation/replication script, not an estimation-search script. The manual likelihood-search loop used during development has been removed from the public entry point.

The code was written for Matlab and uses standard plotting, optimization, and matrix routines. Figures are saved as PDF files.

## License

The original code in this repository is released under the MIT License; see `LICENSE`. Data downloaded by the R scripts remain subject to the terms of their original providers.
