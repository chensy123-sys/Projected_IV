# Projected Instrumental Variables

R code for PIV estimation, simulation studies, and the Australian Longitudinal Survey (ALS) application.

## Setup

Run all commands from the project root. Install the required packages once:

```bash
Rscript -e 'install.packages(c("dplyr", "tibble", "mgcv", "grf", "randomForest", "AER", "DoubleML", "mlr3", "mlr3learners", "ranger", "sisVIVE", "RobustIV", "tidyverse", "mice", "FactoMineR", "missMDA", "ggplot2", "purrr"), repos = "https://cloud.r-project.org")'
mkdir -p logs results results/fig results/tex
```

## Scripts

| Script | Purpose | Terminal command |
| --- | --- | --- |
| `info.R` | Set the number of replications (`500`), true ATE (`2`), and number of instruments (`3`). Automatically loaded by `monteCarlo.R`. | `Rscript script/info.R` (no printed output) |
| `monteCarlo.R` | Run simulations comparing PIV, GE, 2SLS, NUC, sisVIVE, and TSHT. Arguments: DGP (`1`–`6`) and sample size. | `Rscript script/monteCarlo.R 1 5000` |
| `summary.R` | Create LaTeX tables of bias, SD, estimated SE, and 95% coverage. | `Rscript script/summary.R` |
| `violin.R` | Create violin plots of estimation errors. | `Rscript script/violin.R` |
| `test_invalid_IV.R` | Generate the linear and quadratic toy-example figures. | `Rscript script/test_invalid_IV.R` |
| `data_preprocessing_new.R` | Prepare and impute the 1984 and 1985 ALS data. | `Rscript script/data_preprocessing_new.R` |
| `ALS_analysis.R` | Estimate the effect of more than 12 years of education on hourly wages. | `Rscript script/ALS_analysis.R` |

## Simulation workflow

Edit `script/info.R` to change the number of replications; keep `P = 3`. Run the six DGPs at each sample size:

```bash
for n in 2000 5000 10000; do
  for case in 1 2 3 4 5 6; do
    Rscript script/monteCarlo.R "$case" "$n" > "logs/DGP${case}_${n}.log" 2>&1
  done
done
```

The simulation saves `results/<n>_3/RES<case>.rds`. The reporting scripts currently expect different filenames, so prepare copies before running them:

```bash
for n in 2000 5000 10000; do
  mkdir -p "results/$n"
  for case in 1 2 3 4 5 6; do
    cp "results/${n}_3/RES${case}.rds" "results/${n}/RES${case}_5.rds"
  done
done
Rscript script/summary.R
Rscript script/violin.R
```

Tables are saved in `results/tex/`, and violin plots in `results/fig/`. The toy-example script separately saves `results/linear_model.pdf` and `results/quadratic_model.pdf`.

## ALS workflow

Place the raw files `00377-p.csv` and `00413-p.csv` in `results/data_private/`, then run:

```bash
Rscript script/data_preprocessing_new.R
cp ALS_DATA.RData results/ALS_DATA.RData
Rscript script/ALS_analysis.R
```

The copy step connects the preprocessing output to the analysis input. Results are saved as `results/ALS_results_table.tex` and `results/ALS_estimates_CI.pdf`.
