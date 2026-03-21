# Shift Identification Mechanism via Multiple Imputation, for Spatial Time-Series Repair

## Overview

This project implements a workflow for returning a sequence of displaced attribute values to their true time points, within spatially-correlated sensor data that exhibits irregular pattern due to timestamp error. The workflow has been specifically developed for application to multivariate air monitoring data, though may be generalizable to other forms of spatially-correlated time series too.

The computational framework relies on the generation of shifted data representations, and on performing multiple imputation via chained equations using penalised regression models (LASSO or Elastic Net) that depend on time series block validation. The mean absolute scaled error (MASE) relating to each temporal shift is pooled across imputations, and is used to identify the optimum data representation.

To determine whether or not the data has been shifted far enough, one-tailed significance testing of the minimum pooled imputation error is conducted by implementing Rubin's rules, and using Newey-West heteroskedasticity and autocorrelation consistent (HAC) standard errors to account for temporal dependence during estimation of within-imputation variance.

The framework further incorporates a repair mechanism that facilitates correction of potentially multiple (multivariate) series commonly affected by the same error. The pipeline ends with an evaluation of the repair by assessing RV data similarity versus a chosen reference location.

**Input Data Requirements**

The input excel file should:

-   Contain one sheet per univariate spatial set which...
    -   Includes a date/time column
    -   Includes a target series (spatial location) that contains a bounded subsequence displaying irregular pattern
    -   Contains one column per spatial location
-   Employ consistent structure across multiple sheets (if the spatial data is multivariate)

**Pipeline:**

-   Loads data for the chosen working variable (i.e. for one specified excel sheet only)
-   Prepares the working data frame
-   Computes required spatio-temporal window
-   Generates shifted data representations
-   Generates imputed datasets
-   Computes the pooled imputation error for each data shift
-   Performs statistical significance testing on the imputation error
-   Extracts the optimum shift
-   Sequentially repairs period of irregular pattern for each variable within the input data file (i.e. across all excel sheets)
-   Assesses improvement by computing RV similarity versus chosen reference location

------------------------------------------------------------------------

## Author

-   Name: Natalie Benschop
-   GitHub: <https://github.com/Natalie648>

------------------------------------------------------------------------

## Project Structure

``` text
simmi/
│
├── data/                         # Input dataset (Excel format)
├── R/                            # Custom function scripts
├── output/                       # Generated outputs
├── scripts/main.R                # Primary analysis script
├── simmi.Rproj                   # RStudio project file
├── LICENSE                       # MIT License
└── README.md                     # Project documentation
```

## Installation

**1. Install R and RStudio**

Download from: <https://cran.r-project.org/>

**2. Install required packages**

`install.packages(c( "readxl", "dplyr", "Metrics", "data.table", "mice", "glmnet", "sandwich", "lubridate", "writexl", "FactoMineR" ))`

## Usage

If using the example input data provided (`"Additional file 2.xlsx"`), simply open R or RStudio, ensure your working directory is set to the saved location of the `simmi` project folder, and run `source("scripts/main.R")`.

If using your own data, save the excel file in the `data/` folder, open the primary analysis script `"scripts/main.R"` and adjust the code as follows:

**1. Set up input data** -

-   Update the file path in the script: `file_path <- "data/Additional file 2.xlsx"`
-   Specify the following:
    -   `g` → sheet number of your chosen working variable
    -   `d` → column number of series displaying irregular pattern
    -   `k_start,` `k_end` → row numbers that define the irregular subperiod
    -   `f` → column number of chosen reference series for computing RV similarity

**2. Configure model parameters**

Adjust the following default model parameters as necessary:

-   `theta` → periodicity (e.g. 24 for hourly data)
-   `nfolds` → number of folds for time series block validation
-   `m` → number of imputations
-   `p` → number of lagged/leading features to be included
-   `step` → for computing naive forecast
-   `null` → null value for significance testing
-   `method` → "enet_ts_strict" or "lasso_ts"
-   `nw_lag` → number of time lags considered when computing Newey-West standard errors

**3. Adjust code for data preprocessing**

The temporary working data frame must incorporate relevant time covariates as factor variables in the end columns. Review lines 77-82 of the `"main.R"` script and either comment out or tweak code as necessary.

**4. Run the entire main script**

## Outputs

The script generates several outputs in the `output/` folder:

-   `results.csv`\
    → pooled MASE and p-values across shifts

-   `plot.png`\
    → visualization of pooled MASE vs shift

-   `my_spatial_ts_repaired.xlsx`\
    → repaired multivariate time-series data

-   `similarity.csv`\
    → RV similarity versus reference location before and after repair

## Notes

-   Ensure all custom scripts in the `R/` folder are available before running the main script. That is, ensure your working directory is set to the saved location of the `simmi` project folder
-   The algorithm may reduce `k_t` automatically if data is insufficient
-   If results are not satisfactory, increase `k_t` (e.g. to `10*theta`) and re-run
-   Output files will be overwritten each run unless renamed

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

This work utilises open-source R packages and builds upon established methodologies in:

-   Time-series forecasting
-   Penalised regression
-   Multiple imputation via chained equations
-   Statistical hypothesis testing
-   Multivariate similarity analysis
