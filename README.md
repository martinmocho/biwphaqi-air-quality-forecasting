# Baseline Impact-Weighted Philippine Air Quality Index (BIWPhAQI)

**Formulation, Characterization, and Time Series Forecasting**

A health-weighted composite air quality index for the Philippines, built from six pollutants and weighted by long-term respiratory mortality relative risks, then forecast with SARIMA, ETS, and neural network autoregressive models. Written in R.

> Final paper for **STT263A / STT163A**, Department of Mathematics and Statistics, De La Salle University, Manila (December 2025).
> The full write-up is in [`report/`](report/STT263A_BIWPhAQI_Report.pdf).

## Objectives

1. Construct a Philippine air quality index from pollutant concentration data, using transformed mean relative risks for long-term respiratory mortality as pollutant-specific weights.
2. Establish the **BIWPhAQI**, a population-weighted composite time series built from the index and barangay-level population data.
3. Develop and compare ARIMA, ETS, and NNAR forecasting models using a training–test split.
4. Produce a six-year forecast using the best-performing model.

## Approach

1. **Index construction.** NO₂, O₃, SO₂, and CO are converted to µg/m³ (assuming 27 °C and 1 atm). Each pollutant is weighted by β = ln(RR)/ΔC, where RR is the mean relative risk for long-term respiratory mortality from international studies such as HRAPIE-2. The pollutants are combined into a single covariance-adjusted index, so overlapping pollutant signals are not double counted.
2. **Population weighting.** Barangay-level values are averaged into a national daily series weighted by population, then averaged monthly (January 2003 onward). Population values for 2021 and 2022 were filled in with per-barangay linear trends.
3. **Characterization.** Summary statistics, outlier detection, classical and STL decomposition, seasonal plots, and ACF/PACF.
4. **Modeling.** Box–Cox transformed series, trained on data through December 2020 and tested on 2021–2022. Candidate SARIMA models were compared by AIC and residual diagnostics against an ETS model and a neural network autoregressive (NNAR) model.
5. **Forecasting.** Test-period evaluation and a six-year forecast through 2026.

## Key findings

- The selected models were **SARIMA(6,1,1)(1,0,1)[12]**, **ETS(A,N,A)**, and **NNAR(5,1,4)**, all with residuals resembling white noise.
- The neural network had the best short-term accuracy, while the SARIMA model performed best for longer-term projections.
- Forecasts show a generally stable level with a repeating seasonal pattern, suggesting stable air quality from a health-impact perspective.
- Observed fluctuations in the index lined up with real-world atmospheric events, indicating that it responds to actual conditions.

**Limitations noted in the paper:** relative risks come from non-local studies, the data end in 2021 (with later population values extrapolated), national aggregation can mask local pollution episodes, and the index scale has no direct public-health interpretation.

## Repository structure

```
.
├── R/
│   └── biwphaqi_forecasting.R        # annotated analysis code, from raw data to forecasts
├── report/
│   └── STT263A_BIWPhAQI_Report.pdf
├── data/
│   └── README.md                     # where to get the data (not redistributed here)
└── README.md
```

## Reproducing the analysis

1. Download the data as described in [`data/README.md`](data/README.md) and place the files in `data/`.
2. Install the R packages:

   ```r
   install.packages(c("fpp2", "devtools", "GGally", "patchwork", "seasonal", "moments",
                      "tseries", "lubridate", "dplyr", "tidyr", "fracdiff", "urca"))
   ```

3. Run the script from the data folder:

   ```r
   setwd("data")
   source("../R/biwphaqi_forecasting.R")
   ```

   Note that the ARIMA grid search in the script fits many models and can take a while.

## Authors

- Karl Uriel B. Dela Cruz
- Martin Johan M. Ocho
- Arjaye B. Ortiz

<!-- TODO (Martin): add a line describing your own contribution, e.g. "Martin: ..." -->

## Data source

Pollutant concentrations and population estimates come from the Climate and Air Quality data of [Project CCHAIN](https://doi.org/10.34740/KAGGLE/DS/4918229) (Thinking Machines Data Science, Inc., 2024), which draws on CAMS, MERRA-2, ERA5, and WorldPop.
