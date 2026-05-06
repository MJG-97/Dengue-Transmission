# SIR-SI Dengue Transmission Model with Climate Forcing
# Modeling the transmission of dengue fever in Mayote

## Description

This repository contains a compartmental SIR-SI mathematical model (5 compartments) for dengue fever transmission dynamics, incorporating climatic forcing to capture the impact of temperature and precipitation on the transmission rate. Three climate forcing models are compared: constant, Lambrechts (2011), and Gaussian.

## Code Structure

- Data loading and preparation
- Climate functions (Lambrechts, Gaussian with precipitation)
- SIR-SI ODE model (5 equations)
- Weekly incidence prediction functions
- Parameter estimation via Maximum Likelihood Estimation (Negative Binomial)
- Optimization and model comparison (AIC)
- Visualizations (observed/predicted incidences, climate-cases correlations, β(t) dynamics)

## Compared Models

| Model | Description |
|-------|-------------|
| Constant | β(t) constant over time |
| Lambrechts | β(t) depends on temperature via a function from the literature |
| Gaussian | β(t) depends on temperature (Gaussian kernel) and precipitation, with a lag of l weeks |

## Estimation Methods

- **MLE**: Maximum Likelihood Estimation with Negative Binomial distribution (default)
- **WLS**: Weighted Least Squares (available)

## Data

The data used for this project are confidential and not included in this repository. The code expects a CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `n_cases` | Weekly case counts |
| `temp_mean_week` | Weekly mean temperature (°C) |
| `rain_mean_week` | Weekly mean precipitation (mm) |

## Prerequisites
Required R packages:
```r
dplyr, ggplot2, tidyr, deSolve, gridExtra, readr, lubridate, scales, knitr, numDeriv
```

