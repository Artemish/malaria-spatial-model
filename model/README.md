# Climate-Malaria Spatiotemporal Model

## Model Overview

This project implements a spatiotemporal Bayesian hierarchical model to analyze the relationship between climate variables and malaria incidence in Mozambique districts. The model uses a Conditional Autoregressive (CAR) approach to account for both spatial and temporal dependencies in disease transmission patterns.

### Modeling Approach

**Data Structure:**
- **Spatial Units**: 130 districts in Mozambique, with focus on southern region analysis
- **Temporal Coverage**: 60 months of longitudinal data
- **Response Variable**: Monthly malaria case counts per district

**Climate Predictors:**
- **Precipitation (PRCP)**: District-level spatial totals with temporal lags (1-3 months)
- **Temperature**: Maximum daily land surface temperature with temporal lags (0-3 months)
- **Storm Events**: Binary indicator for extreme precipitation events

**Model Features:**

1. **Lagged Climate Effects**: Incorporates 1-4 month temporal lags to capture delayed climate impacts on malaria transmission

2. **Non-linear Relationships**: Uses natural splines with quartile-based knots to model non-linear climate-disease associations

3. **Spatial Autocorrelation**: Employs adjacency-based neighborhood matrices to account for spatial clustering

4. **Temporal Autocorrelation**: Implements AR(2) autoregressive structure to capture temporal dependencies

5. **Population Offset**: Uses log(population) as offset term to model incidence rates

**Statistical Framework:**
- **Family**: Poisson regression for count data
- **Bayesian Implementation**: CARBayesST package with MCMC sampling
- **Model Structure**: `ST.CARar()` with spatial CAR and temporal AR components
- **Sampling**: 500,000 iterations with 100,000 burn-in, thinning every 40 samples

**Key Model Formula:**
```
malaria ~ offset(log(population)) + 
          ns(lag1_PRCP, knots) + ns(lag2_PRCP, knots) +
          ns(TEMPmax, knots) + ns(lag1_TEMPmax, knots) + 
          ns(lag2_TEMPmax, knots) + ns(lag3_TEMPmax, knots)
```

**Output:**
The model produces parameter estimates, fitted values, and model fit statistics for understanding climate-malaria relationships across space and time, with particular focus on precipitation and temperature effects with appropriate lag structures.
