# Zika_IAR
Estimating the ZIKV infection attack rate in the Americas

Data folder contains csv file containing cumulative incidence data for each country included in the manuscript.

IAR_projections folder contains projected IAR for each PAHO country and territory (see text for details).

The Bayesian model used to estimate IAR for each country is provided in the Zika_IAR_estimates_binomial_betapriors_noMprior.stan file. This template file was modified for each country to exclude data types not available (or only available at the national level) for that particular country.
