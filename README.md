### Code package: N. Hauzenberger, M. Pfarrhofer, & L. Rossini (202x). Sparse time-varying parameter VECMs with an application to modeling electricity prices, *International Journal of Forecasting*, conditionally accepted.

[**Publication (open access).**](https://www.dropbox.com/scl/fi/j4yh5t0xs3j62dpjhre2r/HPR_TVP-VECM.pdf?rlkey=b88yuazdmgsgkta8rjggz611m&dl=0)

[**Final working paper version.**](https://www.dropbox.com/scl/fi/j4yh5t0xs3j62dpjhre2r/HPR_TVP-VECM.pdf?rlkey=b88yuazdmgsgkta8rjggz611m&dl=0)

### Electricity price data for the empirical application.
In the empirical application, we focus on modeling European electricity prices. Moreover, in a forecast exercise we consider hourly prices for Germany. 
For the first application, we use daily prices to estimate our model jointly for nine different regional markets: Baltics (BALT), Denmark (DK), Finland (FI), France (FR), Germany (DE), Italy (IT), Norway (NO), Sweden (SE) and Switzerland (CH). Prices for BALT and the Nordic countries (DK, FI, NO, and SE) are obtained from *Nord Pool*; the German, Swiss, and French hourly auction prices are from the power spot market of the *European Energy Exchange (EEX)*; for the Italian prices, we use the single national prices (PUN) from the Italian system operator *Gestore dei Mercati Energetici (GEM)*. As additional exogenous factors, we consider daily prices for coal and fuel. In particular, we use the closing settlement prices for coal (LMCYSPT) and one month forward ICE UK natural gas prices (NATBGAS) due to their importance for the dynamic evolution of electricity prices and potential cointegration relationships. In our second application, the forecast comparison, we choose hourly day-ahead prices for Germany as our main country of interest, and focus on daylight hours and an average of the night hours. For forecasting electricity prices, we also include renewable energy sources (labeled RES) in the form of forecasted average daily demand, forecasted average daily wind generation, and forecasted average daily photovoltaic solar generation as exogenous factors. We provide the data as a .rda file [`MacroUncQ`](./data/MacroUncQ.rda).

### Estimation files to produce a single forecast:
We use a hold-out period of approximately a year and a half (550 hold-out observations). We consider one-step-ahead predictions, which implies that we forecast each individual hour for the following day. We provide replication codes for the a single forecast for the second application. 



we use the macroeconomic uncertainty measure of [Jurado, Ludvigson, and Ng (2015, AER)](https://www.aeaweb.org/articles?id=10.1257/aer.20131193) provided (and regularly updated) on the [web page of Sydney C. Ludvigson](https://www.sydneyludvigson.com/macro-and-financial-uncertainty-indexes). For the different macroeconomic variables, we rely on the popular FRED-QD dataset provided by the *Federal Reserve Bank of St. Louis* and publicly available [here](https://research.stlouisfed.org/econ/mccracken/fred-databases/). Our quarterly sample spans from 1965Q1 to 2019Q4. Sub-section 4.2 in the paper and Table B.1 in the Online Appendix shows the set of variables included for different model sizes. 


**1.) [`A simple univariate example with a GP regression`](!univariate-GPreg.R):** In Sub-section 2.2, we illustrate the GP regression by means of a simple univariate example. We model US GDP growth as a function of the first lag of a macroeconomic uncertainty measure for a sub-sample around the global financial crisis. This stylized example, highlights the role of the kernel and its hyperparameters crucially impacts the posterior estimates of the conditional mean function. The corresponding estimation file replicates Figure 1 of the paper. 

**2.) Conjugate GP-VAR with SV:** Based on the single realization of the Basu & Bundick (2017, ECTA) DSGE model [`BB_realization`](./data/BB_realization.csv), the file [`GPVAR_eqbyeq`](!GPVAR_eqbyeq.R) allows to estimate of a conjugate Gaussian process vector autoregression (GP-VAR) on an equation-by-equation basis. The [`GPVAR_girfs`](!GPVAR_girfs.R) collects the equation-wise estimates and allows to compute generalized impulse response functions (GIRFs). In addition, the folder [`gpvar funcs`](./gpvar_funcs/) contains the MCMC sampler for a GP regression and a C++ function to create a squared exponential kernel:

* [`MCMC sampler for conjugate GP-VAR with SV`](./gpvar_funcs/gp_eqbyeq_mcmc.R) 
* [`Squared exponential kernel function`](./gpvar_funcs/sqexp_kernel.cpp)


Replication codes come without technical support of any kind.
