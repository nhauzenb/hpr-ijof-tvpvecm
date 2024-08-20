### Code package: N. Hauzenberger, M. Pfarrhofer, & L. Rossini (202x). Sparse time-varying parameter VECMs with an application to modeling electricity prices, *International Journal of Forecasting*, conditionally accepted.

[**Publication (open access).**](https://www.dropbox.com/scl/fi/j4yh5t0xs3j62dpjhre2r/HPR_TVP-VECM.pdf?rlkey=b88yuazdmgsgkta8rjggz611m&dl=0)

[**Final working paper version.**](https://www.dropbox.com/scl/fi/j4yh5t0xs3j62dpjhre2r/HPR_TVP-VECM.pdf?rlkey=b88yuazdmgsgkta8rjggz611m&dl=0)

### Electricity price data for the empirical application.
In the empirical application, we focus on modeling European electricity prices. Moreover, in a forecast exercise we consider hourly prices for Germany. 


For the first application, we use daily prices to estimate our model jointly for nine different regional markets: Baltics (BALT), Denmark (DK), Finland (FI), France (FR), Germany (DE), Italy (IT), Norway (NO), Sweden (SE) and Switzerland (CH). Prices for BALT and the Nordic countries (DK, FI, NO, and SE) are obtained from [Nord Pool](https://www.nordpoolgroup.com); the German, Swiss, and French hourly auction prices are from the power spot market of the [European Energy Exchange (EEX)](https://www.eex.com/en/); for the Italian prices, we use the single national prices (PUN) from the Italian system operator [Gestore dei Mercati Energetici (GEM)](https://www.mercatoelettrico.org/it/). We preprocess the data for daylight saving time changes to exclude the 25th hour in October and to interpolate the 24th hour in March. As additional exogenous factors, we consider daily prices for coal and fuel (both available through [ENTSO-E](https://www.entsoe.eu)) and interpolate missing values for weekends and holidays. In particular, we use the closing settlement prices for coal (LMCYSPT) and one month forward ICE UK natural gas prices (NATBGAS) due to their importance for the dynamic evolution of electricity prices and potential cointegration relationships. In our second application, the forecast comparison, we choose hourly day-ahead prices for Germany as our main country of interest, and focus on daylight hours and an average of the night hours. For forecasting electricity prices, we also include renewable energy sources (RES) in the form of forecasted demand, forecasted wind generation, and forecasted photovoltaic solar generation as exogenous factors (also available through [ENTSO-E](https://www.entsoe.eu)). We provide the data for both applications as a .rda file in the folder [`ELP data`](./ELP-data). While [`MC ELP daily`](./ELP-data/MC_elp_daily.rda) provides daily electricity prices for nine different regional markets, [`DE ELP hourly`](./ELP-data/DE_elp_hourly.rda) provides hourly data for Germany used in the forecasting exercise.

### Software information: 
All model estimations were carried out on a computing cluster using parallel processing. The cluster is equipped with 400 Intel E5-2650v3 2.3 GHz cores, managed through a Sun Grid Engine. We used R (version 4.1.3) as our main software. The libraries and packages used for estimation and analysis included `Rcpp`, `forecast`, `GIGrvg`, `MASS`, `Matrix`, `mvtnorm`, `stochvol`, `glasso`, `scoringRules`, and `MCS`. Additional packages used to create figures and tables were `data.table`, `dplyr`, `ggplot2`, `Hmisc`, `reshape2`, `tidyr`, `zoo`, and `scales`.

### Estimation files to produce a single forecast:
**!NOT DONE YET!** 

### Replication codes:
The replication codes reproduces all figures and tables in the manuscript. 

#### Replication of results shown in Sub-Section 4.2: *"Nonlinearities and cointegration in European electricity prices"*

In Sub-Section 4.2, we show and discuss:

* **Figure 1.** The posterior probability of the rank (PPR) over time for the most flexible specification.
* **Figure 2.** The posterior inclusion probability (PIP) over time for autoregressive coefficients (panel (a)) and covariance matrices (panel (b)) of the most flexible specification.
* **Figure 3.** The posterior median of all SV processes, with the red line denoting the first principal component (PCA) of all nine SV processes.
* **Table B.1.** Posterior estimates for the state equations of the error variances.

The file [`Subsection4-2_insample`](Subsection4-2_insample.R) can be used to replicate Figure 1, panels (a) and (b) of Figure 2, panels (a) and (b) of Figure 3, as well as Table B.1 in the Appendix. To do so, the file loads the output data from the folder [`INSAMPLE-data`](./INSAMPLE-data) and stores the figures and tables in a folder labeled `INSAMPLE-res`.


#### Replication of results shown in Sub-Section 4.3: *"Forecast results"*

In Sub-Section 4.3, we show and discuss:

* **Table 1.** Forecast performance for point and density forecasts
* **Table 2.** Model confidence set (MCS) for density forecasts at a 10 percent significance level using the *TR* test statistic.
* **Table B.2.** Model confidence set (MCS) for density forecasts at a 10 percent significance level using the *T-max* test statistic.
* **Table B.3.** Model confidence set (MCS) for density forecasts at a 25 percent significance level using the *TR* test statistic.
* **Table B.4.** Model confidence set (MCS) for density forecasts at a 25 percent significance level using the *T-max* test statistic.

The file [`Subsection4-3_OOS-main`](Subsection4-3_OOS-main.R) can be used to replicate the main Table — Table 1 — of the paper, while [`Subsection4-3_OOS-MCS`](Subsection4-3_OOS-MCS.R) replicates Tables 2, B.2, B.3, and B.4, which show results on the MCS. The .R-files access the output data in the folder [`OOS-data`](./OOS-data) and stores the tables in a folder labeled `OOS-res`.

To make the final .tex tables more accessible, they can be collected as a pdf through the file [`Tables-viewer.tex`](Tables-viewer.tex).
