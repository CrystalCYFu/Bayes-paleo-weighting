# Bayes-paleo-weighting

This repository contains an implementation of the Bayesian framework for weighting climate model ensembles using palaeoclimate proxy data, applied to Arctic sea ice as reported in Fu et al. (in press). When using this framework, please cite:

Fu, C. Y., Osman, M. B., & Aquino-López, M. A. (in press). Palaeoclimate constraints narrow uncertainty in future Arctic winter sea ice projections. _Nature Geoscience_.

## Requirements

Python 3.9 and the following libraries:

-	numpy
-	pandas
-	matplotlib

The t-walk MCMC sampler, available [here](https://www.cimat.mx/~jac/twalk/).

Posterior files from the BaySIC proxy system model, available [here](https://github.com/CrystalCYFu/PyBaySIC).

## Data & Usage

All model and proxy data required to reproduce the primary of the study are included in this repository:

- `proxies_MH.csv`: mid-Holocene (MH) ln(P<sub>D</sub>IP<sub>25</sub>) and ln(P<sub>B</sub>IP<sub>25</sub>) proxy compilation
- `models_MH.xlsx`: corresponding predicted proxies, forward-modelled from PMIP MH SIC output
- `proxies_historical_1950-1999.csv`: historical proxies, forward-modelled from SIC observations (1950-1999 mean)
- `models_historical_1950-1999.xlsx`: corresponding predicted proxies, forward-modelled from CMIP historical SIC output

The model weighting scheme (`Bayes-paleo-weighting.py`) uses the above data and BaySIC posteriors (`server_dino_spavar_7900_2024-09-06_08-56-41.txt` and `server_bras_spavar_7900_2024-09-06_08-56-41.txt`) to derive posterior model weights, which are sampled using t-walk. The weights are then applied to historical simulations and future projections of Arctic winter sea ice:

- `MarchSIE_historical_future.xlsx`: historical and future (SSP1, 2 & 5) March SIE time series calculated from SIC output on each model's native grid, for SIE trajectories
- `MarchSIC_historical_future_climo.nc`: regridded March SIC model output as historical (1990-1999) and SSP1, 2 & 5 (2090-2099) climatologies, for ice cover probability maps
- `MarchSIC_MH_piControl_climo.nc`: as above but with MH and piControl climatologies instead, for reference

Details on all model output used, including the modelling centre, CMIP/PMIP phase, experiment ID, run ID, simulation years, and number of SIC files, are documented in `SIC_model_data.xlsx`. For more information on the data as well as the weighting scheme itself, please refer to the publication and the references therein.

## License

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License](http://creativecommons.org/licenses/by-nc/4.0/).

**Copyright (c) 2026**

## Contact

For questions, please email C. Y. Fu (cyf25@cam.ac.uk).
