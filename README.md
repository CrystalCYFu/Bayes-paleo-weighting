# Bayes-paleo-weighting

This repository contains an implementation of the Bayesian framework for weighting climate model ensembles using palaeoclimate proxy data, applied to Arctic sea ice as reported in Fu et al. (in press).

All model and proxy data required to reproduce the primary of the study are included here. For details on the data and their sources, as well as the weighting scheme itself, please refer to the publication and the references therein.

We welcome further application of this framework; when using it, please cite:

Fu, C. Y., Osman, M. B., & Aquino-López, M. A. (in press). Palaeoclimate constraints narrow uncertainty in future Arctic winter sea ice projections. _Nature Geoscience_.

## Requirements

Python 3.9 and the following libraries:

-	numpy
-	pandas
-	matplotlib

The t-walk MCMC sampler, available [here](https://www.cimat.mx/~jac/twalk/).

Posterior files from BaySIC, available [here](https://github.com/CrystalCYFu/PyBaySIC).

## License

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License](http://creativecommons.org/licenses/by-nc/4.0/).

**Copyright (c) 2026**

## Contact

For questions, please email C. Y. Fu (cyf25@cam.ac.uk).
