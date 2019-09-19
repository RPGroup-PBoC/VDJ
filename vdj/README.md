# `vdj`
---
Welcome to the computational guts of the project. This module contains all
functions used in processing, analysis, and visualization of data in this work.
All functions are thoroughly document and can be viewed by typing
`vdj.module.function?` in a Python prompt. The following is a brief description
of each submodule.

* `__init__.py` \| Standard initiation file for the entire module.
* `bayes.py` \| A collection of functions for performing Bayesian inference via
  Stan using the [Stan probabilistic programming language](http://mc-stan.org)
  and the [PyStan](https://pystan.readthedocs.io/en/latest/) interface. 
* `io.py` \| A collection of functions for file input-output and loading of
  various constants, such as the sequence for each endogenous sequence.
* `stats.py` \| Functions for simple calculation of various summary statistics.
* `viz.py` \| Useful functions for data visualization and setting of the
  plotting theme. 