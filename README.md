# HWA-Dehradun
This is the repository for "Long-term evaluation of a city-scale sterilisation programme in Dehradun shows reduced free-roaming dog density"
Authors: Sarah Erskine et al.


# Overview

This project analyses ~10 years of longitudinal survey data from Dehradun, India, where a large-scale Animal Birth Control (ABC) programme has sterilised over 56,000 dogs since 2016.
The aim is to test whether increasing sterilisation coverage is associated with:

* Reduced dog population density
* Reduced reproduction
* Improved welfare indicators

All analyses are conducted using Bayesian hierarchical models in R.



# Repository contents 

This repository contains four primary model scripts:

Model 1 (Dens).R
→ Primary model: effect of sterilisation on dog population density
Model 2 (Change).R
→ Model of change in density relative to sterilisation
Model 3 (reproductive).R
→ Effects of sterilisation on reproductive indicators
(e.g. % pups, % lactating females)
Model 4 (Skin).R
→ Effects of sterilisation on welfare indicators
(e.g. % dogs with skin conditions)


# Methods

All models are implemented using the brms package (Stan backend), allowing:

* Hierarchical (multi-level) modelling across zones
* Estimation of uncertainty via posterior distributions
* Flexible handling of longitudinal data


  # Data

The data used in these models are derived from repeated street surveys in Dehradun:

* Biannual surveys (breeding and non-breeding seasons)
* Multiple zones across the city
* Recorded variables include:
  * Dog counts (standardised as Total Count per km)
  * % sterilised
  * % pups
  * % lactating females
  * % dogs with skin conditions

Note:
Population density is treated as a relative index rather than an absolute estimate.




# Reproducibility

To run the models:

1. Install required packages

install.packages(c("tidyverse", "brms", "tidybayes"))


2. Run scripts

Each script can be run independently


# Key findings (summary)

Across models, higher sterilisation exposure is consistently associated with:

* Lower dog population density
* Reduced reproductive output
* Improved welfare indicators



# Contact

Sarah Erskine
PhD Researcher, University of Cambridge
