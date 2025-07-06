# ORBayesTtest
Objective and Robust Bayesian \\(t\\)-test

## Authors:
Israel Almodóvar-Rivera

## Introduction

Shiny app for Bayesian hypothesis testing using intrinsic and the Berger Robust prior, and the effective sample size (TESS) Bayesian Information Criterion. Bayes factor and posterior probabilities are reported.

## Installation

ORBayesTtest requires

```
- R version 4.1.0 or higher.
- R packages: shiny, rmarkdown, reaxl, readODS, openxlsx
```

The user can download the .Rmd file and run it on their local machine. This app can be run 

```R
rmarkdown::run("OBayesTtest.Rmd")
```

## How to use

This app enables users to upload a dataset with the variables of their choice or enter their values. The first row will be the names of the variables. The two-sample option allows the user to compare two columns, as well as compare the groups in a single column. I have included a tab for interpreting the Bayes Factors.


## Reference

Almodóvar-Rivera, I.A.; Pericchi-Guerra, L.R. An Objective and Robust Bayes Factor for the Hypothesis Test: One Sample and Two Population Means. Entropy 2024, 26, 88. [https://doi.org/10.3390/e26010088 ](https://www.mdpi.com/1099-4300/26/1/88).
