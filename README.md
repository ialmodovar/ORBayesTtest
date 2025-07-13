# ORBayesTtest
Objective and Robust Bayesian $t$-test

## Authors:
Israel Almodóvar-Rivera

## Introduction

Student's \(t\)-test has been over 100 years since the discovery of one of the most fundamental statistical tests: Student's \(t\)-test. This shiny app employs the objective and robust Bayesian approach for hypothesis testing for one-sample and two-sample mean comparisons for the assumption of equal variances. Bayes factor and posterior probabilities are reported.

## Installation

ORBayesTtest requires

```
- R version 4.1.0 or higher.
- R packages: shiny, rmarkdown, reaxl, readODS, openxlsx
```

This app can be installed as an R package via the devtools package:

```R
library("devtools")
install_github("ialmodovar/ORBayesTtest")
```

To perform the analysis in R, you can use the ORBayesianTtest for inputting data vectors \( x\) and \(y\). You can also use ORBayesianTtest.sm() if you have the summary statistics,i.e., sample mean, sample standard deviation, and sample size.

To run the Shiny App

```R
library("ORBayesTtest")
ORBayesianTtest.app()
## You can also try
ORBayesTtest::ORBayesianTtest.app()
```

## How to use

This app enables users to upload a dataset with the variables of their choice or enter their values. The first row will be the names of the variables. The two-sample option allows the user to compare two columns, as well as compare the groups in a single column. I have included a tab for interpreting the Bayes Factors.


## Reference

Almodóvar-Rivera, I.A.; Pericchi-Guerra, L.R. An Objective and Robust Bayes Factor for the Hypothesis Test: One Sample and Two Population Means. Entropy 2024, 26, 88. [https://doi.org/10.3390/e26010088 ](https://www.mdpi.com/1099-4300/26/1/88).
