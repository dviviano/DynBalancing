% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{DynBalancing_History}
\alias{DynBalancing_History}
\title{Estimate dynamic average treatment effects for different lags effects}
\usage{
DynBalancing_History(
  panel_data,
  covariates_names,
  Time_name,
  unit_name,
  outcome_name,
  treatment_name,
  ds1,
  ds2,
  histories_length,
  fixed_effects = NA,
  pooled = F,
  cluster_SE = NA,
  params = list()
)
}
\arguments{
\item{panel_data}{data.frame format. It must be a panel with each row corresponding to a unit at a given date.}

\item{covariates_names}{vector of strings with column names of the covariates to include in the regression.}

\item{Time_name}{string denoting the column name of the panel corresponding to the date.}

\item{unit_name}{string denoting the column name of the panel corresponding to the individual.}

\item{outcome_name}{string denoting the column name of the panel corresponding to the outcome variable. Rows with missing outcomes will be removed.}

\item{treatment_name}{name of the column of the panel corresponding to the treatment assignment. Rows with missing treatment assignments will be removed.}

\item{ds1}{counterfactual treatment history of interest of the first potential outcome.}

\item{ds2}{counterfactual treatment history of interest of the second potential outcome.}

\item{histories_length}{numeric vector containing the number of lags for which the effect must be computed}

\item{pooled}{boolean indicating whether to run a pooled regression. If pooled = TRUE, the last T periods will not be used, with T equal to the length of the treatment ds1.}

\item{fixed_effect}{vector of strings containing the columns to include in the regression as fixed effects.}

\item{param}{See documentation.}
}
\value{
plots: plots of the main results

all_results:  matrix with ATE, variance, critical quantile, expected value of each potential outcome and corresponding critical quantile;
}
\description{
This function computes the ATE from a panel for a sequence of treatment histories.
Each row of the panel must contain a unit observed at a given date.
The function estimates the ATE using linear projections for potential outcomes
and balancing covariates dynamically. Reference: Viviano, Davide and Bradic Jelena. 2021. Dynamic covariate balancing: estimating treatment effect over time.
}
