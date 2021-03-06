% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Evaluation.R
\name{evaluation}
\alias{evaluation}
\title{evaluation()}
\usage{
evaluation(
  X,
  y,
  genes,
  fit.function = AIC,
  glm.family,
  HallOfFame.genes = NULL,
  HallOfFame.fit = NULL
)
}
\arguments{
\item{X}{m by n matrix of X as covariate}

\item{y}{list of numbers; the variable to be regressed on}

\item{genes}{n by P matrix; P columns of genes with length n (composed of 1s and 0s)}

\item{fit.function}{a user defined function; this function will give the
fitness value of any given gene and X, and y variable. By defaut, it is
the AIC.}
}
\description{
evaluation()
}
\author{
Jiayang Nie
}
