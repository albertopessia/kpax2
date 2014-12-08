\name{InitHyperParameters}
\alias{InitHyperParameters}

\title{Beta Distribution Hyperparameters}
\description{Initizialize Beta distribution hyperparameters}

\usage{InitHyperParameters(n, n1s, r = (log(0.001)/log(0.95)))}

\arguments{
    \item{n}{integer scalar for the sample size}
    \item{n1s}{vector of length m with the total number of ones for each column}
    \item{r}{value of the parameter when the feature is characteristic for some cluster}
}

\details{The default value of argument r should work for most applications. As r goes to infinity, a column will be flagged "characteristic" for a cluster if there are only ones (or zeros) within clusters}

\value{m-by-4-by-2 double array. x[, , 1]: m-by-4 double matrix of Beta distribution alpha parameters. x[, , 2]: m-by-4 double matrix of Beta distribution beta parameters}
