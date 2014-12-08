\name{SimulateDataset}
\alias{SimulateDataset}

\title{Categorical Data Simulation}
\description{Generate a random dataset to be analysed by K-Pax2}

\usage{SimulateDataset(nk, nc, n.attr, op=NA)}

\arguments{
    \item{nk}{integer vector of cluster sizes}
    \item{nc}{cardinality of the alphabet, i.e. number of categories (at least 2)}
    \item{n.attr}{integer vector of length 3}
    \item{op}{probability that a feature with property 3 will be characteristic for the considered cluster. If NA, default to 1 / length(nk)}
}

\details{
It is assumed that the length of 'nk' is equal to the total number of clusters while the sum of nk is equal to the sample size. Floating point numbers are rounded to integers.

'n.attr[j]' should represent the total number of columns classified as j, where j is 1 (noise), 2 (weak signal) and 3 (strong signal). sum(n.attr) is equal to M, the total number of columns.
}

\value{
A list containing the following arguments:
\item{D}{character n-by-M matrix}
\item{partition}{vector of cluster indices}
\item{col.class}{integer K-by-M matrix with column classification}
}
