\name{GetMaxLogPP}
\alias{GetMaxLogPP}

\title{Maximum Log-Posterior}
\description{Given a partition, compute the corresponding maximum of the log-posterior distribution}

\usage{GetMaxLogPP(d, init.part, g.prob, hyper.par)}

\arguments{
    \item{d}{n-by-m logical data matrix}
    \item{init.part}{integer vector of length n or a filepath to the initial partition}
    \item{g.prob}{vector of length 3 with prior probabilities of "properties"}
    \item{hyper.par}{m-by-4-by-2 double array with Beta distribution hyperparameters}
}

\details{
init.part can be either an integer vector or the path to the text file containing the initial partition. The text file can be thought as a CSV file with n rows and 1 column.

Probabilities contained in g.prob represent: uninformative attribute, informative attribute but not characteristic for any cluster, informative attribute and characteristic for a cluster. They should sum to 1.

hyper.par[, , 1] should contain the Beta distribution alpha parameters while hyper.par[, , 2] should contain the Beta distribution beta parameters. The 4 columns of each matrix correspond to the 4 statuses. Use the function 'InitHyperParameters' to obtain an 'hyper.par' variable that should work for most applications.
}

\value{
A list containing the following arguments:
\item{partition}{vector of cluster indices corresponding to the solution. Clusters are sorted in descending order according to their size}
\item{max.logpp}{value of the log posterior probability}
\item{idx.max}{m-by-k.tot integer matrix representing not informative, informative and characteristic attributes, where m is the total number of columns and k.tot is the total number of clusters}
}
