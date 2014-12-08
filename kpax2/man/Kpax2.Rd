\name{Kpax2}
\alias{Kpax2}

\title{Bayesian Cluster Analysis of Categorical Data}
\description{Main function of K-Pax2 algorithm}

\usage{Kpax2(d, init.part, g.prob, hyper.par, max.clust=NULL, t.iter=50, bak.file=NULL)}

\arguments{
    \item{d}{n-by-m logical data matrix}
    \item{init.part}{integer vector of length n or a filepath to the initial partition}
    \item{g.prob}{vector of length 3 with prior probabilities of "properties"}
    \item{hyper.par}{m-by-4-by-2 double array with Beta distribution hyperparameters}
    \item{max.clust}{scalar integer. Maximum number of clusters allowed (used only to preallocate memory)}
    \item{t.iter}{integer scalar. Print a status message every 't.iter' iterations}
    \item{bak.file}{character string representing a path to a backup file. Set 'bak.file = NULL' to disable partition backup}
}

\details{
init.part can be either an integer vector or the path to the text file containing the initial partition. The text file can be thought as a CSV file with n rows and 1 column.

Probabilities contained in g.prob represent: uninformative attribute, informative attribute but not characteristic for any cluster, informative attribute and characteristic for a cluster. They should sum to 1.

hyper.par[, , 1] should contain the Beta distribution alpha parameters while hyper.par[, , 2] should contain the Beta distribution beta parameters. The 4 columns of each matrix correspond to the 4 properties. Use the function 'InitHyperParameters' to obtain an 'hyper.par' variable that should work for most applications.

The unknown number K of clusters is equal (at most) to the sample size n. If K << n (as it is usually the case), fixing max.clust=n would correspond to a waste of computer memory. If max.clust is not specified, a value of 'n' is used.

Set 't.iter <= 0' to disable status messages. Also, every t.iter iterations (if 't.iter > 0') the current partition is saved to bak.file.
}

\value{
A list containing the following arguments:
\item{partition}{vector of cluster indices corresponding to the solution. Clusters are sorted in descending order according to their size}
\item{max.logpp}{value of the log posterior probability}
\item{idx.max}{m-by-k.tot integer matrix representing not informative, informative and characteristic attributes, where m is the total number of columns and k.tot is the total number of clusters}
}
