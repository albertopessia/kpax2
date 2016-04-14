\name{FindInitialPartition}
\alias{FindInitialPartition}

\title{Search a good starting point}
\description{Use 'pam' (Partitioning Around Medoids) to initialize the partition}

\usage{FindInitialPartition(d, g.prob, hyper.par, D=NULL, k.set=NULL, verbose=FALSE)}

\arguments{
    \item{d}{n-by-m logical data matrix}
    \item{g.prob}{vector of length 3 with prior probabilities of "properties"}
    \item{hyper.par}{m-by-4-by-2 double array with Beta distribution hyperparameters}
    \item{D}{pre-computed distance matrix}
    \item{k.set}{integer vector representing the candidate total number of clusters}
    \item{verbose}{logical. Should status messages be printed?}
}

\details{
If D is NULL (the default), it will be set to the output of a call to 'dist' with 'binary' as method.

If k.set is NULL (the default), it will be initialized by scanning cuts of a dendrogram (built from D).
}

\value{
Integer vector with a good starting point for the optimum search.
}
