###############################################################################
#
# K-Pax2 - Bayesian Cluster Analysis of Categorical Data
#
# Copyright (c) 2014 Alberto Pessia <alberto.pessia@gmail.com>
#
# K-Pax2 is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# K-Pax2 is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# K-Pax2. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

###############################################################################
#
# InitializePartition
#
# Description:
#
# Load, check and sanitize the initial partition
#
# Arguments:
#
# init.part: numeric vector of length v representing an initial partition of
#            the data
#            OR
#            the path to the text file with the initial partition
#         v: sample size
#
# Value:
#
# Integer vector of length v, with cluster indices
#
###############################################################################
InitializePartition <- function(init.part, v) {
    if (is.character(init.part)) {
        if (length(init.part) > 1) {
            warning(paste("More than an initial partition has been specified.",
                          "Only the first one will be read."), immediate.=TRUE)

            init.part <- init.part[1]
        }

        if (length(init.part) != 1) {
            stop("Invalid initial partition input file.")
        }

        # keep only the first column
        init.part <- read.table(file=init.part)[[1]]
    }

    if (!is.numeric(init.part)) {
        init.part <- suppressWarnings(as.numeric(init.part))

        if (any(is.na(init.part))) {
            stop("Invalid initial partition.")
        }
    }

    # be sure that partition labels are integers
    partition <- as.integer(init.part + 0.5)

    if (length(partition) != v) {
        stop(paste("Partition length should be equal to the number of sample",
                   "units:", v))
    }

    # remove "gaps" and non-positive values from the partition
    #
    # Example:
    # [1, 1, 3, 1, -2, 4, 3] -> [2, 2, 3, 2, 1, 4, 3]
    partition <- match(partition, sort(unique(partition)))

    return(partition)
}

###############################################################################
#
# ComputeCounts
#
# Description:
#
# Count the total number of ones for each cluster and for each column
#
# Arguments
#
#         d: v-by-m logical data matrix
# partition: integer vector of cluster indices
#
# Value:
#
# m-by-k integer matrix with the counts of ones for each column and for each
# cluster
#
###############################################################################
ComputeCounts <- function(d, partition) {
    f <- function(k) {
        return(as.integer(colSums(d[partition == k, , drop=FALSE])))
    }

    clust <- sort(unique(partition))
    m <- ncol(d)

    return(vapply(clust, FUN=f, FUN.VALUE=integer(length=m), USE.NAMES=FALSE))
}

###############################################################################
#
# FindInitialPartition
#
# Description:
#
# Use pam (Partitioning Around Medoids) to initialize the partition
#
# Arguments
#
#          d: n-by-m logical data matrix
#     g.prob: vector of length 3 with prior probabilities of "properties"
#                 g.prob[1]: probability of an uninformative attribute
#                 g.prob[2]: probability of an informative attribute, but not
#                            characteristic for any cluster
#                 g.prob[3]: probability of an informative attribute,
#                            characteristic for just one cluster
#  hyper.par: m-by-4-by-2 double array. hyper.par[, , 1] should contain the
#             Beta distribution alpha parameters while hyper.par[, , 2] should
#             contain the Beta distribution beta parameters. The 4 columns of
#             each matrix correspond to the 4 properties. Use the function
#             'InitHyperParameters' to obtain an 'hyper.par' variable that
#             should work for most applications
#          D: distance matrix as an object of class "dist". If NULL
#             (the default), it will be set to the output of a call to "dist"
#             function with "binary" method
#      k.set: vector of integers representing the possible number of clusters.
#             If NULL (the default), it will be initialized based on a
#             dendrogram (built from D)
#    verbose: logical. Should status messages be printed?
#
# Value:
#
# Integer vector with a good initial partition
#
###############################################################################
FindInitialPartition <- function(d, g.prob, hyper.par, D=NULL, k.set=NULL,
                                 verbose=FALSE) {
  n <- nrow(d)

  if (is.null(D)) {
    if (verbose) {
      cat("Computing Hamming distance matrix... ")
    }

    # compute a basic distance matrix
    D <- dist(x=d, method="binary")

    if (verbose) {
      cat("done.\n")
    }
  }

  if (is.null(k.set)) {
    if (verbose) {
      cat("Trying to find a proper set of possible number of clusters... ")
    }

    H <- hclust(D)

    v <- max(c(ceiling(sqrt(n)), 100))

    logprobs <- rep(0, v)

    for (k in 1:v) {
      logprobs[k] <- GetMaxLogPP(d, cutree(H, k=k), g.prob,
                                 hyper.par)$max.logpp
    }

    k.n <- which.max(logprobs)

    k.set <- seq(k.n - 10, k.n + 10, 1)
    k.set <- k.set[k.set >= 1 & k.set <= n]

    if (verbose) {
      cat("done.\n")
    }
  }

  max.value <- -Inf
  init.partition <- rep(1, n)

  if (verbose) {
    cat("Searching the best initial partition... ")
  }

  for (i in 1:length(k.set)) {
    partition <- pam(x=D, k=k.set[i], cluster.only=TRUE, do.swap=FALSE)
    logprob <- GetMaxLogPP(d, partition, g.prob, hyper.par)$max.logpp

    if (logprob > max.value) {
      max.value <- logprob
      init.partition <- partition
    }
  }

  if (verbose) {
    cat("done.\n")
  }

  return(init.partition)
}
