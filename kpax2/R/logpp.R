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
# ComputeMaxLogPP
#
# Description:
#
# Compute the maximum log posterior probability log(p(M | D)) for a particular
# partition
#
# Arguments:
#
#      n: integer vector (length k.tot) with the size of each cluster
#    n1s: m-by-k.tot integer matrix with the number of ones for each column and
#         for each cluster
#      a: m-by-4 double matrix of Beta distribution alpha parameters
#      b: m-by-4 double matrix of Beta distribution beta parameters
# log.gp: m-by-4 double matrix of "properties" log probabilities
# log.op: m-by-4 double matrix of "statuses" log probabilities
#
# Value:
#
# A list containing the following arguments:
# max.logpp: maximum log posterior probability for the given arguments
#         z: m-by-k.tot integer matrix with the properties of each column, for
#            each cluster
#
###############################################################################
ComputeMaxLogPP <- function(n, n1s, a, b, log.gp, log.op) {
    m <- nrow(n1s)
    k.tot <- ncol(n1s)

    sum.ab <- a + b
    const <- lgamma(sum.ab) - lgamma(a) - lgamma(b) + (log.gp / k.tot) + log.op

    f <- function(k) {
        return(const + lgamma(a + n1s[, k]) + lgamma(b + n[k] - n1s[, k]) -
               lgamma(sum.ab + n[k]))
    }

    log.tot <- vapply(1:k.tot, FUN=f, FUN.VALUE=matrix(0, nrow=m, ncol=4),
                      USE.NAMES=FALSE)

    # uninformative sites
    h1 <- apply(log.tot[, 1, , drop=FALSE], 1, sum)

    # informative sites
    h2 <- apply(log.tot[, 2, , drop=FALSE], 1, sum)

    # informative and characteristic
    h31 <- log.tot[, 3, , drop=FALSE]
    h32 <- log.tot[, 4, , drop=FALSE]

    test.h3 <- (h32 > h31)

    h3 <- apply(ifelse(test.h3, h32, h31), 1, sum)

    idx.max <- ifelse(h3 > h2, ifelse(h3 > h1, 3, 1), ifelse(h2 > h1, 2, 1))

    max.logpp <- sum(h1[idx.max == 1], h2[idx.max == 2], h3[idx.max == 3])

    z <- matrix(idx.max, nrow=m, ncol=k.tot)

    tmp <- (idx.max == 3)
    if (any(tmp)) {
        z[z == 3 & test.h3[, 1, ]] <- 4
    }

    return(list(max.logpp=max.logpp, z=z))
}

###############################################################################
#
# GetMaxLogPP
#
# Description:
#
# Compute the maximum log posterior probability log(p(M | D)) for a given
# partition
#
# Arguments:
#
#         d: n-by-m logical data matrix
#      part: integer vector of length n representing an initial partition of
#            the data
#            OR
#            the path to the text file with the initial partition
#    g.prob: vector of length 3 with prior probabilities of "properties"
#                g.prob[1]: probability of an uninformative attribute
#                g.prob[2]: probability of an informative attribute, but not
#                           characteristic for any cluster
#                g.prob[3]: probability of an informative attribute,
#                           characteristic for just one cluster
# hyper.par: m-by-4-by-2 double array. hyper.par[, , 1] should contain the
#            Beta distribution alpha parameters while hyper.par[, , 2] should
#            contain the Beta distribution beta parameters. The 4 columns of
#            each matrix correspond to the 4 properties. Use the function
#            'InitHyperParameters' to obtain an 'hyper.par' variable that
#            should work for most applications
#
# Value:
#
# A list containing the following arguments:
# partition: vector of cluster indices corresponding to the (local) optimum
# max.logpp: value of the log posterior probability log(p(M | D))
#   idx.max: m-by-k.tot integer matrix that represents not informative,
#            informative and characteristic attributes, where m is the total
#            number of columns while k.tot is the total number of clusters
#
###############################################################################
GetMaxLogPP <- function(d, init.part, g.prob, hyper.par) {
    if (!is.matrix(d) || !identical(mode(d), "logical")) {
        stop("Expected a logical (binary) matrix, but got something else")
    }

    v <- nrow(d)
    m <- ncol(d)

    if (is.numeric(g.prob)) {
        if (length(g.prob) != 3) {
            stop("Variable 'g.prob' does not have a length of 3")
        }

        if (any(g.prob < 0)) {
            stop("Variable 'g.prob' contains negative probabilities")
        }

        # probabilities must sum to one
        g.prob <- g.prob / sum(g.prob)
    } else {
        stop("Variable 'g.prob' is not of numeric type")
    }

    if (is.numeric(hyper.par)) {
        dim.d <- dim(hyper.par)

        if (length(dim.d) != 3 || !all(dim.d == c(m, 4, 2))) {
            stop("Variable 'hyper.par' does not have the right dimensions")
        }
    } else {
        stop("Variable 'hyper.par' is not of numeric type")
    }

    partition <- InitializePartition(init.part, v)

    n <- tabulate(partition)
    n1s <- ComputeCounts(d, partition)

    k.tot <- max(partition)

    log.gp <- LogGProbs(g.prob, m)
    log.op <- LogOProbs(m, k.tot)

    x <- ComputeMaxLogPP(n=n, n1s=n1s, a=hyper.par[, , 1],
                         b=hyper.par[, , 2], log.gp=log.gp, log.op=log.op)

    return(list(partition=partition, max.logpp=x$max.logpp, idx.max=x$z))
}
