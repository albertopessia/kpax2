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
# ComputeLogCondProb
#
# Description:
#
# Compute the log-probability of a unit, conditioned on its own cluster without
# the unit itself, that is log(p(y | c_{k, -y}))
#
# Arguments:
#
#         d: logical data matrix
# partition: vector of cluster indices
#         n: vector with the size of each cluster
#       n1s: m-by-k integer matrix with the number of ones for each column
#            and for each cluster
#         a: m-by-k double matrix of Beta distribution alpha parameters
#         b: m-by-k double matrix of Beta distribution beta parameters
#
# Value:
#
# double vector (length nrow(d)) with the log conditional probabilities
#
###############################################################################
ComputeLogCondProb <- function(d, partition, n, n1s, a, b) {
    # If nrow(d) and ncol(d) are both large, then the memory to be allocated
    # could be too much for the computer we are using. Let's cross our fingers

    # clusters' size after we remove the unit
    # length(n.vec) == nrow(d)
    n.vec <- n[partition] - 1

    # total number of ones after we remove the unit from its cluster
    # dim(n1s.mat) == c(nrow(d), ncol(d))
    n1s.mat <- t(n1s[, partition, drop=FALSE]) - d

    # Beta distribution hyper parameters
    # dim(aa) == dim(bb) == c(nrow(d), ncol(d))
    aa <- t(a[, partition, drop=FALSE])
    bb <- t(b[, partition, drop=FALSE])

    # dim(x0) == dim(x1) == dim(x2) == c(nrow(d), ncol(d))
    x0 <- log(aa + bb + n.vec)
    x1 <- log(aa + n1s.mat)
    x2 <- log(bb + n.vec - n1s.mat)

    # dim(log.p) == c(nrow(d), ncol(d))
    log.p <- ifelse(d, x1, x2) - x0

    return(apply(log.p, 1, sum))
}

###############################################################################
#
# ComputeLogCondProbW
#
# Description:
#
# Compute the log-probability of a unit conditioned on a cluster 'w', that is
# log(p(y | c_{w})). The difference from 'ComputeLogCondProb' is that cluster w
# is not the cluster the unit belongs to
#
# Arguments:
#
#     d: logical data matrix
#   n.w: integer scalar. cluster 'w' size
# n1s.w: integer vector (length m) with the number of ones, for each column,
#        observed in cluster w
#   a.w: double vector (length m) of Beta distribution alpha parameters for
#        cluster w
#   b.w: double vector (length m) of Beta distribution beta parameters for
#        cluster w
#
# Value:
#
# double vector (length nrow(d)) with the log conditional probabilities
#
###############################################################################
ComputeLogCondProbW <- function(d, n.w, n1s.w, a.w, b.w) {
    la <- log(a.w + n1s.w)
    lb <- log(b.w + n.w - n1s.w)
    lc <- log(a.w + b.w + n.w)

    v <- nrow(d)
    m <- ncol(d)

    result <- numeric(v * m)

    # let X be a v-by-m matrix
    # let idx.lin be a linear index, that is 1 <= idx.lin <= v*m
    #
    # given idx.lin, it is possible to obtain the row and column pair (i, j) as
    #
    # i <- ((idx.lin - 1) %% v) + 1
    # j <- floor((idx.lin - 1) / v) + 1
    idx.lin <- which(d)
    idx.col <- floor((idx.lin - 1) / v) + 1
    result[idx.lin] <- la[idx.col] - lc[idx.col]

    idx.lin <- which(!d)
    idx.col <- floor((idx.lin - 1) / v) + 1
    result[idx.lin] <- lb[idx.col] - lc[idx.col]

    return(apply(matrix(result, nrow=v, ncol=m), 1, sum))
}

###############################################################################
#
# ComputeLogCondProbZw
#
# Description:
#
# Compute the log-probability of a new cluster 'w' when a new candidate unit is
# added to the cluster. The vector Z_w is unknown and has to be found
#
# Arguments:
#
#      n: integer scalar. cluster 'w' size after the unit has been added
#    n1s: m-by-v integer matrix. v is the total number of candidate units while
#         m is the total number of columns of the original data matrix.
#         n1s[j, i] is the total number of ones observed at column 'j', when
#         the candidate unit 'i' is added to cluster 'w'
#      a: double vector (length m) of Beta distribution alpha parameters
#      b: double vector (length m) of Beta distribution beta parameters
# log.gp: m-by-4 double matrix of "properties" log probabilities
# log.op: m-by-4 double matrix of (new) "statuses" log probabilities
#
# Value:
#
# double vector (length v) with the log conditional probabilities
#
###############################################################################
ComputeLogCondProbZw <- function(n, n1s, a, b, z, log.gp, log.op) {
    sum.ab <- a + b
    const <- lgamma(sum.ab) - lgamma(a) - lgamma(b) +
             (log.gp / (ncol(z) + 1)) + log.op

    log.p <- matrix(-Inf, nrow=nrow(n1s), ncol=ncol(n1s))

    # don't change if z == 1 or z == 2
    idx <- apply(z == 1, 1, any)
    log.p[idx, ] <- const[idx, 1] +
                    lgamma(a[idx, 1] + n1s[idx, , drop=FALSE]) +
                    lgamma(b[idx, 1] + n - n1s[idx, , drop=FALSE]) -
                    lgamma(sum.ab[idx, 1] + n)

    idx <- apply(z == 2, 1, any)
    log.p[idx, ] <- const[idx, 2] +
                    lgamma(a[idx, 2] + n1s[idx, , drop=FALSE]) +
                    lgamma(b[idx, 2] + n - n1s[idx, , drop=FALSE]) -
                    lgamma(sum.ab[idx, 2] + n)

    # we need to choose between z = 3 or z = 4
    idx <- apply(z > 2, 1, any)

    h31 <- const[idx, 3] + lgamma(a[idx, 3] + n1s[idx, , drop=FALSE]) +
           lgamma(b[idx, 3] + n - n1s[idx, , drop=FALSE]) -
           lgamma(sum.ab[idx, 3] + n)

    h32 <- const[idx, 4] + lgamma(a[idx, 4] + n1s[idx, , drop=FALSE]) +
           lgamma(b[idx, 4] + n - n1s[idx, , drop=FALSE]) -
           lgamma(sum.ab[idx, 4] + n)

    log.p[idx, ] <- ifelse(h32 > h31, h32, h31)

    return(apply(log.p, 2, sum))
}

###############################################################################
#
# ComputeLogJointProb
#
# Description:
#
# Compute the log-probability for each cluster
#
# Arguments:
#
#   n: vector (length k.tot) with the size of each cluster
# n1s: m-by-k.tot integer matrix with the number of ones for each column and
#      for each cluster
#   a: m-by-k.tot double matrix of Beta distribution alpha parameters
#   b: m-by-k.tot double matrix of Beta distribution beta parameters
#
# Value:
#
# double vector (length k.tot) with the log probabilities of the clusters
#
###############################################################################
ComputeLogJointProb <- function(n, n1s, a, b) {
    sum.ab <- a + b

    log.probs <- lgamma(sum.ab) - lgamma(a) - lgamma(b) + lgamma(a + n1s) +
                 lgamma(t(t(b - n1s) + n)) - lgamma(t(t(sum.ab) + n))

    return(apply(log.probs, 2, sum))
}
