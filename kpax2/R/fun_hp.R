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
# InitHyperParameters
#
# Description:
#
# Initizialize Beta distribution hyperparameters
#
# Arguments:
#
#   n: sample size
# n1s: vector of length m with the total number of ones for each column
#   r: value of the parameter when the feature is characteristic for some
#      cluster. The default value should work for most applications. As r goes
#      to infinity, a column will be flagged "characteristic" for a cluster if
#      there are only ones (or zeros) within clusters
#
# Value:
#
# m-by-4-by-2 double array
# x[, , 1]: m-by-4 double matrix of Beta distribution alpha parameters
# x[, , 2]: m-by-4 double matrix of Beta distribution beta parameters
#
###############################################################################
InitHyperParameters <- function(n, n1s, r=(log(0.001) / log(0.95))) {
    m <- length(n1s)

    hp <- array(0, dim=c(m, 4, 2))

    # uninformative attributes
    # If n > r, these two parameters sum to (r+1), i.e. the theoretical sample
    # size for the characteristic hyperparameters. This is done so because we
    # don't want them to be overwhelmed by the data.
    # The mean is the same to the one obtained with a Jeffreys prior
    if (n > r) {
        hp[, 1, 1] <- (r + 1) * (n1s + 0.5) / (n + 1)
        hp[, 1, 2] <- (r + 1) - hp[, 1, 1]
    } else {
        hp[, 1, 1] <- 0.5 + n1s
        hp[, 1, 2] <- 0.5 + n - n1s
    }

    # informative but not characteristic for any cluster
    hp[, 2, 1] <- 1
    hp[, 2, 2] <- 1

    # informative and characteristic... but for another cluster
    hp[, 3, 1] <- 1
    hp[, 3, 2] <- r

    # informative and characteristic for this cluster
    hp[, 4, 1] <- r
    hp[, 4, 2] <- 1

    return(hp)
}

###############################################################################
#
# LogGProbs
#
# Description:
#
# Augment the vector of prior probabilities g.prob to a m-by-4 matrix with
# log(g.prob) values. We sacrifice some memory to increase computational speed
#
# Arguments:
#
# g.prob: vector of length 3 with prior probabilities of "properties"
#             g.prob[1]: probability of an uninformative attribute
#             g.prob[2]: probability of an informative attribute, but not
#                        characteristic for any cluster
#             g.prob[3]: probability of an informative attribute,
#                        characteristic for just one cluster
#      m: total number of columns
#
# Value:
#
# m-by-4 double matrix with "properties" prior log-probabilities
#
###############################################################################
LogGProbs <- function(g.prob, m) {
    # attributes 3 and 4 are both from property 3 -> use g.prob[3] twice
    g <- c(g.prob, g.prob[3])

    # log(0) is -Inf
    # since we are going to multiply this matrix with a positive scalar
    #(1 / k.tot), no NaN can be produced even if some g_i are zero
    return(matrix(log(g), nrow=m, ncol=4, byrow=TRUE))
}

###############################################################################
#
# LogOProbs
#
# Description:
#
# O(mega) prior log-probabilities
#
# Arguments:
#
#     m: total number of columns
# k.tot: total number of clusters
#
# Value:
#
# m-by-4 double matrix with "attributes" prior log-probabilities
#
###############################################################################
LogOProbs <- function(m, k.tot) {
    o <- matrix(0, nrow=m, ncol=4)

    # informative and characteristic... but for another cluster
    o[, 3] <- log(k.tot - 1) - log(k.tot)

    # informative and characteristic for this cluster
    o[, 4] <- -log(k.tot)

    return(o)
}
