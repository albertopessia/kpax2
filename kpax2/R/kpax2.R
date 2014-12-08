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
# Kpax2
#
# Description:
#
# Main function of K-Pax2 algorithm
#
# Arguments:
#
#          d: n-by-m logical data matrix
#  init.part: integer vector of length n representing an initial partition of
#             the data
#             OR
#             the path to the text file with the initial partition
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
#  max.clust: scalar integer. Maximum number of clusters allowed.
#             0 < max.clust <= n
#             Used only to preallocate memory. The maximum number of clusters K
#             is equal (at most) to the sample size n. If K << n (as it is
#             usually the case), fixing max.clust = n would correspond to a
#             waste of computer memory. If max.clust is not specified, a value
#             of 'n' is used
#     t.iter: integer scalar. Print a status message every 't.iter' iterations.
#             Set 't.iter <= 0' to disable status messages
#   bak.file: character string representing a path to a backup file. Every
#             t.iter iterations (only if 't.iter > 0'), the current
#             partition is saved to bak.file. Set 'bak.file = NULL' to disable
#             partition backup
#
# Value:
#
# A list containing the following arguments:
# partition: vector of cluster indices corresponding to the (local) optimum.
#            Clusters are sorted in descending order according to their size
# max.logpp: value of the log posterior probability log(p(M | D))
#   idx.max: m-by-k.tot integer matrix that represents not informative,
#            informative and characteristic attributes, where m is the total
#            number of columns while k.tot is the total number of clusters
#
###############################################################################
Kpax2 <- function(d, init.part, g.prob, hyper.par, max.clust=NULL, t.iter=50,
                  bak.file=NULL) {
    if (!is.matrix(d) || !identical(mode(d), "logical")) {
        stop("Expected a logical (binary) matrix, but got something else")
    }

    v <- nrow(d)
    m <- ncol(d)

    if (is.null(v) || v < 2) {
        stop("Sample size lesser than 2. Trivial solution")
    }

    if (is.numeric(g.prob)) {
        if (length(g.prob) != 3) {
            stop("Variable 'g.prob' does not have length 3")
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

    if (is.numeric(max.clust)) {
        if (length(max.clust) != 1) {
            stop("Variable 'max.clust' is not a scalar")
        }

        max.clust <- as.integer(max.clust + 0.5)

        if (max.clust > v) {
            max.clust <- v
        } else if (max.clust <= 0) {
            stop("Variable 'max.clust' is lesser than 1")
        }
    } else if (is.null(max.clust)) {
        max.clust <- as.integer(v)
    } else {
        stop("Variable 'max.clust' is not of numeric type")
    }

    if (is.numeric(t.iter)) {
        if (length(t.iter) != 1) {
            stop("Variable 't.iter' is not a scalar")
        }

        t.iter <- as.integer(t.iter + 0.5)

        if (t.iter <= 0) {
            t.iter <- NULL
        }
    } else {
        stop("Variable 't.iter' is not of numeric type")
    }

    if (is.character(bak.file)) {
        if (length(bak.file) != 1) {
            stop("Variable 'bak.file' is not a string")
        }

        finfo <- file.info(bak.file)

        if (!is.na(finfo$isdir) && finfo$isdir) {
            stop(paste("'", bak.file, "' is a directory", sep=""))
        }
    } else if (!is.null(bak.file)) {
        stop("Variable 'bak.file' is not a string")
    }

    result <- FindOptimalPartition(d, init.part, g.prob, hyper.par, max.clust,
                                   t.iter, bak.file)

    return(result)
}
