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
# FindOptimalPartition
#
# Description:
#
# Main search algorithm
#
# Arguments:
#
#          d: v-by-m logical data matrix
#  init.part: integer vector or character string
#     g.prob: double vector of length 3
#  hyper.par: m-by-4-by-2 double array
#  max.clust: integer scalar
#     t.iter: integer scalar or NULL
#   bak.file: character string or NULL
#
# Value:
#
# A list containing the following arguments:
# partition: vector of cluster indices corresponding to the (local) optimum.
#            Clusters are sorted in descending order according to their size
# max.logpp: value of the log posterior probability
#   idx.max: m-by-k.tot integer matrix that represents not informative,
#            informative and characteristic attributes, where m is the total
#            number of columns while k.tot is the total number of clusters
#
###############################################################################
FindOptimalPartition <- function(d, init.part, g.prob, hyper.par, max.clust,
                                 t.iter, bak.file) {
    v <- nrow(d)
    m <- ncol(d)

    partition <- InitializePartition(init.part, v)

    k.tot <- max(partition)

    # To preallocate memory, we must choose an upper bound to the maximum
    # number of clusters. The total number of clusters can't be greater than v,
    # but v could be a very very large number. Choose max.clust as the upper
    # bound unless the current partition has more than max.clust clusters...
    max.clust <- max(max.clust, k.tot)

    # preallocate space
    n.ext <- tabulate(partition, nbins=max.clust)

    # fc = filled.clusters
    fc <- (n.ext > 0)

    n1s.ext <- matrix(0, nrow=m, ncol=max.clust)
    n1s.ext[, fc] <- ComputeCounts(d, partition)

    a <- hyper.par[, , 1]
    b <- hyper.par[, , 2]

    # - n and n1s have to be updated every time we move units between
    #   clusters
    # - g probabilities do not change during the execution of the algorithm
    # - hyperparameters (a and b) don't change either
    # - o probabilities, instead, depend on the total number of clusters. We
    #   have to recompute them every time we change k.tot
    log.gp <- LogGProbs(g.prob, m)
    log.op <- LogOProbs(m, k.tot)

    # compute the maximum log(p(M | D)) for the current (initial) partition
    logpp <- ComputeMaxLogPP(n.ext[fc], n1s.ext[, fc, drop=FALSE], a, b,
                             log.gp, log.op)

    cat("\nInitial status\nTotal number of clusters:", ncol(logpp$z),
        "\nMaximum log(p(M | D)):", logpp$max.logpp, "+ c\n")

    # save the right hyperparameters for this local optimum
    aa <- matrix(0, nrow=m, ncol=max.clust)
    bb <- matrix(0, nrow=m, ncol=max.clust)
    
    # If array A has dimension c(d_{1}, ..., d_{l}, ..., d_{L}), to access
    # element A[i_{1}, ..., i_{l}, ..., i_{L}] it is possible to use the
    # following linear index
    # linear.idx <- i_{1} + d_{1} * (i_{2} - 1) + ...
    #               + (d_{1} * ... * d_{l-1}) * (i_{l} - 1) + ...
    #               + (d_{1} * ... * d_{L-1}) * (i_{L} - 1)
    #
    # A[i_{1}, ..., i_{l}, ..., i_{L}] == A[linear.idx]
    idx.lin <- as.numeric((1:m) + m * (logpp$z - 1))
    aa[, fc] <- a[idx.lin]
    bb[, fc] <- b[idx.lin]

    # as the empty cluster, choose the first one with 0 units
    if (all(fc)) {
        new.cluster <- 0
        w.set <- sample(which(fc))
    } else {
        new.cluster <- which(!fc)[1]
        w.set <- sample(c(which(fc), new.cluster))
    }

    w <- w.set[w.idx <- 1]

    if (w == new.cluster) {
        keep.splitting <- TRUE
    } else {
        keep.splitting <- FALSE
    }

    # algorithm levels
    #
    # we start from level 1 and we move to the next one only when there are no
    # changes in the logpp. If there is a change in a level i > 1, we fall
    # back to level 1. We iterate the process until no changes are found.
    #
    # level 1: collect an unknown number of units from a single cluster
    # level 2: collect an unknown number of units from any cluster
    # level 3: switch units

    # -- TODO -- a complete rework is required, maybe a completely different
    #            approach
    al.levels <- 3
    al <- 1

    counter <- 1
    r <- NULL
    has.changed <- FALSE

    done <- FALSE
    while (!done) {
        if (al == 1) {
            # try to move units from cluster 'k' (to be found) to cluster 'w'
            if (!keep.splitting) {
                r <- OperatorCollectUnitsTypeA(w, d, partition, n.ext, n1s.ext,
                                               aa, bb, logpp$z, log.gp, log.op,
                                               0)
            } else {
                r <- OperatorCreateClusterTypeA(w, d, partition, n.ext,
                                                n1s.ext, a, b, aa, bb, logpp$z,
                                                log.gp, log.op, 0)
            }
        } else if (al == 2) {
            # try to move units from any cluster to cluster 'w'
            if (!keep.splitting) {
                r <- OperatorCollectUnitsTypeB(w, d, partition, n.ext, n1s.ext,
                                               aa, bb, logpp$z, log.gp, log.op,
                                               0)
            } else {
                r <- OperatorCreateClusterTypeB(w, d, partition, n.ext,
                                                n1s.ext, a, b, aa, bb, logpp$z,
                                                log.gp, log.op, 0)
            }
        } else if (al == 3) {
            # try to switch units
            r <- OperatorSwitchUnits(d, partition, n.ext, n1s.ext, aa, bb)
        }

        if (!is.null(r)) {
            # found a better solution!
            partition <- r

            n.ext <- tabulate(partition, nbins=max.clust)

            fc <- (n.ext > 0)

            n1s.ext[,  fc] <- ComputeCounts(d, partition)
            n1s.ext[, !fc] <- 0

            k.new <- sum(fc)

            if (k.new != k.tot) {
                k.tot <- k.new
                log.op <- LogOProbs(m, k.tot)
            }

            # compute the new log posterior maximum
            logpp <- ComputeMaxLogPP(n.ext[fc], n1s.ext[, fc, drop=FALSE], a,
                                     b, log.gp, log.op)

            idx.lin <- as.numeric((1:m) + m * (logpp$z - 1))
            aa[, fc] <- a[idx.lin]
            bb[, fc] <- b[idx.lin]

            # new.cluster position in w.set
            nc.pos <- w.set == new.cluster

            if (all(fc)) {
                # we reached the upper bound
                new.cluster <- 0
                keep.splitting <- FALSE
            } else {
                # find the new index for the empty cluster
                new.cluster <- which(!fc)[1]
            }

            # remove merged clusters
            w.set[!(w.set %in% unique(partition))] <- NA
            w.set[nc.pos] <- new.cluster

            has.changed <- TRUE
        } else if (keep.splitting) {
            keep.splitting <- FALSE
        }

        if (w.idx == length(w.set)) {
            fc <- (n.ext > 0)

            # as the empty cluster, choose the first one with 0 units
            if (all(fc)) {
                new.cluster <- 0
                w.set <- sample(which(fc))
            } else {
                new.cluster <- which(!fc)[1]
                w.set <- sample(c(which(fc), new.cluster))
            }

            w <- w.set[w.idx <- 1]

            if (w == new.cluster) {
                keep.splitting <- TRUE
            } else {
                keep.splitting <- FALSE
            }

            # if no changes have been found, move to next status
            # if a change has been found, go back to status 1
            if (has.changed) {
                al <- 1
                has.changed <- FALSE
            } else if (al < al.levels) {
                al <- al + 1
            } else {
                cat("\nDid not find a better solution in the last iteration.",
                    "Quitting.\n")

                done <- TRUE
            }
        } else if (!keep.splitting) {
            while (w.idx < length(w.set) &&
                   is.na(w <- w.set[w.idx <- w.idx + 1])) {}

            if (!is.na(w) && w == new.cluster) {
                keep.splitting <- TRUE
            } else if (is.na(w) && w.idx == length(w.set)) {
                # it could happen that the last element of w.set is NA. In this
                # case we have to reset w.set
                fc <- (n.ext > 0)

                # as the empty cluster, choose the first one with 0 units
                if (all(fc)) {
                    new.cluster <- 0
                    w.set <- sample(which(fc))
                } else {
                    new.cluster <- which(!fc)[1]
                    w.set <- sample(c(which(fc), new.cluster))
                }

                w <- w.set[w.idx <- 1]

                if (w == new.cluster) {
                    keep.splitting <- TRUE
                } else {
                    keep.splitting <- FALSE
                }

                # if no changes have been found, move to next status
                # if a change has been found, go back to status 1
                if (has.changed) {
                    al <- 1
                    has.changed <- FALSE
                } else if (al < al.levels) {
                    al <- al + 1
                } else {
                    cat("\nDid not find a better solution in the last",
                        "iteration. Quitting.\n")

                    done <- TRUE
                }
            }
        } else {
            w <- new.cluster
        }

        if (!is.null(t.iter) && (counter %% t.iter) == 0) {
            cat("\nStatus after", counter, "iterations.\nTotal number of",
                "clusters:", ncol(logpp$z), "\nMaximum log(p(M | D)):",
                logpp$max.logpp, "+ c\n")

            if (!is.null(bak.file)) {
                try(write.table(partition, row.names=FALSE, col.names=FALSE,
                                file=bak.file))
            }
        }

        counter <- counter + 1
    }

    # order the clusters with respect to decreasing size
    partition <- match(partition, sort(unique(partition)))
    k.set <- order(tabulate(partition), decreasing=TRUE)

    final.part <- numeric(v)
    final.idx.max <- matrix(0, nrow=m, ncol=ncol(logpp$z))

    for (i in 1:length(k.set)) {
        final.part[partition == k.set[i]] <- i
        final.idx.max[, i] <- logpp$z[, k.set[i]]
    }

    return(list(partition=final.part, max.logpp=logpp$max.logpp,
                idx.max=final.idx.max))
}
