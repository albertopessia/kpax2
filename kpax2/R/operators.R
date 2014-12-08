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
# OperatorCollectUnitsTypeA
#
# Description:
#
# Try to move units from cluster 'k' (to be found) to cluster 'w'
#
# Arguments:
#
#         w: integer scalar, representing the cluster receiving the units
#         d: logical data matrix
# partition: integer vector (length n) of cluster indices
#         n: integer vector (length max.clust) with the size of each cluster
#       n1s: m-by-max.clust integer matrix with the number of ones for each
#            column and for each cluster
#        aa: m-by-max.clust matrix of Beta distribution alpha parameters
#        bb: m-by-max.clust matrix of Beta distribution beta parameters
#         z: m-by-k.tot integer matrix with "properties" and "statuses"
#    log.gp: m-by-4 double matrix of "properties" log probabilities
#    log.op: m-by-4 double matrix of "statuses" log probabilities
#         T: integer scalar, representing the fixed number of units to move.
#            Set 'T=0' to find the units in a "dynamic" approach
#
# Value:
#
# Integer vector with the new partition or NULL if no better solution has been
# found
#
###############################################################################
OperatorCollectUnitsTypeA <- function(w, d, partition, n, n1s, aa, bb, z,
                                      log.gp, log.op, T) {
    v <- nrow(d)
    m <- ncol(d)
    k.tot <- ncol(z)

    # to avoid floating point precision issues, use a tolerance of 1e-10 to
    # test for equality
    eps <- 10

    if (k.tot == 1) {
        # there is only one cluster! we cannot collect units anywhere!!
        return(NULL)
    }

    clusters <- which(n > 0)

    # this constant is used only if a cluster is merged with 'w'
    # k.tot -> k.tot - 1
    # We compute it here because it is the same for every cluster
    #
    # psi_{j,h,s} = omega_{j,h,s}^{(new)} / omega_{j,h,s}^{(old)}
    log.psi <- LogOProbs(m, k.tot - 1) - log.op

    f <- function(k) {
        # x[1] is the increment/decrement in the log posterior
        # x[-1] represents the units to be moved
        x <- rep(0, v+1)

        idx.c <- which(clusters == k)
        idx.k <- partition == k
        idx.k.moved <- rep(FALSE, v)

        # we save n[k] and n1s[, k] because they will change their values
        # during the loop
        n.k <- n[k]
        n1s.k <- n1s[, k, drop=FALSE]
        n.w <- n[w]
        n1s.w <- n1s[, w, drop=FALSE]

        log.phi <- 0
        log.prob.tot <- 0

        increment.best <- 0
        increment.old <- 0
        ac.old <- -Inf

        counter <- 1

        done <- FALSE
        while (!done) {
            log.p.k <- ComputeLogCondProb(d[idx.k, , drop=FALSE], rep(1, n.k),
                                          n.k, n1s.k, aa[, k, drop=FALSE],
                                          bb[, k, drop=FALSE])

            log.p.w <- ComputeLogCondProbW(d[idx.k, , drop=FALSE],
                                           n.w, n1s.w, aa[, w, drop=FALSE],
                                           bb[, w, drop=FALSE])

            log.probs <- log.p.w - log.p.k

            i.max <- which.max(log.probs)

            if (n.k == 1) {
                # moving "i" from cluster "k" to cluster "w" will decrease the
                # total number of clusters! We have to compute the bonus (or
                # penalty) factor for moving from k.tot to k.tot - 1
                idx.lin.not.merged <- as.numeric((1:m) + m *
                                                 (z[, -idx.c, drop=FALSE] - 1))

                t1 <- log.gp[idx.lin.not.merged] / (k.tot * (k.tot - 1)) +
                      log.psi[idx.lin.not.merged]

                idx.lin.merged <- as.numeric((1:m) + m *
                                             (z[, idx.c, drop=FALSE] - 1))

                t2 <- (log.gp[idx.lin.merged] / k.tot) + log.op[idx.lin.merged]

                log.phi <- sum(t1) - sum(t2)
            }

            log.prob.tot <- log.prob.tot + log.probs[i.max]
            increment <- log.phi + log.prob.tot

            # update info
            i <- which(idx.k)[i.max]

            idx.k[i] <- FALSE
            idx.k.moved[i] <- TRUE

            if (round(increment - increment.best, eps) > 0) {
                # we got a better solution!
                x[1] <- increment.best <- increment
                x[-1] <- idx.k.moved
            }

            counter <- counter + 1

            if (any(idx.k)) {
                n.k <- n.k - 1
                n.w <- n.w + 1

                n1s.k <- n1s.k - t(d[i, , drop=FALSE])
                n1s.w <- n1s.w + t(d[i, , drop=FALSE])

                if (T == 0) {
                    # if the increment is already positive (better solution
                    # than the original), keep going no matter what
                    #
                    # stop the search if the "curve is not convex"
                    # if ac.new > 0 then we are going "up", which is ok
                    # if ac.new < 0 then we are going "down". In this case we
                    #                    have to check if ac.new > ac.old
                    # if ac.new = 0 then we want to go on and see what happens
                    #                    next
                    ac.new <- round(increment - increment.old, eps)

                    if (increment < 0 && ac.new < 0 && ac.new < ac.old) {
                        done <- TRUE
                    } else {
                        increment.old <- increment
                        ac.old <- ac.new
                    }
                } else if (counter > T) {
                    done <- TRUE
                }
            } else {
                done <- TRUE
            }
        }

        return(x)
    }

    r <- vapply(clusters[clusters != w], FUN=f, FUN.VALUE=rep(0, v+1),
                USE.NAMES=FALSE)

    best.solution <- which.max(r[1, ])

    if (r[1, best.solution] > 0) {
        result <- partition
        result[r[-1, best.solution] > 0] <- w
    } else {
        result <- NULL
    }

    return(result)
}

###############################################################################
#
# OperatorCreateClusterTypeA
#
# Description:
#
# Try to create a new cluster collecting units from cluster 'k' (to be found).
# This is also equivalent to splitting cluster 'k'
#
# Arguments:
#
#         w: integer scalar, representing the (new) cluster receiving the units
#         d: logical data matrix
# partition: integer vector (length n) of cluster indices
#         n: integer vector (length max.clust) with the size of each cluster
#       n1s: m-by-max.clust integer matrix with the number of ones for each
#            column and for each cluster
#         a: m-by-4 double matrix of Beta distribution alpha parameters
#         b: m-by-4 double matrix of Beta distribution beta parameters
#        aa: m-by-max.clust matrix of Beta distribution alpha parameters
#        bb: m-by-max.clust matrix of Beta distribution beta parameters
#         z: m-by-k.tot integer matrix with "properties" and "statuses"
#    log.gp: m-by-4 double matrix of "properties" log probabilities
#    log.op: m-by-4 double matrix of "statuses" log probabilities
#         T: integer scalar, representing the fixed number of units to move.
#            Set 'T=0' to find the units in a "dynamic" approach
#
# Value:
#
# Integer vector with the new partition or NULL if no better solution has been
# found
#
###############################################################################
OperatorCreateClusterTypeA <- function(w, d, partition, n, n1s, a, b, aa, bb,
                                       z, log.gp, log.op, T) {
    # main differences from OperatorCollectUnitsTypeA are:
    # - we cannot merge any cluster with w, otherwise we would end up with an
    #   identical partition
    # - we don't know the instance of z for cluster 'w', so we need to compute
    #   it every time we move a unit
    # - k.new is always equal to (k.tot + 1)
    v <- nrow(d)
    m <- ncol(d)
    k.tot <- ncol(z)

    # to avoid floating point precision issues, use a tolerance of 1e-10 to
    # test for equality
    eps <- 10

    # remove clusters of size 1 from computation
    clusters <- which(n > 1)

    if (length(clusters) == 0) {
        return(NULL)
    }

    log.op.new <- LogOProbs(m, k.tot + 1)

    # this constant is the same for each cluster
    idx.lin <- as.numeric((1:m) + m * (z - 1))
    const <- sum(log.op.new[idx.lin] - log.op[idx.lin] -
                 (log.gp[idx.lin] / (k.tot * (k.tot + 1))))

    f <- function(k) {
        # x[1] is the increment/decrement in the log posterior
        # x[-1] represents the units to be moved
        x <- rep(0, v+1)

        idx.k <- partition == k
        idx.k.moved <- rep(FALSE, v)

        # we save n[k] and n1s[, k] because they will change their values
        # during the loop
        n.k <- n[k]
        n1s.k <- n1s[, k, drop=FALSE]
        n.w <- n[w]
        n1s.w <- n1s[, w]

        # keep track of the sum of log conditional probabilities
        lcp <- 0

        increment.best <- 0
        increment.old <- 0
        ac.old <- -Inf

        counter <- 1

        done <- FALSE
        while (!done) {
            log.p.k <- lcp + ComputeLogCondProb(d[idx.k, , drop=FALSE],
                                                rep(1, n.k), n.k, n1s.k,
                                                aa[, k, drop=FALSE],
                                                bb[, k, drop=FALSE])

            log.p.w <- ComputeLogCondProbZw(n.w + 1,
                                            t(d[idx.k, , drop=FALSE]) + n1s.w,
                                            a, b, z, log.gp, log.op.new)

            increment <- const + log.p.w - log.p.k

            # update info
            i.max <- which.max(increment)
            i <- which(idx.k)[i.max]

            lcp <- log.p.k[i.max]

            idx.k[i] <- FALSE
            idx.k.moved[i] <- TRUE

            if (round(increment[i.max] - increment.best, eps) > 0) {
                # we got a better solution!
                x[1] <- increment.best <- increment[i.max]
                x[-1] <- idx.k.moved
            }

            counter <- counter + 1

            if (sum(idx.k) > 1) {
                n.k <- n.k - 1
                n.w <- n.w + 1

                n1s.k <- n1s.k - t(d[i, , drop=FALSE])
                n1s.w <- n1s.w + d[i, ]

                if (T == 0) {
                    # if the increment is already positive (better solution
                    # than the original), keep going no matter what
                    #
                    # stop the search if the "curve is not convex"
                    # if ac.new > 0 then we are going "up", which is ok
                    # if ac.new < 0 then we are going "down". In this case we
                    #                    have to check if ac.new > ac.old
                    # if ac.new = 0 then we want to go on and see what happens
                    #                    next
                    ac.new <- round(increment[i.max] - increment.old, eps)

                    if (increment[i.max] < 0 && ac.new < 0 && ac.new<ac.old) {
                        done <- TRUE
                    } else {
                        increment.old <- increment[i.max]
                        ac.old <- ac.new
                    }
                } else if (counter > T) {
                    done <- TRUE
                }
            } else {
                done <- TRUE
            }
        }

        return(x)
    }

    r <- vapply(clusters, FUN=f, FUN.VALUE=rep(0, v+1), USE.NAMES=FALSE)

    best.solution <- which.max(r[1, ])

    if (r[1, best.solution] > 0) {
        result <- partition
        result[r[-1, best.solution] > 0] <- w
    } else {
        result <- NULL
    }

    return(result)
}

###############################################################################
#
# OperatorCollectUnitsTypeB
#
# Description:
#
# Try to move units from any cluster to cluster 'w'
#
# Arguments:
#
#         w: integer scalar, representing the cluster receiving the units
#         d: logical data matrix
# partition: integer vector (length n) of cluster indices
#         n: integer vector (length max.clust) with the size of each cluster
#       n1s: m-by-max.clust integer matrix with the number of ones for each
#            column and for each cluster
#        aa: m-by-max.clust matrix of Beta distribution alpha parameters
#        bb: m-by-max.clust matrix of Beta distribution beta parameters
#         z: m-by-k.tot integer matrix with "properties" and "statuses"
#    log.gp: m-by-4 double matrix of "properties" log probabilities
#    log.op: m-by-4 double matrix of "statuses" log probabilities
#         T: integer scalar, representing the fixed number of units to move.
#            Set 'T=0' to find the units in a "dynamic" approach
#
# Value:
#
# Integer vector with the new partition or NULL if no better solution has been
# found
#
###############################################################################
OperatorCollectUnitsTypeB <- function(w, d, partition, n, n1s, aa, bb, z,
                                      log.gp, log.op, T) {
    v <- nrow(d)
    m <- ncol(d)
    k.tot <- ncol(z)

    # to avoid floating point precision issues, use a tolerance of 1e-10 to
    # test for equality
    eps <- 10

    result <- NULL

    if (k.tot == 1) {
        # there is only one cluster! we cannot collect units anywhere!!
        return(result)
    }

    clusters <- which(n > 0)

    # set of merged clusters
    m.set <- rep(FALSE, k.tot)

    idx.w <- partition == w

    log.phi <- 0
    log.prob.tot <- 0

    increment.best <- 0
    increment.old <- 0
    ac.old <- -Inf

    k.new <- k.tot
    counter <- 1

    done <- FALSE
    while (!done) {
        log.p.k <- ComputeLogCondProb(d[!idx.w, , drop=FALSE],
                                      partition[!idx.w], n, n1s, aa, bb)

        log.p.w <- ComputeLogCondProbW(d[!idx.w, , drop=FALSE],
                                       n[w], n1s[, w, drop=FALSE],
                                       aa[, w, drop=FALSE],
                                       bb[, w, drop=FALSE])

        log.probs <- log.p.w - log.p.k

        i.max <- which.max(log.probs)
        i <- which(!idx.w)[i.max]

        # find the cluster of i.max
        k <- partition[i]

        if (n[k] == 1) {
            # moving "i" from cluster "k" to cluster "w" will decrease the
            # total number of clusters! We have to compute the bonus (or
            # penalty) factor for moving from k.tot to k.new
            k.new <- k.new - 1

            m.set[clusters == k] <- TRUE

            # psi_{j,h,s} = omega_{j,h,s}^{(new)} / omega_{j,h,s}^{(old)}
            log.psi <- LogOProbs(m, k.new) - log.op

            idx.lin.not.merged <- as.numeric((1:m) +
                                             m * (z[, !m.set, drop=FALSE] - 1))

            t1 <- ((k.tot - k.new) * log.gp[idx.lin.not.merged]) /
                  (k.tot * k.new) + log.psi[idx.lin.not.merged]

            idx.lin.merged <- as.numeric((1:m) +
                                         m * (z[, m.set, drop=FALSE] - 1))

            t2 <- (log.gp[idx.lin.merged] / k.tot) + log.op[idx.lin.merged]

            log.phi <- sum(t1) - sum(t2)
        }

        log.prob.tot <- log.prob.tot + log.probs[i.max]
        increment <- log.phi + log.prob.tot

        # update info
        partition[i] <- w
        idx.w[i] <- TRUE

        if (round(increment - increment.best, eps) > 0) {
            # we got a better solution!
            result <- partition
        }

        counter <- counter + 1

        if (!all(idx.w)) {
            n[k] <- n[k] - 1
            n[w] <- n[w] + 1

            n1s[, k] <- n1s[, k] - d[i, ]
            n1s[, w] <- n1s[, w] + d[i, ]

            if (T == 0) {
                # if the increment is already positive (better solution
                # than the original), keep going no matter what
                #
                # stop the search if the "curve is not convex"
                # if ac.new > 0 then we are going "up", which is ok
                # if ac.new < 0 then we are going "down". In this case we
                #                    have to check if ac.new > ac.old
                # if ac.new = 0 then we want to go on and see what happens
                #                    next
                ac.new <- round(increment - increment.old, eps)

                if (increment < 0 && ac.new < 0 && ac.new < ac.old) {
                    done <- TRUE
                } else {
                    increment.old <- increment
                    ac.old <- ac.new
                }
            } else if (counter > T) {
                done <- TRUE
            }
        } else {
            done <- TRUE
        }
    }

    return(result)
}

###############################################################################
#
# OperatorCreateClusterTypeB
#
# Description:
#
# Try to create a new cluster 'w' by collecting units from any other cluster
#
# Arguments:
#
#         w: integer scalar, representing the (new) cluster receiving the units
#         d: logical data matrix
# partition: integer vector (length n) of cluster indices
#         n: integer vector (length max.clust) with the size of each cluster
#       n1s: m-by-max.clust integer matrix with the number of ones for each
#            column and for each cluster
#         a: m-by-4 double matrix of Beta distribution alpha parameters
#         b: m-by-4 double matrix of Beta distribution beta parameters
#        aa: m-by-max.clust matrix of Beta distribution alpha parameters
#        bb: m-by-max.clust matrix of Beta distribution beta parameters
#         z: m-by-k.tot integer matrix with "properties" and "statuses"
#    log.gp: m-by-4 double matrix of "properties" log probabilities
#    log.op: m-by-4 double matrix of "statuses" log probabilities
#         T: integer scalar, representing the fixed number of units to move.
#            Set 'T=0' to find the units in a "dynamic" approach
#
# Value:
#
# Integer vector with the new partition or NULL if no better solution has been
# found
#
###############################################################################
OperatorCreateClusterTypeB <- function(w, d, partition, n, n1s, a, b, aa, bb,
                                       z, log.gp, log.op, T) {
    v <- nrow(d)
    m <- ncol(d)
    k.tot <- ncol(z)

    # to avoid floating point precision issues, use a tolerance of 1e-10 to
    # test for equality
    eps <- 10

    result <- NULL

    # remove clusters of size 1 from computation
    idx <- partition %in% which(n > 1)

    if (!any(idx)) {
        return(result)
    }

    idx.w <- rep(FALSE, v)

    log.op.new <- LogOProbs(m, k.tot + 1)

    idx.lin <- as.numeric((1:m) + m * (z - 1))
    const <- sum(log.op.new[idx.lin] - log.op[idx.lin] -
                 (log.gp[idx.lin] / (k.tot * (k.tot + 1))))

    # keep track of the sum of log conditional probabilities
    lcp <- 0

    increment.best <- 0
    increment.old <- 0
    ac.old <- -Inf

    counter <- 1

    done <- FALSE
    while (!done) {
        log.p.k <- lcp + ComputeLogCondProb(d[idx, , drop=FALSE],
                                            partition[idx], n, n1s, aa, bb)

        log.p.w <- ComputeLogCondProbZw(n[w] + 1,
                                        t(d[idx, , drop=FALSE]) + n1s[, w],
                                        a, b, z, log.gp, log.op.new)

        increment <- const + log.p.w - log.p.k

        # update info
        i.max <- which.max(increment)
        i <- which(idx)[i.max]

        lcp <- log.p.k[i.max]

        # find the cluster of i.max
        k <- partition[i]

        n[k] <- n[k] - 1
        n[w] <- n[w] + 1

        n1s[, k] <- n1s[, k] - d[i, ]
        n1s[, w] <- n1s[, w] + d[i, ]

        partition[i] <- w
        idx[i] <- FALSE
        idx.w[i] <- TRUE

        if (round(increment[i.max] > increment.best, eps) > 0) {
            # we got a better solution!
            increment.best <- increment[i.max]
            result <- partition
        }

        if (n[k] == 1) {
            idx[partition == k] <- FALSE
        }

        counter <- counter + 1
        
        if (any(idx)) {
            if (T == 0) {
                # if the increment is already positive (better solution
                # than the original), keep going no matter what
                #
                # stop the search if the "curve is not convex"
                # if ac.new > 0 then we are going "up", which is ok
                # if ac.new < 0 then we are going "down". In this case we
                #                    have to check if ac.new > ac.old
                # if ac.new = 0 then we want to go on and see what happens
                #                    next
                ac.new <- round(increment[i.max] - increment.old, eps)

                if (increment[i.max] < 0 && ac.new < 0 && ac.new < ac.old) {
                    done <- TRUE
                } else {
                    increment.old <- increment[i.max]
                    ac.old <- ac.new
                }
            } else if (counter > T) {
                done <- TRUE
            }
        } else {
            done <- TRUE
        }
    }

    return(result)
}

###############################################################################
#
# OperatorSwitchUnits
#
# Description:
#
# Try to 'switch' units from different clusters
#
# Arguments:
#
#         d: logical data matrix
# partition: integer vector (length n) of cluster indices
#         n: integer vector (length max.clust) with the size of each cluster
#       n1s: m-by-max.clust integer matrix with the number of ones for each
#            column and for each cluster
#        aa: m-by-max.clust matrix of Beta distribution alpha parameters
#        bb: m-by-max.clust matrix of Beta distribution beta parameters
#
# Value:
#
# Integer vector with the new partition or NULL if no better solution has been
# found
#
# Note:
#
# --TODO-- To be optimized
#
###############################################################################
OperatorSwitchUnits <- function(d, partition, n, n1s, aa, bb) {
    v <- nrow(d)
    m <- ncol(d)

    # to avoid floating point precision issues, use a tolerance of 1e-10 to
    # test for equality
    eps <- 10

    idx.c <- n > 0

    result <- NULL

    if (sum(idx.c) == 1) {
        # there is only one cluster! we cannot switch units!!
        return(result)
    }

    f <- function(i) {
        idx0 <- partition[-(1:i)] != partition[i]
        idx1 <- which(idx0) + i

        tmp <- matrix(0, nrow=length(idx1), ncol=m)
        tmp[,  d[i, ]] <- log( a[idx1,  d[i, ], drop=FALSE] +
                              n1[idx1,  d[i, ], drop=FALSE])
        tmp[, !d[i, ]] <- log( b[idx1, !d[i, ], drop=FALSE] + nn[idx1] -
                              n1[idx1, !d[i, ], drop=FALSE])

        p.iw <- rowSums(tmp - log( a[idx1, , drop=FALSE] +
                                   b[idx1, , drop=FALSE] +
                                  nn[idx1]))

        x1 <- log(a[i, ] + n1[i, ])
        x2 <- log(b[i, ] + nn[i] - n1[i, ])

        D <- d[idx1, , drop=FALSE]

        tmp <- matrix(0, nrow=length(idx1), ncol=m)

        idx2 <- floor((which(D) - 1) / nrow(tmp)) + 1
        tmp[idx2] <- x1[idx2]

        idx2 <- floor((which(!D) - 1) / nrow(tmp)) + 1
        tmp[idx2] <- x2[idx2]

        p.jk <- rowSums(tmp) - sum(log(a[i, ] + b[i, ] + nn[i]))

        res <- numeric(length(idx0))
        res[ idx0] <- p.jk + p.iw - dia[idx1] - dia[i]
        res[!idx0] <- -Inf
        names(res) <- (i+1):v

        return(res)
    }

    increment <- 0
    increment.best <- 0
    increment.old <- 0

    done <- FALSE
    while (!done) {
        # v-by-m matrices
        a <- t(aa[, partition, drop=FALSE])
        b <- t(bb[, partition, drop=FALSE])

        # total number of units in each cluster, after a unit has been removed
        # length(nn) = v
        nn <- n[partition] - 1

        # total number of ones for each column and for each cluster, after a
        # unit has been removed. v-by-m matrix
        n1 <- t(n1s[, partition, drop=FALSE]) - d

        # diagonal
        tmp <- matrix(0, nrow=v, ncol=m)
        tmp[ d] <- log(a[d] + n1[d])
        tmp[!d] <- log(b[!d] + (nn - n1)[!d])
        dia <- rowSums(tmp - log(a + b + nn))

        # lower triangular matrix
        ltm <- lapply(X=1:(v-1), FUN=f)

        i.unit <- which.max(unlist(lapply(X=ltm, FUN=max)))
        j.unit <- as.integer(names(which.max(ltm[[i.unit]])))

        n1s[, partition[i.unit]] <- n1s[, partition[i.unit]] - d[i.unit, ]
        n1s[, partition[i.unit]] <- n1s[, partition[i.unit]] + d[j.unit, ]

        n1s[, partition[j.unit]] <- n1s[, partition[j.unit]] - d[j.unit, ]
        n1s[, partition[j.unit]] <- n1s[, partition[j.unit]] + d[i.unit, ]

        tmp <- partition[i.unit]
        partition[i.unit] <- partition[j.unit]
        partition[j.unit] <- tmp

        increment <- increment + max(ltm[[i.unit]])

        # to avoid floating point precision issues, use an epsilon of 1e-10 to
        # compare the numbers. If the difference is truly of 1e-10 on
        # logarithmic scale... is it really a better solution to care about?
        if (round(increment - increment.best, eps) > 0) {
            # we got a better solution!
            increment.best <- increment
            result <- partition
        }

        # continue only if we get a strict increase
        if (round(increment - increment.old, eps) > 0) {
            increment.old <- increment
        } else {
            done <- TRUE
        }
    }

    return(result)
}
