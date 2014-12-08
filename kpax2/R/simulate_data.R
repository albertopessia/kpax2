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
# SimulateDataset
#
# Description:
#
# Generate a random dataset to be analysed by K-Pax2
#
# Arguments:
#
#     nk: integer vector of cluster sizes. It is assumed that the length of nk
#         is equal to the total number of clusters while the sum of nk is equal
#         to the sample size. Floating point numbers are rounded to integers
#     nc: cardinality of the alphabet, i.e. number of categories (at least 2)
# n.attr: integer vector of length 3. n.attr[j] should represent the total
#         number of columns classified as j, where j is 1 (noise), 2 (weak
#         signal) and 3 (strong signal). sum(n.attr) is equal to M, the total
#         number of columns
#     op: probability that a feature with property 3 will be characteristic for
#         the considered cluster. If NA it will be used the default value 1/K,
#         where K is the total number of clusters
#
# Value:
#
# A list containing the following arguments:
#         D: character n-by-M matrix
# partition: vector of cluster indices
# col.class: integer K-by-M matrix with column classification
#
###############################################################################
SimulateDataset <- function(nk, nc, n.attr, op=NA) {
    if (is.numeric(nk)) {
        if (length(nk) == 0) {
            stop("Variable 'nk' has length 0")
        }

        nk <- as.integer(nk + 0.5)

        if (any(nk <= 0)) {
            stop("Variable 'nk' must contain only positive numbers")
        }
    } else {
        stop("Variable 'nk' is not of numeric type")
    }

    n <- sum(nk)
    K <- length(nk)

    if (is.numeric(nc)) {
        if (length(nc) != 1) {
            stop("Variable 'nc' is not a scalar")
        }

        nc <- as.integer(nc + 0.5)

        if (nc < 2) {
            stop("Variable 'nc' is lesser than 2")
        }
    } else {
        stop("Variable 'nc' is not of numeric type")
    }

    if (is.numeric(n.attr)) {
        if (length(n.attr) != 3) {
            stop("Variable 'n.attr' does not have length 3")
        }

        n.attr <- as.integer(n.attr + 0.5)

        if (any(n.attr < 0)) {
            stop("Variable 'n.attr' cannot contain negative numbers")
        }

        if (all(n.attr == 0)) {
            stop("Variable 'n.attr' must contain at least a positive number")
        }

        if (n.attr[3] > 0 & (n.attr[3] * (nc - 1)) < K) {
            stop(paste("By definition, each cluster should have at least a",
                       "characteristic feature. Increase either n.attr[3] or",
                       "nc."))
        }
    } else {
        stop("Variable 'n.attr' is not of numeric type")
    }

    M <- sum(n.attr)

    if (is.numeric(op)) {
        if (length(op) != 1) {
            stop("Variable 'op' is not a scalar")
        }

        if (op <= 0 || op >= 1) {
            stop("Variable 'op' must be in the set (0, 1)")
        }
    } else if (!is.na(op)) {
        stop("Variable 'op' is not of numeric type, nor NA")
    }

    partition <- rep(1:K, times=nk)
    alphabet <- as.character(1:nc)
    a <- log(0.001) / log(0.95)

    if (is.na(op)) {
        op <- 1 / K
    }

    # legend
    # 1: noise
    # 2: weak signal
    # 3: strong signal but not characteristic for this cluster
    # 4: strong signal and characteristic for this cluster
    col.class <- matrix(0, nrow=M, ncol=K)

    # columns are considered independent from each other, so it does not really
    # matter if column classification is random or not, as long as their number
    # is equal to the one given by n.attr
    lb <- 1

    if (n.attr[1] > 0) {
        idx.1 <- lb:n.attr[1]
        lb <- n.attr[1] + 1
    } else {
        idx.1 <- integer(0)
    }

    if (n.attr[2] > 0) {
        idx.2 <- lb:(n.attr[1] + n.attr[2])
        lb <- n.attr[1] + n.attr[2] + 1
    } else {
        idx.2 <- integer(0)
    }

    if (n.attr[3] > 0) {
        idx.3 <- lb:M
    } else {
        idx.3 <- integer(0)
    }

    col.class[idx.1, ] <- 1
    col.class[idx.2, ] <- 2
    col.class[idx.3, ] <- 3

    D <- matrix("", nrow=n, ncol=M)

    f <- function(x) {
        p <- rgamma(x, 1, 1)
        return(p / sum(p))
    }

    g <- function(p, n) {
        nv <- round(p * n)
        v.set <- 1:length(nv)

        dif <- n - sum(nv)
        while (dif != 0) {
            if (dif > 0) {
                tmp <- sample(v.set, size=dif, replace=TRUE)
            } else {
                # we need to be careful here to not go below zero
                # easy but dirty, we do it one by one (shouldn't abs(dif) be
                # one anyway?)
                tmp <- sample(v.set[nv > 0], size=1)
            }

            nv <- nv + sign(dif) * tabulate(bin=tmp, nbins=length(nv))
            dif <- n - sum(nv)
        }

        return(sample(rep(alphabet[1:length(nv)], times=nv)))
    }

    # simulate probabilities from a Dirichlet distribution Dir(1, 1, ..., 1)
    if (n.attr[1] > 0) {
        obs.noise <- sample(2:nc, size=n.attr[1], replace=TRUE)
        p.noise <- lapply(obs.noise, f)

        # instead of using p.noise[[j]] as probabilities for column j, we
        # assume they are the observed proportion for that column
        #
        # This will guarantee that the number of categories for column j
        # will be exactly equal to obs.noise[j], as requested
        #
        # Sampling from a multinomial with prob p.noise[[j]] doesn't
        # guarantee it
        D[, idx.1] <- unlist(lapply(p.noise, g, n))
    }

    if (n.attr[2] > 0) {
        obs.weak <- sample(2:nc, size=n.attr[2], replace=TRUE)

        for (k in 1:K) {
            idx.units <- (partition == k)
            p.weak <- lapply(obs.weak, f)
            D[idx.units, idx.2] <- unlist(lapply(p.weak, g, nk[k]))
        }
    }

    if (n.attr[3] > 0) {
        # first of all, we need to guarantee that each cluster has at least
        # one characteristic attribute (by definition)
        got.fingerprint <- rep(FALSE, K)

        column.counter <- 1
        category.counter <- 1
        for (k in 1:K) {
            idx.k <- (partition == k)

            while (!got.fingerprint[k]) {
                # do we still have enough categories available?
                if (category.counter < nc) {
                    p.charac <- rbeta(1, a, 1)
                    s.set <- alphabet[alphabet != alphabet[category.counter]]
                    c.idx <- rbinom(nk[k], size=1, prob=p.charac) == 1

                    v <- sample(s.set, size=nk[k], replace=TRUE)
                    v[c.idx] <- alphabet[category.counter]

                    col.class[idx.3[column.counter], k] <- 4
                    D[idx.k, idx.3[column.counter]] <- v

                    got.fingerprint[k] <- TRUE
                    category.counter <- category.counter + 1
                } else {
                    # move on

                    # units which have not yet been assigned in this column
                    idx.u <- D[, idx.3[column.counter]] == ""

                    p.charac <- rbeta(1, a, 1)
                    s.set <- alphabet[alphabet != alphabet[nc]]
                    c.idx <- rbinom(sum(idx.u), size=1, prob=p.charac) == 1

                    v <- sample(s.set, size=sum(idx.u), replace=TRUE)
                    v[c.idx] <- alphabet[nc]

                    D[idx.u, idx.3[column.counter]] <- v

                    # it could happen that all these units are in the last
                    # cluster, in which case we have just assigned a
                    # characteristic category also for the last group
                    if (all(idx.u & (partition == K))) {
                        col.class[idx.3[column.counter], K] <- 4
                        got.fingerprint[K] <- TRUE
                    }

                    column.counter <- column.counter + 1
                    category.counter <- 1

                    if (column.counter > length(idx.3)) {
                        stop(paste("Not enough strong signal columns nor",
                                   "categories to fullfill the requirement",
                                   "of at least a characteristic site per",
                                   "cluster"))
                    }
                }
            }
        }

        # finish to fill in the column
        if (category.counter > 1) {
            idx.u <- D[, idx.3[column.counter]] == ""

            pr <- c(rbeta(category.counter - 1, 1, a),
                    rbeta(nc - category.counter + 1, a, 1))
            pr <- pr / sum(pr)

            v <- sample(alphabet, size=sum(idx.u), replace=TRUE, prob=pr)

            D[idx.u, idx.3[column.counter]] <- v

            column.counter <- column.counter + 1
            category.counter <- 1
        }

        idx.3 <- idx.3[-(1:(column.counter - 1))]

        # the total number of characteristic properties is distributed as a
        # binomial distribution (sum of K independent bernoulli r.v.)
        # we want to have at least 1 for each column classified as "3", i.e.
        # a truncated binomial
        p.fp <- dbinom(1:K, size=K, prob=op)
        p.fp <- p.fp / sum(p.fp)

        fp <- sample(1:K, size=length(idx.3), replace=TRUE, prob=p.fp)

        for (j in 1:length(idx.3)) {
            nfp <- fp[j]

            p.charac <- rbeta(1, a, 1)

            pr <- c(rbeta(1, 1, a), rbeta(nc - 1, a, 1))
            pr <- pr / sum(pr)

            v <- sample(alphabet, size=n, replace=TRUE, prob=pr)

            clust.set <- sort(sample(1:K, size=nfp, replace=FALSE))
            col.class[idx.3[j], clust.set] <- 4

            for (i in clust.set) {
                c.idx <- rbinom(nk[i], size=1, prob=p.charac) == 1
                v[partition == i][c.idx] <- alphabet[1]
            }

            D[, idx.3[j]] <- v
        }
    }

    return(list(D=D, partition=partition, col.class=col.class))
}

###############################################################################
#
# SaveSimulatedDataset
#
# Description:
#
# Write the values returned by SimulateDataset function into text files
#
# Arguments:
#
#         sim: list object as returned by SimulateDataset function
# output.file: character string representing a path to an output file. Several
#              output files are produced, so no file extension is required and
#              the file name should represent the common base name
#
# Value:
#
# Text files containing the simulated dataset
#
###############################################################################
SaveSimulatedDataset <- function(sim, output.file) {
    out.file.attributes <- paste(output.file, "_attributes.csv", sep="")
    write.table(sim$col.class, file=out.file.attributes, sep=",",
                row.names=FALSE, col.names=1:ncol(sim$col.class))

    out.file.partition <- paste(output.file, "_partition.csv", sep="")
    write.table(sim$partition, file=out.file.partition, row.names=FALSE,
                col.names=FALSE)

    out.file.data <- paste(output.file, "_data.csv", sep="")
    write.table(sim$D, file=out.file.data, sep=",", quote=TRUE,
                row.names=FALSE, col.names=FALSE)
}
