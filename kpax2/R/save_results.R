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
# SaveResults
#
# Description:
#
# Write the result into several text files
#
# Arguments:
#
#   full.data: a list containing the following variables
#                    data: original data encoded as a logical (binary) matrix
#                  values: character vector with unique values per column
#                    keys: numeric vector with indices of each values
#                     uid: n-by-1 character vector with the sample units' ids
#                    poly: M-by-1 boolean vector indicating where the
#                          polymorphic variables where found in the original
#                          dataset
#      result: a list containing the following variables
#                  partition: vector of cluster indices corresponding to the
#                             (local) optimum. Clusters are sorted in
#                             descending order according to their size
#                  max.logpp: value of the log posterior probability
#                    idx.max: m-by-k.tot integer matrix that represents not
#                             informative, informative and characteristic
#                             attributes, where m is the total number of
#                             columns while k.tot is the total number of
#                             clusters
# output.file: character string representing a path to an output file. Several
#              output files are produced (depending on the value of variable
#              'what'), so no file extension is required and the file name
#              should represent the common base name
#        what: integer scalar.
#                  1: basic summary and partition only
#                  2: basic summary, partition and column classification
#                  3: basic summary, partition, column classification and
#                     characteristic features
#                  4: basic summary, partition, column classification,
#                     characteristic features and original data
#
# Value:
#
# Text files containing K-Pax2 output
#
###############################################################################
SaveResults <- function(full.data, result, output.file, what=1) {
    if (is.character(output.file)) {
        if (length(output.file) != 1) {
            stop("Variable 'output.file' is not a string")
        }

        finfo <- file.info(output.file)

        if (!is.na(finfo$isdir) && finfo$isdir) {
            stop(paste("'", output.file, "' is a directory", sep=""))
        }
    } else {
        stop("Variable 'output.file' is not a string")
    }

    if (is.numeric(what)) {
        if (length(what) != 1) {
            stop("Variable 'what' is not a scalar")
        }

        what <- as.integer(what + 0.5)

        if (what <= 0) {
            what <- 1
        }
    } else {
        stop("Variable 'what' is not of numeric type")
    }

    clusters <- sort(unique(result$partition))

    v <- nrow(full.data$data)
    m <- ncol(full.data$data)
    M <- length(full.data$poly)

    k.tot <- length(clusters)

    out.file.1.a <- paste(output.file, "_partition.csv", sep="")
    out.file.1.b <- paste(output.file, "_summary.txt", sep="")

    write.table(result$partition, file=out.file.1.a, row.names=FALSE,
                col.names=FALSE)

    write(format(result$max.logpp, nsmall=6), file=out.file.1.b, append=FALSE)

    # -- TODO --
    # We should also add the hyperparameters that have been used to run the
    # analysis

    cat("Partition has been written to file:", out.file.1.a, "\n")
    cat("Summary has been written to file:", out.file.1.b, "\n")

    if (what >= 2) {
        f <- function(j) {
            tmp <- result$idx.max[full.data$keys == j, , drop=FALSE] == 4

            x <- apply(tmp, 2, any)
            y <- rep(1, k.tot)

            if (any(x)) {
                y <- ifelse(x, 4, 3)
            } else {
                tmp <- result$idx.max[full.data$keys == j, , drop=FALSE] == 2
                x <- apply(tmp, 2, any)

                if (any(x)) {
                    y <- rep(2, k.tot)
                }
            }

            return(y)
        }

        complete.idx.max <- matrix(1, nrow=M, ncol=k.tot)
        complete.idx.max[full.data$poly, ] <- t(vapply(1:max(full.data$keys),
                                                 FUN=f,
                                                 FUN.VALUE=numeric(k.tot),
                                                 USE.NAMES=FALSE))

        out.file.2 <- paste(output.file, "_attributes.csv", sep="")
        write.table(complete.idx.max, file=out.file.2, sep=",",
                    row.names=FALSE, col.names=1:ncol(complete.idx.max))

        cat("Attributes have been written to file:", out.file.2, "\n")
    }

    if (what >= 3) {
        out.file.3 <- paste(output.file, "_characteristic.csv", sep="")

        tmp <- result$idx.max == 4

        idx.lin <- which(tmp)
        i.idx <- ((idx.lin - 1) %% nrow(tmp)) + 1
        j.idx <- floor((idx.lin - 1) / nrow(tmp)) + 1

        col.idx <- full.data$keys[i.idx]
        unique.idx <- sort(unique(col.idx))

        to.save <- matrix("", nrow=length(unique.idx), ncol=k.tot)

        for (i in 1:length(unique.idx)) {
            j <- which(col.idx == unique.idx[i])

            values <- toupper(full.data$values[i.idx[j]])
            values[values == "-"] <- "+"

            to.save[i, j.idx[j]] <- values
        }

        write.table(cbind(as.character(which(full.data$poly)[unique.idx]),
                    to.save), file=out.file.3, sep=",", row.names=FALSE,
                    col.names=c("Site", paste("Cluster", 1:k.tot)))
    }

    if (what >= 4) {
        out.file.4 <- paste(output.file, "_info.txt", sep="")

        orig.data <- matrix(" ", nrow=M, ncol=v)

        p.idx <- which(full.data$poly)
        values <- tolower(full.data$values)
        for (i in 1:v) {
            # there could be missing values within different sequences. This is
            # why we handle each sequence separately
            un <- full.data$data[i, ]
            r.idx <- p.idx[full.data$keys[un]]
            orig.data[r.idx, i] <- values[un]
        }

        for (k in clusters) {
            idx.k <- result$partition == k

            attr.str <- rep(" ", length=M)
            d.k <- t(orig.data[, idx.k, drop=FALSE])

            # characteristic sites for cluster k
            cs.k <- result$idx.max[, k] == 4

            # characteristic values observed at these sites
            cv.k <- full.data$values[cs.k]

            # what are their positions in the data?
            cp.k <- full.data$keys[cs.k]

            # check if a column has more than a characteristic attribute
            co.k <- tabulate(cp.k, nbins=max(full.data$keys))

            # only one characteristic attribute
            idx.a <- co.k == 1
            idx.b <- cp.k %in% which(idx.a)

            idx.c <- which(full.data$poly)[idx.a]

            idx.d <- (d.k[, idx.c, drop=FALSE] ==
                      rep(cv.k[idx.b], each=nrow(d.k)))

            d.k[, idx.c][idx.d] <- toupper(d.k[, idx.c][idx.d])
            d.k[, idx.c][idx.d][d.k[, idx.c][idx.d] == "-"] <- "+"

            attr.str[idx.c] <- toupper(cv.k[idx.b])

            # more than a characteristic attribute
            idx.a <- co.k > 1

            if (any(idx.a)) {
                cset <- which(idx.a)
                idx.b <- cp.k %in% cset

                for (j in cset) {
                    idx.j <- full.data$keys == j

                    cs.jk <- cs.k & idx.j
                    cv.jk <- cv.k[idx.b]

                    idx.c <- which(full.data$poly)[j]
                    idx.d <- d.k[, idx.c, drop=FALSE] %in% cv.jk

                    d.k[, idx.c][idx.d] <- toupper(d.k[, idx.c][idx.d])
                    d.k[, idx.c][idx.d][d.k[, idx.c][idx.d] == "-"] <- "+"

                    attr.str[idx.c] <- as.character(co.k[idx.c])
                }
            }

            w <- matrix(pmax(nchar(attr.str),
                             apply(d.k, 2, function(x) max(nchar(x)))),
                        nrow=nrow(d.k), ncol=ncol(d.k), byrow=TRUE)

            new.data <- format(d.k, justify="right", width=w)
            new.data <- apply(new.data, 1, paste, collapse="")
            new.data <- cbind(new.data, full.data$uid[idx.k])
            new.data <- apply(new.data, 1, paste, collapse=" ; ")

            write(paste(format(attr.str, justify="right", width=w[1, ]),
                        collapse=""), file=out.file.4, append=TRUE)
            write(new.data, file=out.file.4, append=TRUE)
            write("", file=out.file.4, append=TRUE)
        }

        cat("Output has been written to file:", out.file.4, "\n")
    }
}
