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
