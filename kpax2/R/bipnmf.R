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
# BiPNMF
#
# Description:
#
# The Bi-PNMF algorithm
#
# Arguments:
#
#            X: data matrix (nonnegative)
#            W: initial W matrix
#            H: initial H matrix
# row.clusters: number of row clusters
# col.clusters: number of column clusters
#     max.iter: maximum number of iterations
#       t.iter: integer scalar. Print a status message every 't.iter'
#               iterations. Set 't.iter <= 0' to disable status messages
#          tol: convergence tolerance
#
# Value:
# A list containing updated matrices W and H
#
# Details:
#
# To get row clusters, just pick the maximum entry of each row of W
# To get column clusters, just pick the maximum entry of each row of H
#
# Bi-PNMF solves the following approximation problem
# minimize_{W >= 0, H >= 0}  ||X - W W' X H H' ||_F
#
###############################################################################
BiPNMF <- function(X, W=NULL, H=NULL, row.clusters=ncol(W),
                   col.clusters=ncol(H), max.iter=1000, t.iter=100, tol=1e-8) {
  if (is.null(W)) {
    W <- matrix(runif(nrow(X) * row.clusters), nrow=nrow(X))
  }

  if (is.null(H)) {
    H <- matrix(runif(ncol(X) * col.clusters), nrow=ncol(X))
  }

  trXX <- norm(X, "F")^2

  WtW <- t(W) %*% W
  WtXH <- t(W) %*% X %*% H
  HtH <- t(H) %*% H
  XH <- X %*% H

  # ||X - W W' X H H'||_{F}^{2}
  new.value <- trXX -
               2 * norm(WtXH, "F")^2 +
               sum(diag((HtH %*% t(WtXH) %*% WtW %*% WtXH)))

  if (t.iter > 0) {
    cat("Initial difference =", new.value, "\n")
  }

  eps <- .Machine$double.eps
  keepgoing <- TRUE
  iter <- 0

  while (keepgoing) {
    iter <- iter + 1
    old.value <- new.value

    W <- W * ((2 * XH %*% t(WtXH)) /
              (W %*% WtXH %*% HtH %*% t(WtXH) +
               XH %*% HtH %*% t(WtXH) %*% WtW +
               eps))^(0.25)

    XtW <- t(X) %*% W
    WtXH <- t(XtW) %*% H
    WtW <- t(W) %*% W

    ## update H
    H <- H * ((2 *XtW %*% WtXH) /
              (XtW %*% WtW %*% WtXH %*% HtH +
               H %*% t(WtXH) %*% WtW %*% WtXH +
               eps))^(0.25)

    XH <- X %*% H
    WtXH <- t(W) %*% XH
    HtH <- t(H) %*% H

    new.value <- trXX -
                 2 * norm(WtXH, "F")^2 +
                 sum(diag((HtH %*% t(WtXH) %*% WtW %*% WtXH)))

    if (abs(new.value - old.value) < tol * abs(old.value)) {
      if (t.iter > 0) {
        cat("Converged after", iter, "steps.\n")
      }

      keepgoing <- FALSE
    } else if (t.iter > 0 && (iter %% t.iter) == 0) {
      cat("Iteration", iter, "done. Current difference =", new.value, "\n")
    }

    if (iter == max.iter) {
      keepgoing <- FALSE
    }
  }

  return(list("W"=W, "H"=H))
}

###############################################################################
#
# SVDInit
#
# Details:
#
# Bi-PNMF initialization using SVD. This function is especially needed when
# the number of rows is much less than the number of columns
#
# Arguments:
#
#            X: data matrix
# row.clusters: number of row clusters
# col.clusters: number of column clusters
#
# Value:
#
# A list containing W and H matrices
#
###############################################################################
SVDInit <- function(X, row.clusters, col.clusters) {
  svd.W <- svd(X, nu=row.clusters, nv=row.clusters)
  ord <- order(svd.W$d, decreasing=TRUE)
  W <- Discretize(svd.W$u[, ord[1:row.clusters], drop=FALSE])

  svd.H <- svd(t(X), nu=col.clusters, nv=col.clusters)
  ord <- order(svd.H$d, decreasing=TRUE)
  H <- Discretize(svd.H$u[, ord[1:col.clusters], drop=FALSE])

  return(list(W=W, H=H))
}

###############################################################################
#
# Discretize
#
# Description:
#
# Project a non-negative matrix to a cluster indicator matrix, i.e. each row
# has only one 1-entry and zero elsewhere.
#
# Arguments:
# Y: Matrix to be discretized
#
# Value:
# Z: Discretized matrix (i.e. cluster indicator matrix)
#
# Details:
#
# Finds a cluster indicator matrix Z that (locally) minimizes ||Z - Y||_F
# See "Multiclass Spectral Clustering" by Stella Yu and Jianbo Shi, 2003
#
###############################################################################
Discretize <- function(Y) {
  n <- nrow(Y)
  m <- ncol(Y)

  # normalize each row of Y
  # sqrt(diag(Y %*% t(Y)))
  vm <- sqrt(rowSums(Y^2))
  idx <- vm > 0
  Y[idx, ] <- Y[idx, ] / vm[idx]

  R <- matrix(0, m, m)
  R[, 1] <- Y[sample(1:n, size=1), ]

  if (m > 1) {
    c <- rep(0, n)
    for (j in 2:m) {
      c <- c + abs(Y %*% R[, j - 1])
      R[, j] <- Y[which.min(c), ]
    }
  }

  Z <- matrix(0, nrow=n, ncol=m)

  old.value <- 0
  keepgoing <- TRUE
  iter <- 0
  iter.max <- 1000

  while (keepgoing) {
    iter <- iter + 1

    tmp1 <- Y %*% R

    idx <- max.col(tmp1, ties.method="first")

    Z[, ] <- 0
    Z[(1:n) + n * (idx - 1)] <- 1

    tmp2 <- svd(t(Z) %*% Y)

    new.value <- sum(tmp2$d)

    if (abs(new.value - old.value) < .Machine$double.eps ||
        iter > iter.max) {
      keepgoing <- FALSE
    } else {
      old.value <- new.value
      R <- tmp2$v %*% t(tmp2$u)
    }
  }

  return(Z)
}
