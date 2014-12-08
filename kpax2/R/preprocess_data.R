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
# LoadData
#
# Description:
#
# Read data in FASTA or CSV format and convert it to K-Pax2 format. Only
# categorical data is assumed. Sequences (if FASTA format) are required to be
# aligned. Data (if CSV format) must be representable as a regular data matrix
# (same number of variables for each unit). Only polymorphic columns are stored
# while their attributes are converted to integer numbers (if not already
# integers)
#
# Arguments:
#
#    input.file: path to input data file
#   output.file: path to output file. Use NULL if you don't want to write the
#                output to a Rdata file. Default value: NULL
#      filetype: string representing the data format to read. Use either "csv"
#                or "fasta"
#        header: logical value. Used only if filetype="csv". This indicates
#                whether the file contains the names of the variables as its
#                first line
#           ids: logical value. Used only if filetype="csv". This indicates
#                whether the file contains the ids of the units as its first
#                column. If ids=FALSE, they will be labeled by the numbers
#                1 to n
#          miss: character vector with the values to be considered missing.
#                Use NULL if no missing value is expected/considered. For a
#                FASTA file, only single characters (length 1) are allowed.
#                Typical values for 'miss' are:
#                    - DNA data
#                      miss <- c("-", "?", "*", "#", "b", "d", "h", "k", "m",
#                                "n", "r", "s", "v", "w", "x", "y")
#                    - Protein data
#                      miss <- c("-", "?", "*", "#", "b", "j", "x", "z")
# lines.to.read: The (maximal) number of lines to read at each iteration.
#                Should be used when reading a very large dataset. Negative
#                values indicate that one should read up to the end of input on
#                the connection in one single step.
#
# Details:
#
# Standard conversion tables for FASTA data
#
#            Conversion table (DNA)
# ------------------------------------------
# |          Nucleotide |  Code |  Integer |
# ------------------------------------------
# |           Adenosine |   A   |     1    |
# |            Cytosine |   C   |     2    |
# |             Guanine |   G   |     3    |
# |           Thymidine |   T   |     4    |
# |              Uracil |   U   |     5    |
# |     Purine (A or G) |   R   |     6    |
# | Pyrimidine (C or T) |   Y   |     7    |
# |                Keto |   K   |     8    |
# |               Amino |   M   |     9    |
# |  Strong Interaction |   S   |    10    |
# |    Weak Interaction |   W   |    11    |
# |               Not A |   B   |    12    |
# |               Not C |   D   |    13    |
# |               Not G |   H   |    14    |
# |          Not T or U |   V   |    15    |
# |                 Any |   N   |    16    |
# |              Masked |   X   |    17    |
# |                 Gap |   -   |    18    |
# ------------------------------------------
#
#              Conversion table (PROTEIN)
# --------------------------------------------------
# |                  Amino Acid |  Code |  Integer |
# --------------------------------------------------
# |                     Alanine |   A   |     1    |
# |                    Arginine |   R   |     6    |
# |                  Asparagine |   N   |    16    |
# |               Aspartic acid |   D   |    13    |
# |                    Cysteine |   C   |     2    |
# |                   Glutamine |   Q   |    25    |
# |               Glutamic acid |   E   |    19    |
# |                     Glycine |   G   |     3    |
# |                   Histidine |   H   |    14    |
# |                  Isoleucine |   I   |    21    |
# |                     Leucine |   L   |    22    |
# |                      Lysine |   K   |     8    |
# |                  Methionine |   M   |     9    |
# |               Phenylalanine |   F   |    20    |
# |                     Proline |   P   |    24    |
# |                 Pyrrolysine |   O   |    23    |
# |              Selenocysteine |   U   |     5    |
# |                      Serine |   S   |    10    |
# |                   Threonine |   T   |     4    |
# |                  Tryptophan |   W   |    11    |
# |                    Tyrosine |   Y   |     7    |
# |                      Valine |   V   |    15    |
# | Asparagine or Aspartic acid |   B   |    12    |
# |  Glutamine or Glutamic acid |   Z   |    26    |
# |       Leucine or Isoleucine |   J   |    28    |
# |                         Any |   X   |    17    |
# |            Translation stop |   *   |    27    |
# |                         Gap |   -   |    18    |
# --------------------------------------------------
#
# Value:
#
# A list (and a Rdata file if output.file is not NULL) containing the following
# variables:
#   data: original data matrix encoded as a logical (binary) matrix
# values: character vector with unique values per column
#   keys: numeric vector with indices of each values
#    uid: n-by-1 character vector with the sample units' ids
#   poly: M-by-1 boolean vector indicating where the polymorphic variables
#         where found in the original dataset. Note that it doesn't refer to
#         the binary data matrix stored in variable 'data'!
#
# Example:
#
# If original data is made of just the following three rows
#     0 1 2 5
#     2 3 4 2
#     1 2 0 3
# then
#       data = | 0  0  1  0  0  1  0  0  0  1 |
#              | 0  1  0  0  1  0  1  1  0  0 |
#              | 1  0  0  1  0  0  0  0  1  0 |
#     values = c(1, 2, 1, 2, 3, 2, 4, 2, 3, 5) (i.e. 12 123 24 235)
#       keys = c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4) (i.e. 11 222 33 444)
#
# Missing data (0 values) are discarded!
#
# Note:
#
# For FASTA data, it is assumed that the only possible characters are those
# listed in the previous two conversion tables. Unpredictable things could
# happen if unknown non-standard characters are found within the FASTA file.
# Be sure the file is a standard FASTA file. If not, consider to convert the
# file into a CSV file
#
# Thymidine (T) and Uracil (U) are given different codes, 4 and 5 respectively,
# because T and U characters represent different amino acids. If you find both
# in a DNA/RNA dataset, consider to convert the number 5 to 4
#
###############################################################################
LoadData <- function(input.file, output.file=NULL, filetype="fasta",
                     header=FALSE, ids=FALSE, miss=NULL, lines.to.read=-1) {
    if (!is.character(input.file)) {
        stop("Variable 'input.file' is not of type 'character'")
    }

    if (length(input.file) > 1) {
        warning(paste("'input.file' length is greater than one. Only the",
                      "first element will be read"), immediate.=TRUE)

        input.file <- input.file[1]
    }

    if (length(input.file) != 1) {
        stop("'input.file' must be a character vector of length 1")
    }

    finfo <- file.info(input.file)

    if (is.na(finfo$size)) {
        stop("Could not find the input file")
    }

    if (finfo$isdir) {
        stop("The specified path is a directory, not a text file")
    }

    if (finfo$size == 0) {
        stop("The input file is empty")
    }

    if (!(is.character(output.file) || is.null(output.file))) {
        stop("Invalid type for variable 'output.file'")
    }

    filetype <- tolower(filetype)

    if (!(identical(filetype, "fasta") || identical(filetype, "csv"))) {
        stop("'filetype' must be either 'fasta' or 'csv'")
    }

    if (!is.null(miss)) {
        miss <- sort(unique(tolower(miss)))

        if (identical(filetype, "fasta") && any(nchar(miss) != 1)) {
            warning(paste("Variable 'miss' contains non-character missing",
                          "values. Ignoring them"), immediate.=TRUE)

            miss <- miss[nchar(miss) == 1]
        }

        if (length(miss) == 0) {
           miss <- NULL
        }
    }

    if (!is.double(lines.to.read)) {
        lines.to.read <- -1
    } else {
        lines.to.read <- as.integer(lines.to.read + 0.5)

        if (lines.to.read < 1) {
            lines.to.read <- -1
        }
    }

    if (identical(filetype, "fasta")) {
        full.data <- ReadFastaData(input.file, miss, lines.to.read)
    } else {
        if (!is.logical(header)) {
            stop("'header' variable is not of type logical")
        }

        if (!is.logical(ids)) {
            stop("'ids' variable is not of type logical")
        }

        header <- all(header)
        ids <- all(ids)

        full.data <- ReadCsvData(input.file, header, ids, miss, lines.to.read)
    }

    cat("Conversion successful\n")

    if (!is.null(output.file)) {
        save(full.data, file=output.file, ascii=FALSE, compress=TRUE,
             compression_level=6)
    }

    cat("Data has been loaded\n")

    return(full.data)
}

###############################################################################
#
# ReadFastaData
#
# Description:
#
# Read data in FASTA format
#
# Arguments:
#
#    input.file: path to input data file
#          miss: character vector with the letters to be considered missing
#                values. Use NULL if no missing value is expected/considered.
#                Since we are reading a FASTA file, only single characters
#                (length 1) are allowed.
#                Typical values for 'miss' are:
#                    - DNA data
#                      miss <- c("-", "?", "*", "b", "d", "h", "k", "m", "n",
#                                "r", "s", "v", "w", "x", "y")
#                    - Protein data
#                      miss <- c("-", "?", "*", "b", "j", "x", "z")
# lines.to.read: The (maximal) number of lines to read. Negative values
#                indicate that one should read up to the end of input on the
#                connection
#
# Value:
#
# A list containing the following variables:
#   data: original data matrix encoded as a logical (binary) matrix
# values: character vector with unique values per column
#   keys: numeric vector with indices of each values
#    uid: n-by-1 character vector with the sample units' ids
#   poly: M-by-1 boolean vector indicating where the polymorphic variables
#         where found in the original dataset. Note that it doesn't refer to
#         the binary data matrix stored in variable 'data'!
#
###############################################################################
ReadFastaData <- function(input.file, miss, lines.to.read) {
    cat("Loading FASTA data...\n")

    fid <- file(input.file, "rt")
    on.exit(close(fid))

    # we read the file twice
    # the first time we check the total number of sequences, if each sequence
    # has the same length and where are the SNPs
    # the second time we store the data, keeping only the SNPs

    # Note: This function has been written for reading a huge dataset

    # skip blank lines
    done <- FALSE
    while (!done) {
        txt <- Trim(readLines(con=fid, n=1))

        if (!(length(txt) > 0 && identical(txt, ""))) {
            done <- TRUE
        }
    }

    if (length(txt) == 0) {
        stop("Invalid FASTA file. No sequence has been found.")
    }

    if (!identical(substr(txt, 1, 1), ">")) {
        stop(paste("Invalid FASTA file. The first non empty row should be a",
                   "sequence id."))
    }

    # read the first sequence, so that we can pre-allocate enough memory
    # For the first sequence, we preallocate space for a sequence of 10^7 sites
    # Memory required is approximately 80 MB, so it is not really a big problem
    # 10^7 should be enough for almost any dataset
    # We will still check within the loop if it is enough
    seq.len <- 0
    seq.ref <- character(length=1e+07)

    done <- FALSE
    while (!done) {
        txt <- Trim(readLines(con=fid, n=1))

        if (length(txt) > 0 && !identical(substr(txt, 1, 1), ">")) {
            nc <- nchar(txt)

            if (nc > 0) {
                w <- seq.len + nc

                # do we have enough space to store the first sequence?
                if (w > length(seq.ref)) {
                    W <- floor(w / length(seq.ref)) * 1e+07
                    tmp <- seq.ref
                    seq.ref <- character(length=length(tmp)+W)
                    seq.ref[1:length(tmp)] <- tmp
                }

                seq.ref[(seq.len+1):w] <- strsplit(tolower(txt), split="")[[1]]
                seq.len <- w
            }
        } else {
            done <- TRUE
        }
    }

    if (seq.len == 0) {
        stop("Invalid FASTA file. No sequence has been found.")
    }

    seq.ref <- seq.ref[1:seq.len]

    miss.seq.ref <- seq.ref %in% miss

    # at least a sequence has been found
    if (identical(substr(txt, 1, 1), ">")) {
        n <- 2
    } else {
        n <- 1
    }

    poly <- logical(length=seq.len)
    s.tmp <- ""

    while (length(txt <- Trim(readLines(con=fid, n=lines.to.read))) > 0) {
        sidx <- nchar(txt) > 0
        stot <- sum(sidx)
        seq.delim <- which(substr(txt[sidx], 1, 1) == ">")
        nseq <- length(seq.delim)

        if (nseq > 0) {
            # usually the last sequence is broken...
            if (seq.delim[1] > 1) {
                u <- 1:(seq.delim[1] - 1)
                str <- unlist(strsplit(tolower(c(s.tmp, txt[sidx][u])), "",
                                       fixed=TRUE))
            } else {
                str <- unlist(strsplit(tolower(s.tmp), "", fixed=TRUE))
            }

            if (length(str) != seq.len) {
                stop(paste("Sequences have to be of the same length.",
                           "Error at sequence number", n))
            }

            miss.seq.s <- str %in% miss

            # maybe this sequence has a non-missing value where it is missing
            # in the reference sequence
            j <- (miss.seq.ref - miss.seq.s) > 0

            if (any(j)) {
                seq.ref[j] <- str[j]
                miss.seq.ref[j] <- FALSE
            }

            # to be compared, both values must be non-missing
            j <- !(miss.seq.ref | miss.seq.s)
            poly[j] <- poly[j] | (str[j] != seq.ref[j])

            if (nseq > 1) {
                for (i in 2:nseq) {
                    u <- (seq.delim[i-1] + 1):(seq.delim[i] - 1)
                    str <- unlist(strsplit(tolower(txt[sidx][u]), "",
                                           fixed=TRUE))

                    if (length(str) != seq.len) {
                        stop(paste("Sequences have to be of the same length.",
                                   "Error at sequence number", n+i-1))
                    }

                    miss.seq.s <- str %in% miss

                    j <- (miss.seq.ref - miss.seq.s) > 0

                    if (any(j)) {
                        seq.ref[j] <- str[j]
                        miss.seq.ref[j] <- FALSE
                    }

                    j <- !(miss.seq.ref | miss.seq.s)
                    poly[j] <- poly[j] | (str[j] != seq.ref[j])
                }
            }

            if (seq.delim[nseq] < stot) {
                s.tmp <- txt[sidx][(seq.delim[nseq] + 1):stot]
            } else {
                s.tmp <- ""
            }

            n <- n + nseq

            cat(paste(n - 1, "sequences have been pre-processed\n"))
        } else {
            tmp <- s.tmp
            s.tmp <- character(length=length(tmp)+stot)
            s.tmp[1:length(tmp)] <- tmp
            s.tmp[(length(tmp)+1):(length(tmp)+stot)] <- txt[sidx]
        }
    }

    # by construction, the last sequence has not been pre-processed
    str <- unlist(strsplit(tolower(s.tmp), "", fixed=TRUE))

    if (length(str) != seq.len) {
        stop(paste("Sequences have to be of the same length. Error at ",
                   "sequence number", n))
    }

    miss.seq.s <- str %in% miss

    j <- (miss.seq.ref - miss.seq.s) > 0

    if (any(j)) {
        seq.ref[j] <- str[j]
        miss.seq.ref[j] <- FALSE
    }

    j <- !(miss.seq.ref | miss.seq.s)
    poly[j] <- poly[j] | (str[j] != seq.ref[j])

    if (n < 2) {
        stop("Found only a sequence. Not enough data.")
    }

    m <- sum(poly)

    cat(paste("Found", n, "sequences:", m, "SNPs out of", seq.len,
              "total sites\nProcessing data...\n"))

    # ? A C G T U R Y K M  S  W  B  D  H  V  N  X  -  E  F  I  L  O  P  Q  Z  *
    # 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
    # -------------------------------------------------------------------------
    #  J
    # 28
    enc <- 0:28
    names(enc) <- c("NA", "a", "c", "g", "t", "u", "r", "y", "k", "m", "s",
                    "w", "b", "d", "h", "v", "n", "x", "-", "e", "f", "i", "l",
                    "o", "p", "q", "z", "*", "j")

    d <- matrix(0, nrow=n, ncol=m)
    uid <- character(length=n)

    # 'c0' counts the current unit index
    # 'c1' and 'c2' are two counter variables, used for 'broken' sequences
    c0 <- 1
    c1 <- 0
    c2 <- 0

    # go back at the beginning of the file and start again
    a <- seek(con=fid, where=0)

    # skip blank lines
    done <- FALSE
    while (!done) {
        txt <- Trim(readLines(con=fid, n=1))

        if (!(length(txt) > 0 && identical(txt, ""))) {
            done <- TRUE
        }
    }

    # we already checked that the first non-empty element is the id
    uid[1] <- Trim(substr(txt, 2, nchar(txt)))

    done <- FALSE
    while (!done) {
        txt <- Trim(readLines(con=fid, n=1))

        if (length(txt) > 0 && !identical(substr(txt, 1, 1), ">")) {
            nc <- nchar(txt)

            if (nc > 0) {
                str <- strsplit(tolower(txt),
                                split="")[[1]][poly[(c1+1):(c1+nc)]]
                str[str %in% miss] <- "NA"

                z <- length(str)

                if (z > 0) {
                    d[c0, (c2+1):(c2+z)] <- enc[str]
                }

                c1 <- c1 + nc
                c2 <- c2 + z
            }
        } else {
            done <- TRUE
        }
    }

    if (identical(substr(txt, 1, 1), ">")) {
        uid[c0 <- c0 + 1] <- Trim(substr(txt, 2, nchar(txt)))
    }

    c1 <- 0
    c2 <- 0
    while (length(txt <- Trim(readLines(con=fid, n=lines.to.read))) > 0) {
        sidx <- nchar(txt) > 0
        stot <- sum(sidx)
        seq.delim <- which(substr(txt[sidx], 1, 1) == ">")
        nseq <- length(seq.delim)

        if (nseq > 0) {
            # usually the last sequence is broken...
            if (seq.delim[1] > 1) {
                # save sequence -> increase counter
                u <- 1:(seq.delim[1] - 1)

                nc <- sum(nchar(txt[sidx][u]))

                str <- unlist(strsplit(tolower(txt[sidx][u]), "",
                                       fixed=TRUE))[poly[(c1+1):(c1+nc)]]
                str[str %in% miss] <- "NA"

                z <- length(str)

                d[c0, (c2+1):(c2+z)] <- enc[str]

                c0 <- c0 + 1
            } else {
                c0 <- c0 + 1
            }

            c1 <- 0
            c2 <- 0

            if (nseq > 1) {
                for (i in 2:nseq) {
                    uid[c0] <- Trim(substr(txt[sidx][seq.delim[i-1]], 2,
                                           nchar(txt[sidx][seq.delim[i-1]])))

                    u <- (seq.delim[i-1] + 1):(seq.delim[i] - 1)

                    str <- unlist(strsplit(tolower(txt[sidx][u]), "",
                                           fixed=TRUE))[poly]
                    str[str %in% miss] <- "NA"

                    d[c0, ] <- enc[str]

                    c0 <- c0 + 1
                }
            }

            w <- seq.delim[nseq]
            uid[c0] <- Trim(substr(txt[sidx][w], 2 , nchar(txt[sidx][w])))

            if (w < stot) {
                u <- (w + 1):stot

                nc <- sum(nchar(txt[sidx][u]))

                str <- unlist(strsplit(tolower(txt[sidx][u]), "",
                                       fixed=TRUE))[poly[(c1+1):(c1+nc)]]
                str[str %in% miss] <- "NA"

                z <- length(str)

                d[c0, (c2+1):(c2+z)] <- enc[str]

                c1 <- c1 + nc
                c2 <- c2 + z
            }

            cat(paste(c0 - 1, "sequences have been processed\n"))
        } else {
            nc <- sum(nchar(txt[sidx]))

            str <- unlist(strsplit(tolower(txt[sidx]), "",
                                   fixed=TRUE))[poly[(c1+1):(c1+nc)]]
            str[str %in% miss] <- "NA"

            z <- length(str)

            d[c0, (c2+1):(c2+z)] <- enc[str]

            c1 <- c1 + nc
            c2 <- c2 + z
        }
    }

    cat(paste("All", c0, "sequences have been processed\n"))

    bd <- Data2Binary(d, enc)

    return(list(data=bd$binary.data, values=bd$values, keys=bd$keys,
                uid=uid, poly=poly))
}

###############################################################################
#
# ReadCsvData
#
# Description:
#
# Read data in CSV format
#
# Arguments:
#
#    input.file: path to data file
#        header: logical value indicating whether the file contains the names
#                of the variables as its first line
#           ids: logical value indicating whether the file contains the ids of
#                the units as its first column. If ids=FALSE, they will be
#                labeled by the numbers 1 to n
#          miss: character vector with the values to be considered missing. Use
#                NULL if no missing value is expected/considered. Use the empty
#                string "" if an empty space should be considered missing
# lines.to.read: The (maximal) number of lines to read. Negative values
#                indicate that one should read up to the end of input on the
#                connection
#
# Value:
#
# A list containing the following variables:
#   data: original data matrix encoded as a logical (binary) matrix
# values: character vector with unique values per column
#   keys: numeric vector with indices of each values
#    uid: n-by-1 character vector with the sample units' ids
#   poly: M-by-1 boolean vector indicating where the polymorphic variables
#         where found in the original dataset. Note that it doesn't refer to
#         the binary data matrix stored in variable 'data'!
#
###############################################################################
ReadCsvData <- function(input.file, header, ids, miss, lines.to.read) {
    cat("Loading CSV data...\n")

    fid <- file(input.file, "rt")
    on.exit(close(fid))

    # The main difference from "readFASTAdata" function is that we don't know
    # what kind of data we are going to read. This means that we could have
    # characters, numbers or both. To us, it doesn't matter. Of course, we
    # cannot build a standard conversion table

    # we read the file twice
    # the first time we check the total number of units and where are the
    # polymorphic columns. The second time we store the data

    # Note: This function has been written for reading a huge dataset

    n <- 0
    keep.reading <- FALSE

    # read the first unit. we skip the header if present
    # header==TRUE => skip = 1; header==FALSE => skip = 0
    skip <- as.integer(header)

    txt <- scan(file=fid, what=character(0), nlines=1, skip=skip,
                na.strings=miss, sep=",", quote="\"", dec=".", comment.char="",
                fill=FALSE, strip.white=TRUE, quiet=TRUE,
                blank.lines.skip=TRUE, multi.line=FALSE)

    if (length(txt) > 0) {
        n <- 1
        len <- length(txt)
        poly <- logical(length=len)
        unit.ref <- tolower(txt)
        what <- rep.int(list(character(0)), len)

        if (!ids) {
            M <- len
            idx.set <- 1:len
        } else {
            M <- len - 1
            idx.set <- 2:len
            unit.ref[1] <- NA
        }

        keep.reading <- TRUE
    }

    miss.unit.ref <- is.na(unit.ref)
    al <- unique(unit.ref[!miss.unit.ref])

    if (length(al) > 0) {
        W <- floor(length(al) / 1e+03) * 1e+03
        alphabet <- character(length=max(c(W, 1e+03)))
        alphabet[1:length(al)] <- al
        al.idx <- length(al)
    } else {
        # only missing data? does it make sense? just to be sure, handle this
        # case anyway
        alphabet <- character(length=1e+03)
        al.idx <- 0
    }

    # we use the two options 'fill=FALSE' and 'multi.line=FALSE' to check if
    # the units have the same number of variables
    while (keep.reading) {
        txt <- scan(file=fid, what=what, nmax=lines.to.read, na.strings=miss,
                    sep=",", quote="\"", dec=".", comment.char="", fill=FALSE,
                    strip.white=TRUE, quiet=TRUE, blank.lines.skip=TRUE,
                    multi.line=FALSE)

        if (!all(sapply(txt, length) == 0)) {
            n <- n + length(txt[[1]])

            for (i in idx.set) {
                r <- tolower(txt[[i]])
                miss.r <- is.na(r)

                # NOTE: miss.unit.ref represents all the columns (without ids)
                #       while miss.r is just one column
                if (miss.unit.ref[i] && any(!miss.r)) {
                    unit.ref[i] <- r[!miss.r][1]
                    miss.unit.ref[i] <- FALSE
                }

                if (!miss.unit.ref[i]) {
                    poly[i] <- poly[i] | any(r[!miss.r] != unit.ref[i])
                }

                al <- unique(r[!miss.r])
                al.new <- !(al %in% alphabet)

                if (any(al.new)) {
                    w <- al.idx + sum(al.new)

                    if (w > length(alphabet)) {
                        W <- floor(w / length(alphabet)) * 1e+03
                        tmp <- alphabet
                        alphabet <- character(length=length(tmp)+W)
                        alphabet[1:length(tmp)] <- tmp
                    }

                    alphabet[(al.idx+1):w] <- al[al.new]
                    al.idx <- w
                }
            }

            cat(paste(n, "units have been pre-processed\n"))
        } else {
            keep.reading <- FALSE
        }
    }

    if (n < 2) {
        stop("Found only a unit. Not enough data.")
    }

    alphabet <- alphabet[alphabet != ""]

    # try to convert the alphabet to numeric type
    alphabet.num <- suppressWarnings(as.numeric(alphabet))

    if (any(is.na(alphabet.num))) {
        alphabet <- sort(alphabet)
    } else {
        alphabet <- alphabet[order(alphabet.num)]
    }

    m <- sum(poly)

    cat(paste("Found a total of", n, "units:", m, "polymorphic variables out",
              "of", M, "total variables\nProcessing data...\n"))

    enc <- 0:length(alphabet)
    names(enc) <- c("?", alphabet)

    d <- matrix(0, nrow=n, ncol=m)
    uid <- character(length=n)

    # go back at the beginning of the file and start again
    a <- seek(con=fid, where=0)

    # skip the header (if any) and read the first sequence
    txt <- scan(file=fid, what=character(0), nlines=1, skip=skip,
                na.strings=miss, sep=",", quote="\"", dec=".", comment.char="",
                fill=FALSE, strip.white=TRUE, quiet=TRUE,
                blank.lines.skip=TRUE, multi.line=FALSE)

    c0 <- 1

    if (!ids) {
        uid[c0] <- as.character(c0)
    } else {
        uid[c0] <- txt[1]
    }

    d[c0, ] <- match(txt[poly], alphabet, nomatch=0)

    # we now know that each line has the same number of variables
    while(length(txt <- scan(file=fid, what=character(0), nlines=lines.to.read,
                             na.strings=miss, sep=",", quote="\"", dec=".",
                             comment.char="", fill=FALSE, strip.white=TRUE,
                             quiet=TRUE, blank.lines.skip=TRUE)) > 0) {
        # txt is now a vector of length lines.read * len
        # length(txt) / len is then the total number of lines read
        # We are computing this because the last iteration could read less than
        # 'lines.to.read' lines
        lr <- as.integer(length(txt) / len)

        # values from the same column are found every 'len' positions
        idx.units <- (c0+1):(c0+lr)

        if (!ids) {
            uid[idx.units] <- as.character(idx.units)
        } else {
            uid[idx.units] <- txt[1 + ((0:(lr-1))*len)]
        }

        idx <- which(matrix(poly, nrow=length(poly), ncol=lr))

        d[idx.units, ] <- matrix(match(txt[idx], alphabet, nomatch=0),
                                 nrow=lr, ncol=m, byrow=TRUE)

        c0 <- c0 + lr

        cat(paste(c0, "units have been processed\n"))
    }

    cat(paste("All", c0, "units have been processed\n"))

    bd <- Data2Binary(d, enc)

    if (ids) {
        poly <- poly[-1]
    }

    return(list(data=bd$binary.data, values=bd$values, keys=bd$keys,
                uid=uid, poly=poly))
}

###############################################################################
#
# Data2Binary
#
# Description:
#
# Convert raw data to binary format
#
# Arguments:
#
#   d: n-by-m integer matrix. Only integer numbers are allowed. A value of 0 is
#      considered missing
# enc: encoding table for integer values stored in 'd'
#
# Value:
#
# binary.data: raw.data as a binary matrix
#      values: integer vector with unique values per column
#        keys: numeric vector with indices of each values
#
###############################################################################
Data2Binary <- function(d, enc) {
    cat("Converting data to binary data...\n")

    n <- nrow(d)
    m <- sum(apply(d, 2, function(x) sum(unique(x) > 0)))

    binary.data <- matrix(FALSE, nrow=n, ncol=m)

    values <- character(length=m)
    keys <- integer(length=m)

    counter <- 1
    for (j in 1:ncol(d)) {
        v.col <- d[, j]

        vl <- unique(v.col)
        vl <- sort(vl[vl > 0])

        nv <- length(vl)

        m1 <- matrix(v.col, nrow=n, ncol=nv)
        m2 <- matrix(vl, nrow=n, ncol=nv, byrow=TRUE)
        v.tmp <- (m1 == m2)

        idx <- counter:(counter + nv - 1)

        values[idx] <- names(enc[vl+1])
        keys[idx] <- j
        binary.data[, idx] <- v.tmp

        counter <- counter + nv
    }

    return(list(binary.data=binary.data, values=values, keys=keys))
}

###############################################################################
#
# Trim
#
# Description:
#
# Strip leading and trailing blanks
#
# Arguments:
#
# x: character vector
#
# Value:
#
# A vector of the same length of x, where each element does not have leading or
# trailing blanks
#
###############################################################################
Trim <- function(x) {
    return(sub("\\s+$", "", sub("^\\s+", "", x, perl=TRUE), perl=TRUE))
}
