\name{LoadData}
\alias{LoadData}

\title{K-Pax2 Input Tool Facility}
\description{Read data in FASTA or CSV format and convert it to K-Pax2 format}

\usage{LoadData(input.file, output.file=NULL, filetype="fasta", header=FALSE,
                ids=FALSE, miss=NULL, lines.to.read=-1)}

\arguments{
    \item{input.file}{path to input data file}
    \item{output.file}{path to a Rdata output file. Use NULL if output should not be saved on disk}
    \item{filetype}{string representing the data format to read. Use either "csv" or "fasta"}
    \item{header}{logical value. Used only if filetype="csv". This indicates whether the file contains the names of the variables as its first line}
    \item{ids}{logical value. Used only if filetype="csv". This indicates whether the file contains the ids of the units as its first column. If ids=FALSE, they will be labeled by the numbers 1 to n}
    \item{miss}{character vector with the values to be considered missing. Use NULL if no missing value is expected/considered}
    \item{lines.to.read}{The (maximal) number of lines to read at each iteration. Should be used when reading a very large dataset. Negative values indicate that one should read up to the end of input on the connection in one single step.}
}

\details{
For argument miss, only single characters (length 1) are allowed when using FASTA data files.

For FASTA data, it is assumed that the only possible characters are the standard ones (see for example http://en.wikipedia.org/wiki/FASTA_format ). Unpredictable things could happen if unknown non-standard characters are found within the FASTA file. Be sure the file is a standard FASTA file. If not, consider to convert the file into a CSV file.

Thymidine (T) and Uracil (U) are given different codes (because T and U characters are used for different amino acids) and for this reason, handled as different. If you find both in a DNA/RNA dataset, consider to convert one letter to the other.
}

\value{
A list (and a Rdata file if output.file is not NULL) containing the following variables
\item{data}{original data matrix encoded as a logical (binary) matrix}
\item{values}{character vector with unique values per column}
\item{keys}{numeric vector with indices of each values}
\item{uid}{n-by-1 character vector with the sample units' ids}
\item{poly}{M-by-1 boolean vector indicating where the polymorphic variables where found in the original dataset. Note that it doesn't refer to the binary data matrix stored in variable 'data'!}
}

\examples{
# Typical values for 'miss' are (considering gaps as missing values is problem specific):
# - DNA data
miss <- c("-", "?", "*", "#", "b", "d", "h", "k", "m", "n", "r", "s", "v", "w", "x", "y")
# - Protein data
miss <- c("-", "?", "*", "#", "b", "j", "x", "z")
}
