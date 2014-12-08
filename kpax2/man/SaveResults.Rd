\name{SaveResults}
\alias{SaveResults}

\title{K-Pax2 Output Tool Facility}
\description{Write the result of K-Pax2 into several text files}

\usage{SaveResults(full.data, result, output.file, what=1)}

\arguments{
    \item{full.data}{a list object, output of LoadData function}
    \item{result}{a list object, output of 'Kpax2' or 'GetMaxLogPP' functions}
    \item{output.file}{character string representing a path to the output file}
    \item{what}{integer scalar representing what to print}
}

\details{
Several output files are produced by SaveResults (depending on the value of variable 'what'). No file extension is required for argument output.file and the file name should only represent a common base name.

Possible values for argument 'what' are the following. 1: basic summary and partition only. 2: basic summary, partition and column classification. 3: basic summary, partition, column classification and characteristic features. 4: basic summary, partition, column classification, characteristic features and original data.
}

\value{Text files containing K-Pax2 output}
