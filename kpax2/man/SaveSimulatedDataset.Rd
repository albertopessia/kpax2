\name{SaveSimulatedDataset}
\alias{SaveSimulatedDataset}

\title{K-Pax2 Output Tool Facility}
\description{Write the values returned by SimulateDataset function into several text files}

\usage{SaveSimulatedDataset(sim, output.file)}

\arguments{
    \item{sim}{a list object, output of SimulatedDataset function}
    \item{output.file}{character string representing a path to the output file}
}

\details{
Several output files are produced by SaveSimulatedDataset. No file extension is required for argument output.file and the file name should only represent a common base name.
}

\value{Text files representing the simulated dataset}
