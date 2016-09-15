# K-Pax2 - Bayesian Cluster Analysis of Categorical Data

### Note: development of this package is now discontinued.

K-Pax2 has been superseded by
[Kpax3.jl](https://github.com/albertopessia/Kpax3.jl), a package for the Julia
programming language. Kpax3.jl provides the same functionalities of
K-Pax2 but in a faster and better way. It also add a MCMC algorithm for a
proper Bayesian analysis. Check it out!

## About

K-Pax2 is a R package, written with the purpose of clustering (big) datasets
of categorical statistical variables. Main application of K-Pax2 is with
genetic datasets, such as dna/protein multiple sequence alignments. Being a
general method, it can be easily applied to any kind of categorical dataset.
K-Pax2 output consists of a classification of both the rows (statistical units)
and columns (statistical variables) of the provided data matrix

## Reference

To know more about the method and how to properly setup the parameters, please
read the following paper:

Pessia A., Grad Y., Cobey S., Puranen J., Corander J. (2015). _K-Pax2: Bayesian
identification of cluster-defining amino acid positions in large sequence
datasets_. MGen 1(1). doi:
[10.1099/mgen.0.000025](https://dx.doi.org/10.1099/mgen.0.000025)

Visit also the [BSG group website](http://www.helsinki.fi/bsg/) for other
interesting projects.

## Installation

### devtools

The easiest way to install K-Pax2 is using
[devtools](https://github.com/hadley/devtools):

    devtools::install_github("albertopessia/kpax2/kpax2")
    library(kpax2)

### Standard approach

Download the latest version of K-Pax2 from the
[releases page](https://github.com/albertopessia/kpax2/releases) and extract the
content of the compressed file into a folder of your choice. Let
"PathToPackageFile" denote the complete path to the package file (having
extension .tar.gz).

Decide if the package should be made available only to the current user or
to every user on the same machine. In the latter case, run all the following
commands with administrative privileges.

Note for Windows users: within R, file paths should be written using the slash
'/' character instead of the usual backslash '\' character. If you want to use
the backslash, you have to escape it. Examples of valid file paths:

    "C:/Users/CurrentUser/Documents/kpax2/kpax2_x.x.x.tar.gz"
    "C:\\Users\\CurrentUser\\Documents\\kpax2\\kpax2_x.x.x.tar.gz"

#### Default system library

If you want to install the package into the default system library, run R and
type the command

    install.packages("PathToPackageFile", repos=NULL, type="source")

An alternative approach is to open the terminal and issue the command

    R CMD INSTALL "PathToPackageFile"

If it fails to find R, then it means that it is not in the system path. Ask
your machine administrator for help.

You can now load the package in R using

    library(kpax2)

#### Custom library

If you want to install the package into a custom library, different from the
default system one, choose its location and denote with "PathToLibrary" its
complete path. Open the terminal and issue the command

    R CMD INSTALL --library="PathToLibrary" "PathToPackageFile"

You can now load the package in R using

    library(kpax2, lib.loc="PathToLibrary")

## Usage

For a tutorial on how to properly use this package, follow the instructions
written in the file [kpax2_tutorial.R](tutorial/kpax2_tutorial.R).

## License

See [LICENSE.md](LICENSE.md)
