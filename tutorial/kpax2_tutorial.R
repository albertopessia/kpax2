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
# INTRODUCTION
#
###############################################################################
#
# Welcome to K-Pax2 tutorial!
#
# In this tutorial you will learn how to:
# - load a multiple sequence alignment (MSA) or a more generic dataset
# - setup the parameters
# - run the algorithm
# - save the results
#
###############################################################################

###############################################################################
#
# NOTE FOR WINDOWS USERS
#
###############################################################################
#
# File paths should be written using the slash '/' character instead of the
# usual backslash '\' character. If you want to use the backslash, then you
# have to escape it first (i.e. use '\\'). Examples of valid file paths:
#
# "C:/Program Files/R/R-3.0.0/bin/R.exe"
# "C:\\Program Files\\R\\R-3.0.0\\bin\\R.exe"
#
###############################################################################

###############################################################################
#
# LOAD THE PACKAGE
#
###############################################################################
#
# First of all, we need to load the package "kpax2". I will assume the package
# it's already installed somewhere into your system. If it is not, please
# follow the installation instructions written in the README file.
#
# To load the package, open a R session and run the command
library(kpax2)

# If you installed the package into a non-standard library folder, then use:
# library(kpax2, lib.loc="path/to/nonstandard/library/folder")

###############################################################################
#
# PRELIMINARY SETUP
#
###############################################################################
#
# The tutorial has been designed in such a way that you should just copy/paste
# the commands into your R environment. To make it possible, we need to change
# our working directory to a common one.
#
# Extract folder kpax2_tutorial from the zip file wherever your prefer (be sure
# to have write permissions at that location).
#
# Now change the following command according to its path
setwd("path_to_folder_kpax2_tutorial")

###############################################################################
#
# IMPORT DATA
#
###############################################################################
#
# The data we will use in this tutorial consist of 200 aligned protein
# sequences of the HIV-1 "Env" gene
#
# Data has been randomly selected from the good quality sequences (i.e. no
# frameshifts and stop codons within the sequence) of the HIV Databases 2013
# Filtered Web Alignments
# http://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html
#
# HIV has been chosen for its high mutation and recombination rate (as can be
# easily seen from the MSA). We will check how K-Pax2 will handle such a
# peculiar dataset
#
# First of all, we need to convert the data into the K-Pax2 format. To do so,
# we need to specify which characters are to be considered missing data.
# We are going to analyze protein sequences, so let's just use the usual
# missing characters for a protein FASTA format.
miss.char <- c("?", "*", "#", "-", "b", "x", "z")

# Notice that we also put the gap symbol among the missing data. Sometimes we
# would like to consider it a valid character. It's your choice.

# It's a good idea to save the converted dataset, so that we won't need to
# convert it again in a later analysis. This will also save some time when
# analyzing, in parallel in a cluster environment, a huge dataset.

# IMPORTANT: the default name for the data variable, saved in the output file,
#            is "full.data". Keep this in mind when writing automated scripts.

# Choose the output file name (and path maybe)
data.output.file <- "HIV1_dataset.Rdata"

# Finally load the data
full.data <- LoadData(input.file="HIV1_dataset.fasta",
                      output.file=data.output.file,
                      miss=miss.char)

# Check its structure just to get familiar with it
str(full.data)

# In this tutorial we are analyzing protein sequences, but it's good to know
# that it's possible to load any kind of categorial data (as csv file format).
# To know more about loading data, type
?LoadData

###############################################################################
#
# SETUP PARAMETERS
#
###############################################################################
#
# K-Pax2 requires the definition of 3 main input parameters:
#
# - an initial starting partition
# - a vector of prior probabilities for the column classification
# - an array of parameters for the beta distributions as used by the model
#
# plus some variables used to customize a K-Pax2 run

# Initial partition
# The dataset is not very big, so we can start our search by placing each unit
# into its own cluster. This should make it easier for the algorithm to find
# the global optimum (or a close solution to it).
#
# When clustering thousands of units, this approach is not optimal from a
# memory point of view. In this case, we could choose a random partition or
# start already from a reasonable one. For example, we could cut a hierarchical
# cluster tree.
#
# It is also possible to use an external file to define the initial partition.
# The file must have only one column while the number of rows must be equal to
# the sample size. Each row 'i' represents the cluster of unit 'i'. Check file
# 'HIV1_best_partition.csv' to see how to format a partition file.
#
# IMPORTANT: For the algorithm is easier to merge clusters rather than
#            splitting them. It is always a good idea to start from a partition
#            that has many more clusters than the optimal solution.
#            This requires some initial guessing, but we can refine it later
#            after few runs.
#initial.partition <- "best_partition.csv"
initial.partition <- 1:200

# Column status prior probabilities
# This is a non-negative numeric vector of length 3. If it does not sum to 1,
# it is going to be automatically normalized.
# Probabilities contained in this variable represent an initial guess for the
# proportion of sites possessing the following statuses:
# - uninformative (noise)
# - informative but not characteristic for any cluster (weak signal)
# - informative and characteristic for a cluster (strong signal)
#
# These default values have proven to work well for many different datasets.
status.prior.probs <- c(0.6, 0.35, 0.05)

# Beta distribution hyper-parameters
# These are a bit more difficult to explain and to properly setup. You can do
# it on your own, of course, but it requires some technical knowledge about the
# R package code and a good understanding of the model.
# There is a function that will setup a proper array for you, using the default
# approach explained in the paper. This also has  proven to work well for many
# different datasets.
beta.hyper.par <- InitHyperParameters(n=200,
                                      n1s=apply(full.data$data, 2, sum),
                                      r=(log(0.001) / log(0.95)))

# To know more about these hyperparameters and what the value of r is, please
# read the paper (and the help file)
?InitHyperParameters

# When analyzing a dataset made of thousands of units, memory allocation
# becomes an extremely important factor. We could try to save some memory by
# fixing in advance a maximum number of clusters we are (almost) sure we will
# never reach.
# Note that if the optimal partition has more than this maximum number of
# clusters, the algorithm will automatically allocate more memory on the fly.
# Memory allocation at running time is a very expensive operation, so don't
# define this variable to a small number
# If memory allocation is not an issue, we can just set this variable to n
max.K <- 200

# Print a status message every 50 iterations
iter.print <- 50

# Save a backup partition every "iter.print" iterations, in case something bad
# happen at running time (crash, hang or manual stop)
backup.file <- "backup_partition.csv"

###############################################################################
#
# RUN K-PAX2
#
###############################################################################
#
# We are now finally ready to run K-Pax2!
#
# There are 3 main different approaches:
# a) Run it sequentially (e.g. using a for loop)
# b) Run it in parallel (e.g. using multiple R sessions at the same time)
# c) A mixture of (a) and (b)
#
# In any case, we will pick the best solution among all the runs.
#
# We are not going to do any of them, as it is not the purpose of this
# tutorial. We will just run it once, to show how it works, and then we will
# load a previous (best) solution from an external file.

# Beware, this is going to take a (long) while. Maybe you want to skip this and
# move to the next section. Check the help file for more details
?Kpax2

res <- Kpax2(d=full.data$data,
             init.part=initial.partition,
             g.prob=status.prior.probs,
             hyper.par=beta.hyper.par,
             max.clust=max.K,
             t.iter=iter.print,
             bak.file=backup.file)

# Check the structure of K-Pax2 results
# res$partition is the optimal partition found by the algorithm
# res$max.logpp is its corresponding posterior probability
# res$idx.max is a matrix with the columns classification
str(res)

###############################################################################
#
# SAVE RESULTS
#
###############################################################################
#
# To save the results, first choose a path to a file (multiple files will be
# created using this same root)
res.output.file <- "kpax2_results"

# and then (obviously you can't do this at this moment if you skipped the
# previous Kpax2 command. Don't worry, skip also this and keep reading)
SaveResults(full.data=full.data,
            result=res,
            output.file=res.output.file,
            what=4)

# You can choose how much information you want to get by changing parameter
# 'what', starting from 1 (basic information) to 4 (most information)
?SaveResults

# Output description:
# filename_summary.txt        = Text file with the log posterior probability.
#                               This can be used to compare different parallel
#                               runs.
# filename_partition.csv      = Best partition found by the algorithm.
# filename_attributes.csv     = Complete columns classification matrix.
# filename_characteristic.csv = Characteristic sites and their corresponding
#                               values.
# filename_info.txt           = Original dataset with highlighted clusters and
#                               their corresponding characteristic values.

# Now, suppose we are not satisfied with the current best partition and we
# think that units should be moved around, or a cluster should be split or two
# clusters merged. This could happen when metadata is available, so that we
# can immediately see if something is not right.
#
# After manually editing the partition file, we can load it again and see
# if we arrived to a better solution or not. In this tutorial, we are going to
# load the best partition known for this dataset.
my.partition <- "HIV1_best_partition.csv"

# No need to run the algorithm again. There is a function which has been wrote
# just to fulfill this task.
?GetMaxLogPP

my.result <- GetMaxLogPP(d=full.data$data,
                         init.part=my.partition,
                         g.prob=status.prior.probs,
                         hyper.par=beta.hyper.par)

# If you happen to find a better solution, please let me know!

# You can explore results associated with the best solution to understand
# better K-Pax2 output. This command will overwrite the previous (probably
# worse) results.
SaveResults(full.data=full.data,
            result=my.result,
            output.file=res.output.file,
            what=4)

# The best partition is made of 9 clusters. By looking at the metadata, we can
# easily understand why they have been grouped this way:
# 1) HIV-1 Subtype B
# 2) HIV-1 Subtype C
# 3) HIV-1 Subtype A1
# 4) HIV-1 Subtype CRF_AE
# 5) HIV-1 Subtype G
# 6) HIV-1 Subtype D
# 7) HIV-1 Subtype BF
# 8) HIV-1 Recombinant / Unclear Subtype
# 9) HIV-1 Group O and Group N
#
# The added value of K-Pax2, when comparing it to existing clustering methods,
# is that we now also know which amino acids (check
# 'kpax2_results_characteristics.csv' file) are responsible for the difference
# between these subtypes