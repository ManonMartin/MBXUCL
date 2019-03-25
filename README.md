[![Build Status](https://travis-ci.org/ManonMartin/MBXUCL.svg?branch=master)](https://travis-ci.org/ManonMartin/MBXUCL)

# MBXUCL
R package for the AHDIDA Team (ISBA, UCLouvain).

## Installation 
### Install the package from GitHub (recommended)



#### Solution 1: using package `devtools`
```R
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("manonmartin/mbxucl")
```

####  Solution 2: using package `remotes`
```R
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github("manonmartin/mbxucl")
```



### Install the package from a local repository:

####  Solution 3: install from a local copy 
```R

#A. Install necessary dependencies ===========
install.packages(c("clValid", "phyclust", "proxy", "pls", 
    "pander", "stats","ggplot2", "reshape2", "spls", "plyr", 
    "gridExtra", "clusterSim", "modeest"))


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ropls")

#B. Install MBXUCL ===========

# Instructions: 
# 1. Download the MBXUCL library from https://github.com/ManonMartin/MBXUCL
# 2. Unzip the folder
# 3. Change package_path to the folder location
package_path <- "/Users/manon/Documents/MBXUCL"
# 4. Install MBXUCL
install.packages(pkgs = package_path, repos = NULL, 
                 type = "source")
```

/!\\ package dependencies should be already available if installed 
from a local repository.

## Load the package in R

Once installed, the package can be loaded as usual in R with `library(MBXUCL)`.
