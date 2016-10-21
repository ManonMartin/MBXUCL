## ----setup, include = FALSE----------------------------------------------
require(knitr)
require(pander)
require(MBXUCL)


## ----install1, tidy=TRUE,  eval=FALSE------------------------------------
#  
#  install.packages(file.path("path to MBXUCL-package", "MBXUCL-package/MBXUCL"), repos = NULL, type="source")
#  require(MBXUCL)
#  

## ----install2, tidy=TRUE,  eval=FALSE------------------------------------
#  require(devtools)
#  install_github("ManonMartin/MBXUCL")
#  require(MBXUCL)
#  

## ----inertia-------------------------------------------------------------
data("HumanSerum")
Inertia.res = MBXUCL::Inertia(x = HumanSerumSpectra, y = ClassHS, print = TRUE)

kable(Inertia.res[["Between_within"]])

kable(Inertia.res[["Per_group"]])

## ----pca-----------------------------------------------------------------
data("HumanSerum")
PCA.res = MBXUCL::SVDforPCA(HumanSerumSpectra, ncomp=4)

## ----eigPCA--------------------------------------------------------------
pander(PCA.res[["eigval"]][1:4])

## ----Scores--------------------------------------------------------------
pander(PCA.res[["scores"]][1:10,])

