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
pander(PCA.res[["pcs"]][1:10,])

## ---- out.width='70%', fig.width=10, fig.height=10-----------------------
DrawPCA(PCA.res, drawNames=TRUE,
 createWindow=FALSE, main = "PCA score plot for HumanSerum dataset",
   class = ClassHS, axes =c(1,2), type ="scores")

## ----Loadings------------------------------------------------------------
pander(PCA.res[["pcv"]][1:10,])

## ---- out.width='100%', fig.width=12, fig.height=8-----------------------
DrawPCA(PCA.res, drawNames=TRUE,
 createWindow=FALSE, main = "PCA loadings plot for HumanSerum dataset",
    axes = 1, type ="loadings", loadingstype="l")

## ----ClustMIC, out.width='100%', fig.width=12, fig.height=12-------------
data("HumanSerum")
ClustMIC.res = MBXUCL::ClustMIC(Intensities = HumanSerumSpectra, nClust = 4, Trcl = ClassHS, Dendr = TRUE)

pander(ClustMIC.res)

## ----binClustMIC, out.width='100%', fig.width=12, fig.height=12----------
Pos = t(GPL[, seq(4,dim(GPL)[2], 2)])
Trcl = c(rep(1,8), rep(2,8), rep(3,8))
binClustMIC.res = MBXUCL::binClustMIC(Positions = Pos, Distance = "Jaccard", nClust = 3, Trcl, Dendr = TRUE)

pander(binClustMIC.res)

## ----PLSDA, out.width='100%', fig.width=10, fig.height=10----------------
data("HumanSerum")
PLSDA.res = PLSDA(x = HumanSerumSpectra, y = ClassHS, nLV = NULL, drawRMSEP = TRUE)

perf.plsda = PLSDA.res[4:6]
pander(perf.plsda)

