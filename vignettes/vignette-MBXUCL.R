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
#  install_github("ManonMartin/MBXUCL", dependencies = TRUE)
#  require(MBXUCL)
#  

## ------------------------------------------------------------------------
#---- Simulated spectra
DataSimul$x[1:10, 1:10]
dim(DataSimul$x)

LinePlot(DataSimul$x, rows = 1)

#---- Class
table(DataSimul$y)

## ------------------------------------------------------------------------
#---- Serum spectra
HumanSerumSpectra[1:10, 1:10]
dim(HumanSerumSpectra)

LinePlot(HumanSerumSpectra, rows = 1)

#---- Class
table(ClassHS)

## ----inertia-------------------------------------------------------------
Inertia.res = MBXUCL::Inertia(x = HumanSerumSpectra, y = ClassHS, print = TRUE)

kable(Inertia.res[["Between_within"]])

kable(Inertia.res[["Per_group"]])

## ----pca-----------------------------------------------------------------
PCA.res = MBXUCL::SVDforPCA(HumanSerumSpectra, ncomp=4)

## ----eigPCA--------------------------------------------------------------
pander(PCA.res[["eigval"]][1:4])

## ----Scores--------------------------------------------------------------
pander(PCA.res[["scores"]][1:10,])

## ----scoreplot, out.width='70%', fig.width=10, fig.height=10-------------
DrawScores(PCA.res, type.obj = "PCA", drawNames=TRUE,
 createWindow=FALSE, main = "PCA scores plot for HumanSerum dataset",
   color = ClassHS, pch = ClassHS,  axes =c(1,2), drawEllipses = TRUE)

## ----Loadings------------------------------------------------------------
pander(PCA.res[["loadings"]][1:10,])

## ---- out.width='100%', fig.width=12, fig.height=8-----------------------
DrawLoadings(PCA.res, type.obj = "PCA",
 createWindow=FALSE, main = "PCA loadings plot for HumanSerum dataset",
    axes = c(1:2),  loadingstype="s")

## ----ClustMIC, out.width='100%', fig.width=12, fig.height=12-------------
ClustMIC.res = MBXUCL::ClustMIC(Intensities = HumanSerumSpectra, nClust = 4, Trcl = ClassHS, Dendr = TRUE)

pander(ClustMIC.res)

## ----binClustMIC, out.width='100%', fig.width=12, fig.height=12----------
Pos = t(GPL[, seq(4,dim(GPL)[2], 2)])
Trcl = c(rep(1,8), rep(2,8), rep(3,8))
binClustMIC.res = MBXUCL::binClustMIC(Positions = Pos, Distance = "Jaccard", nClust = 3, Trcl, Dendr = TRUE)

pander(binClustMIC.res)

## ----PLSDA, out.width='100%', fig.width=10, fig.height=10----------------
Class = ClassHS
PLSDA.res = PLSDA(x = HumanSerumSpectra, y = Class, nLV = NULL, drawRMSEP = TRUE)

perf.plsda = PLSDA.res[4:6]
pander(perf.plsda)

## ---- out.width='100%', fig.width=12, fig.height=12----------------------
DrawScores(PLSDA.res, type.obj = "PLSDA", drawNames=TRUE,
 createWindow=FALSE, main = "PLSDA scores plot for HumanSerum dataset",
   color = ClassHS, pch = ClassHS,  axes =c(1,2), drawEllipses = TRUE)

DrawLoadings(PLSDA.res, type.obj = "PLSDA",
 createWindow=FALSE, main = "PLSDA loadings plot for HumanSerum dataset",
    axes = c(1:2),  loadingstype="s")


## ----OPLSDA--------------------------------------------------------------
x = DataSimul[["x"]]
y = DataSimul[["y"]]
m = dim(x)[2]
no=3
nb = 15
oplsda.res = OPLSDA(x=x, y=y, impT = FALSE,impG = FALSE, no=no, nb = nb)

## ---- out.width='70%',fig.width=8, fig.height=8--------------------------
COL=rep("gray93",no)
mp=barplot(oplsda.res$CV,
           axes=F, axisnames=F, border=1,col=COL)
axis(1,at=mp, labels=c(1:no))
axis(2)
title(main="OPLS: Choice of the n[orthog. Components]", xlab="Orthogonal OPLS-DA Components",   ylab="||Wortho|| / ||p||", cex.main = 0.8)


## ---- out.width='70%', fig.width=10, fig.height=10, eval=FALSE-----------
#  # not evaluated
#  cvOPLSDA.res = cvOPLSDA(x = x, y = y, k_fold = 10, NumOrtho = 10)
#  plot(1:length(cvOPLSDA.res$RMSECV), cvOPLSDA.res$RMSECV, main = "RMSECV", lty = 19)
#  

## ---- out.width='100%', fig.width=12, fig.height=12----------------------
Class = y
DrawScores(oplsda.res, type.obj = "OPLSDA", drawNames=TRUE,
 createWindow=FALSE, main = "OPLSDA scores plot for HumanSerum dataset",
   color = Class, pch = Class,  axes =c(1,2), drawEllipses = TRUE)

DrawLoadings(oplsda.res, type.obj = "OPLSDA",
 createWindow=FALSE, main = "OPLSDA loadings plot for HumanSerum dataset",
    axes = c(1:3),  loadingstype="s")


## ---- out.width='80%', fig.width=12, fig.height=12-----------------------

# S-plot
###########################
## correlation and covariance matrices
# covariance
s = as.matrix(x, ncol = ncol(x))
p1 = c()
for (i in 1:ncol(s)) {
  scov = cov(s[, i], oplsda.res$Tp)
  p1 = matrix(c(p1, scov), ncol = 1) # covariance x-T
}

# correlation
pcorr1 = c()
Tno=as.matrix(oplsda.res$Tp, ncol=1)
for (i in 1:nrow(p1)) {
  den = apply(Tno, 2, sd, na.rm = TRUE) * sd(s[, i])
  corr1 = p1[i, ]/den
  pcorr1 = matrix(c(pcorr1, corr1), ncol = 1) # correlation
}

# plot

plot(p1, pcorr1, xlab = "p(cov)[1]", ylab = "p(corr)[1]",
     main = "S-plot (OPLS-DA)", ylim=c(min(pcorr1, na.rm = T)*1.1,
                                       max(pcorr1, na.rm = T)*1.1), xlim=c(min(p1, na.rm = T)*1.1, 
                                                                           max(p1, na.rm = T)*1.1))

sel = p1*pcorr1
sel = order(sel, decreasing = TRUE)[1:nb]
text(p1[sel], pcorr1[sel], labels = colnames(s)[sel], cex = 0.7, pos = 1)
abline(v=0,lty = 2)
abline(h=0,lty = 2)




## ---- out.width='100%', out.width='70%', fig.width=12, fig.height=5------

delta=mean(sort(abs(oplsda.res$b))[m-nb+c(0,1)])
xax = as.numeric(names(oplsda.res$b))

par(mar=c(4,2,2,1))
plot(oplsda.res$b,type="l",xaxt = "n", yaxt = "n", main="OPLS: Vector of descriptors' rank",xlab="ppm" )
abline(h=0)
abline(h=delta*c(-1,1),lty=2)
axis(side=2,cex.axis=0.7)
axis(side=1,at=c(1,seq(50,m,50)),
     labels=xax[c(1,seq(50,m,50))],cex.axis=0.7)


## ------------------------------------------------------------------------
xtrain = x[-c(1:5),]
ytrain = y[-c(1:5)]
xnew = x[c(1:5),]
oplsda.res = OPLSDA(x=xtrain, y=ytrain,
                    impT = FALSE,impG = FALSE, no=2, nb = 15, out.path = ".")
OPLSDA_pred(ropls = oplsda.res, x.new = xnew)



## ---- fig.width=8, fig.height=8------------------------------------------
x <- c(1:20)
y <- sample(20)
colors <- c(rep(1,10), rep(2,10))
points_pch <-  c(rep(17,10), rep(18,10))
labels <- as.character(c(1:20))

# decide on the color or pch values
# while adding a personalized legend
legend_color = c(rep("A", 10), rep("B", 10))
legend_pch = c(rep("C", 10), rep("D", 10))
ScatterPlot(x = x, y = y, points_labs=labels, createWindow=FALSE,
            main = 'Scatter plot', color = colors, pch = points_pch,
            xlab = "x axis", ylab = "y axis",
            legend_pch = legend_pch,legend_color=legend_color)

## ---- fig.width=8, fig.height=8------------------------------------------
LinePlot(X = DataSimul$x, createWindow = FALSE, main = "Line plot",  rows = c(1, 2),
              type = "p", num.stacked = 2, xlab = "x-axis", ylab = "y-axis",
              ang = "0", xaxis_type = "numerical", nxaxis = 10)

