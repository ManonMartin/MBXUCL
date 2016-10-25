#===============================================================
######################
# fonction OPLSDA
######################

# OPLS-DA for a two-class problem
#' @export OPLSDA
#' @title OPLS-DA for a two-class problem
#' @description OPLS-DA for a two-class problem. Only one response is allowed and the number of predictive components is fixed to 1 accordingly.
#'
#' @param x  A data matrix on which will be based the analysis.
#' @param y  A numerical vector representing the class of individuals.
#' @param impT If \code{TRUE}, prints the results.
#' @param impG If \code{TRUE}, save results' graphics in pdf format.
#' @param no Number of orthogonal components to keep.
#' @param out.path Path to output the results' graphics
#' @param nb Number of biomarkers to select.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{b.coef}}{Model coefficients}
#'   \item{\code{Tp}}{Predictive scores}
#'   \item{\code{Tortho}}{Orthogonal scores}
#'   \item{\code{Pp}}{Predictive loadings}
#'   \item{\code{Portho}}{Orthogonal loadings}
#'   \item{\code{Wortho}}{Orthogonal weights matrix}
#'   \item{\code{Selected.biomarkers}}{Vector of identified biomarkers}
#'   \item{\code{colmeans}}{Values of x column means}
#' }
#'
#' @examples
#' data("DataSimul")
#' x = DataSimul[["x"]]
#' y = DataSimul[["y"]]
#' oplsda.res = OPLSDA(x=x, y=y, impT = FALSE,impG = FALSE, no=2, nb = 15, out.path = ".")
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline barplot legend text title
#' @importFrom stats cov sd


OPLSDA <-function(x, y, impT=FALSE,impG=FALSE, no=2,
                                    nb=15, out.path = ".") {

# checks
if (sum(!y %in% c(0,1)) >0) {
  warning("y is not a vector of 1 and 0s")
}

checkArg(impT, "bool", can.be.null=FALSE)
checkArg(impG, "bool", can.be.null=FALSE)
checkArg(no, "int", can.be.null=FALSE)
checkArg(nb, "int", can.be.null=FALSE)
checkArg(out.path, "str", can.be.null=FALSE)

colmeans = apply(x, 2, mean)
x=x-matrix(apply(x,2,mean),nrow=dim(x)[1],
           ncol=dim(x)[2],byrow=T) # centrage de x sur les colonnes (par descripteur)


xoriginal=x
varnames=dimnames(x)[[2]]
obsnames = dimnames(x)[[1]]
xtrain = x
ytrain = y

m=dim(x)[2] # nombre de descripteurs
n=dim(x)[1]

options(warn=-1)
xax = as.numeric(colnames(x))
if (sum(is.na(xax))>0) {xax = c(1:m)}
options(warn=1)

# Initialisation des matrices de sortie de l'algorithme it\'eratif :

Tp = c()
Pp = c()
C = c()
W = c()
CV = c()
Tortho = c()
Portho = c()
Wortho = c()
VarXortho=c()
VarX=c()


# Debut de l'algorithme de l'OPLS-DA avec le calcul des poids de X :
#1.
w = (t(xtrain) %*% ytrain) %*% solve(t(ytrain) %*% ytrain) # X-weights w
w = as.vector(w)
#2.
norm.w = sqrt(w%*%w)                        # norm of w (vector)
w = w %*% solve(norm.w)                     # normalised w


# creation de la boucle pour extraire les composantes orthogonales
for (j in 1:no) {
  #3.
  t = (xtrain %*% w) %*% solve(t(w) %*% w)                   # X-scores
  #4.
  c = t(ytrain) %*% t %*% solve(t(t) %*% t)                 # y-weights
  #5.
  u = ytrain %*% c %*% solve(t(c) %*% c)                    # y-scores
  #6.
  p = (t(xtrain) %*% t) %*% solve(t(t) %*% t)                # 1st PLS loading
  p= as.vector(p)
  norm.p = sqrt(p%*%p)
  #7.
  wo = p - (w %*% (t(w) %*% p) %*% solve(t(w) %*% w))       # Xortho-weights
  wo = as.vector( wo)
  #8.
  norm.wo = sqrt(wo%*%wo)
  wo = wo %*% solve(norm.wo)                                # normalised wortho
  cv = norm.wo %*% solve(norm.p)              # criterion for the number of orthogonal components to keep
  #9.
  to = xtrain %*% wo %*% solve(t(wo) %*% wo)                 # Xortho-scores
  #10.
  po = t(xtrain) %*% to %*% solve(t(to) %*% to)              # ortho. loading

  varortho = to %*% t(po)
  var = t %*% t(p)

  #11.
  xtrain = xtrain - to %*% t(po)

  # matrice var-cov
  covo=cov(varortho) # otho

  covp=cov(xtrain) # pred
  # pls=mvr(ytrain ~ xtrain, ncomp=1)
  # pls$Xtotvar
  covorig=cov(cbind(xoriginal))

  #  total variation X
  varxortho=100*(sum(diag(covo))/sum(diag(covorig))) # var orthog retirée, en %
  varx=100*(sum(diag(covp))/sum(diag(covorig))) # var non orthog, en %

  #12.
  Tp = matrix(c(Tp, t))
  Pp = matrix(c(Pp, p))
  C = matrix(c(C, c))
  W = matrix(c(W, w))
  CV = matrix(c(CV, cv))
  Tortho = matrix(c(Tortho, to))
  Portho = matrix(c(Portho, po))
  Wortho = matrix(c(Wortho, wo))
  VarXortho = matrix(c(VarXortho,varxortho))
  VarX = matrix(c(VarX,varx))
}



Tp = matrix(Tp, ncol = no, dimnames = list(obsnames, NULL))
Pp = matrix(Pp, ncol = no, dimnames = list(varnames, NULL))
C = matrix(C, ncol = no)
W = matrix(W, ncol = no)

Tortho = matrix(Tortho, ncol = no, dimnames = list(obsnames, NULL))
Portho = matrix(Portho, ncol = no, dimnames = list(varnames, NULL))
Wortho = matrix(Wortho, ncol = no)

VarXortho = matrix(VarXortho, ncol = no)
VarX = matrix(VarX, ncol = no)

#13.
xortho = Tortho %*% t(Portho)



# Calcul des coefficients de regression b
b.coef=t(C[,no]%*%t(W[,no]))
b.coef = as.vector(b.coef)
names(b.coef)=varnames


# Recherche des nb biomarkeurs ayant les b.coef les plus grands
###############################################################
namindbiom=order(abs(b.coef))[((m+1)-(1:nb))]
indbiom = b.coef[namindbiom]


# sortie
#
ropls<-list(b.coef=b.coef,Tp = Tp[, no], Tortho = Tortho, Pp = Pp[, no],
            Portho = Portho, Wortho = Wortho, Selected.biomarkers =indbiom, colmeans = colmeans,
            CV = CV, original.dataset = xoriginal)


# Sorties graphiques
##############################

if(impG == TRUE) {

  # Validation plot
  ###########################
  pdf(file.path(out.path, "OPLS_Validation.pdf"), width = 10, height = 6)

  COL=rep("gray93",no)
  #par(mar=c(4,4,4,4))
  mp=barplot(t(CV), axes=F, axisnames=F, border=1,col=COL)
  axis(1,at=mp, labels=c(1:no))
  axis(2)
  title(main="OPLS: Choice of the n[orthog. Components]", xlab="Orthogonal OPLS-DA Components", ylab="||Wortho|| / ||p||")

  dev.off()


  # Xvar explique
  ###########################
  pdf(file.path(out.path, "OPLS_Xvar.pdf"), width = 6, height = 6)

  COL=rep("gray93",no+1)
  COL[1]="cyan1"
  mp1=graphics::barplot(c(VarX[no],VarXortho), axes=F, axisnames=F, width=rep(0.7, no), border=1, col=COL)
  axis(1,at=mp1, tick=F,labels=c("PC", paste0("OC",1:no)))
  axis(2)
  mtext(paste(round(c(VarX[no],VarXortho),1),"%", sep=" "),col="gray40", at =mp1,  side =1)
  title(main="OPLS: Percentage of X-variance \n explained by each component", xlab="OPLS Components", ylab="% X-variance")

  dev.off()

  # OPLS-DA score scatter plot
  ###########################
  index=ceiling(no/2)

  pdf(file.path(out.path,"OPLS_Score.pdf"), width = 12, height = 6*index)
  col = numeric()
  col[ytrain==0] = 3
  col[ytrain==1] = 2

  par(mfrow=c(ceiling(no/2),2))
  for (i in 1:no) {
    max.pc1 = 1.1 * (max(abs(Tp[, i])))
    max.pc2 = 1.1 * (max(abs(Tortho[, i])))

    lim = c()
    if (max.pc1 > max.pc2) {
      lim = c(-max.pc1, max.pc1)
    }else {
      lim = c(-max.pc2, max.pc2)
    }

    plot(Tp[, no], Tortho[, i], col=col, pch = 19, xlim =lim,
         ylim = lim, xlab = paste0("Predictive T score [",1,"]"), ylab = paste0("Orthogonal T score [",i,"]"),
         main = "OPLS-DA score scatter plot")
    abline(h=0, v=0, lty = 2, col = "gray")
    legend(x="topright", border=T, y.intersp=0.7, cex=1, pch = 19,
           legend=list("Group1", "Group0"), col=c(2,3))
  }

  dev.off()


  # S-plot
  ###########################
  ## correlation and covariance matrices
  # covariance
  s = as.matrix(x, ncol = ncol(x))
  p1 = c()
  for (i in 1:ncol(s)) {
    scov = cov(s[, i], Tp[,no])
    p1 = matrix(c(p1, scov), ncol = 1) # covariance x-T
  }
  # correlation
  pcorr1 = c()
  Tno=as.matrix(Tp[,no], ncol=1)
  for (i in 1:nrow(p1)) {
    den = apply(Tno, 2, sd, na.rm = TRUE) * sd(s[, i])
    corr1 = p1[i, ]/den
    pcorr1 = matrix(c(pcorr1, corr1), ncol = 1) # correlation
  }

  # plot
  pdf(file.path(out.path,"OPLS_Splot.pdf"), width = 10, height = 8)

  plot(p1, pcorr1, xlab = "p(cov)[1]", ylab = "p(corr)[1]",
       main = "S-plot (OPLS-DA)", ylim=c(min(pcorr1, na.rm = T)*1.1, max(pcorr1, na.rm = T)*1.1), xlim=c(min(p1, na.rm = T)*1.1, max(p1, na.rm = T)*1.1))
  sel = p1*pcorr1
  sel = order(sel, decreasing = TRUE)[1:nb]
  text(p1[sel], pcorr1[sel], labels = colnames(s)[sel], cex = 0.7, pos = 1)
  abline(v=0,lty = 2)
  abline(h=0,lty = 2)

  dev.off()


  ########################################
  ## plot of the biomarker coefficients ##
  ########################################

  # OPLS pretreated PLS coefficients
  delta=mean(sort(abs(b.coef))[m-nb+c(0,1)])

  pdf(file.path(out.path,"OPLS_coef.pdf"), width = 10, height = 6)

  par(mar=c(4,2,2,1))

  plot(b.coef,type="l",xaxt = "n", yaxt = "n", main="OPLS: Vector of descriptors' rank",xlab="ppm" )
  abline(h=0)
  abline(h=delta*c(-1,1),lty=2)
  axis(side=2,cex.axis=0.7)
  axis(side=1,at=c(1,seq(50,m,50)),labels=xax[c(1,seq(50,m,50))],cex.axis=0.7)

  dev.off()



  # loadings plot
  ##########################
  pdf(file.path(out.path,"OPLS_LoadingPred.pdf"), width = 10, height = 6)

  plot(Pp[,no], xaxt = "n", type="l", main="OPLS pretreated PLS coefficients (loadings)")
  axis(side=2,cex.axis=0.7)
  axis(side=1,at=c(1,seq(50,m,50)),labels=xax[c(1,seq(50,m,50))],cex.axis=0.7)

  dev.off()

}

###########################
# text impression
###########################
if(impT==TRUE)
{
  #  cat("\n Dataset",  dataname)
  print((ropls))
}
return(ropls)
}



#===============================================================
#===============================================================
##############################
### Fonction OPLSDA_pred
##############################
#' @export OPLSDA_pred
#' @title Prediction for OPLS-DA
#'
#' @description
#' Prediction of a new set of observations based on a built model.
#'
#' @param ropls Result from OPLS-DA analysis with the OPLSDA function.
#' @param x.new A vector or matrix of new observations.
#'
#' @return A vector with predicted y values.
#'
#' @examples
#' data("DataSimul")
#' x = DataSimul[["x"]]
#' y = DataSimul[["y"]]
#' oplsda.res = OPLSDA(x=x[-c(1:5),], y=y[-c(1:5)],
#'          impT = FALSE,impG = FALSE, no=2, nb = 15, out.path = ".")
#' OPLSDA_pred(ropls = oplsda.res, x.new = x[c(1:5),])


OPLSDA_pred = function(ropls, x.new) {


  no = dim(ropls$Wortho)[2]
  # centrage de etnew
  x.new = x.new-ropls$colmeans

  l = dim(x.new)[1]

  if (is.null(l)) {
    for (i in 1:no) {
      tnew.ortho = t(x.new) %*% ropls$Wortho[,i] %*% solve(t(ropls$Wortho[,i]) %*% ropls$Wortho[,i])
      x.new = as.vector(t(x.new) - tnew.ortho %*% t(ropls$Portho[,i]))
    }
    #Prediction
    y.pred[k] = ropls$b.coef%*%x.new
  } else {
    y.pred =c()
    for(k in 1:l) {
      # Application du filtre sur la nouvelle observation
      for (i in 1:no) {
        tnew.ortho = t(x.new[k,]) %*% ropls$Wortho[,i] %*% solve(t(ropls$Wortho[,i]) %*% ropls$Wortho[,i])
        x.new[k,] = as.vector(t(x.new[k,]) - tnew.ortho %*% t(ropls$Portho[,i]))
      }
      #Prediction
      y.pred[k] = ropls$b.coef%*%x.new[k,]
    }
    names(y.pred) = rownames(x.new)
  }

  y.pred

}




#===============================================================

##############################
### Fonction cvOPLSDA
##############################
#' @export cvOPLSDA
#' @title K-fold cross-validation for OPLS-DA.
#'
#' @description
#' K-fold cross-validation for OPLS-DA based on the RMSE criterion and stratified according to y.
#'
#' @param x A data matrix on which will be based the analysis.
#' @param y A numerical vector representing the class of individuals.
#' @param k_fold The number of sub-datasets to create for cross-validation.
#' @param NumOrtho The maximum number of orthogonal components allowed in OPLSDA.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{RMSECV}}{Result from OPLS-DA analysis with the OPLSDA function.}
#'   \item{\code{folds_i}}{A vector indicating the group of observations for cross-validation.}
#'   \item{\code{Prop1}}{Table of proportions of 1's in y for each subgroup.}
#' }
#'
#' @examples
#' data("DataSimul")
#' x = DataSimul[["x"]]
#' y = DataSimul[["y"]]
#' cvOPLSDA.res = cvOPLSDA(x = x, y = y, k_fold = 10, NumOrtho = 5)
#' plot(cvOPLSDA.res$RMSECV)
#'
#' @importFrom plyr ddply

cvOPLSDA = function(x, y, k_fold = 10, NumOrtho = 1){

  # checks
  if (sum(!y %in% c(0,1)) >0) {
    warning("y is not a vector of 1 and 0s")
  }

  checkArg(k_fold, "int", can.be.null=FALSE)
  checkArg(NumOrtho, "int", can.be.null=FALSE)



  library("plyr")
  df = data.frame(rowname = rownames(x), class = y)
  n = dim(x)[1]

  createFolds <- function(x,k){
    n <- nrow(x)
    x$folds <- rep(1:k,length.out = n)[sample(n,n)]
    x
  }



  folds <- plyr::ddply(.data = df, .variables = .(class),.fun = plyr::here(createFolds),k = k_fold)
  folds_i = folds$folds

  Prop1 = plyr::ddply(.data = folds,.variables = .(folds),.fun = plyr::here(plyr::summarise),prop = sum(class)/length(class))



  NumOrtho = 1:NumOrtho

  index = vector("list", k_fold)
  Xtrain = vector("list", k_fold)
  Ytrain = vector("list", k_fold)
  Xnew = vector("list", k_fold)
  for (k in 1:k_fold) {
    index[[k]] = which(folds_i==k)
    Xtrain[[k]]= x[-index[[k]],]
    Ytrain[[k]] = y[-index[[k]]]
    Xnew[[k]] = x[index[[k]],]
  }

  RMSECV = c()
  for (j in 1:length(NumOrtho)) {

    ropls.train = mapply(OPLSDA, x = Xtrain,
           y = Ytrain, MoreArgs = list(impT=FALSE,
           impG=FALSE, no=j, out.path = "."), SIMPLIFY = FALSE,
           USE.NAMES = TRUE)


    ypred = mapply(OPLSDA_pred, ropls = ropls.train, x.new = Xnew,
                 SIMPLIFY = FALSE, USE.NAMES = TRUE)

    Ypred = unlist(ypred)
    sorted = order(as.numeric(names(Ypred)))
    Ypred = Ypred[sorted]

    RMSECV[j] = sqrt(sum((y-Ypred)^2)/n)

  }
  rcvOPLSDA = list(RMSECV = RMSECV, folds_i = folds_i, Prop1 = Prop1)
  return(rcvOPLSDA)

}



