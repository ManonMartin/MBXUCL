
#' @export PLSDA
#' @title PLS-DA
#'
#' @description
#' PLS-DA model for multiple classes (2 or more).
#'
#' @param x A numeric matrix of spectral intensities.
#' @param y The real groups' memberships of the samples.
#' @param nLV The number of latent variables for PLS-DA.
#'
#' @param drawRMSEP Logical, if \code{"TRUE"}, will draw plot of RMSEP values to decide on the number of latent variables.
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{original.dataset}}{Original matrix of spectral intensities}
#'   \item{\code{Y}}{Dummy matrix representing the classes}
#'   \item{\code{nLV}}{Number of latent variables kept. If \code{NULL}, will be set to 3 in the finla model.}
#'   \item{\code{RMSEP}}{Root mean squared error of prediction }
#'   \item{\code{R2}}{Coefficient of multiple determination}
#'   \item{\code{Q2}}{Cross-validated coefficient of multiple determination}
#'   \item{\code{Q2cum}}{Cumulative cross-validated coefficient of multiple determination (over all response variables)}
#'   \item{\code{ExpVarY}}{Proportion of Y variance explained by the model}
#'   \item{\code{ExpVarX}}{Proportion of X variance explained by the model}
#'   \item{\code{coefficients}}{PLSDA coefficients}
#'   \item{\code{scores}}{PLSDA scores}
#'   \item{\code{loadings}}{PLSDA loadings}
#'    }
#' @author Manon Martin
#'
#' @examples
#' data("HumanSerum")
#' PLSDA = function(x = HumanSerumSpectra, y = ClassHS, nLV = NULL, drawRMSEP = TRUE)
#'
#'@importFrom pls cppls mvrValstats RMSEP R2 MSEP
#'@importFrom stats model.matrix
#'@importFrom graphics axis mtext par plot



PLSDA = function(x, y, nLV=NULL, drawRMSEP = TRUE){


checkArg(nLV, "num", can.be.null=TRUE)
checkArg(drawRMSEP, "bool", can.be.null=FALSE)


#########################
#########################
#
#  I. Initialisations   #
#
#########################
#########################

if (is.null(nLV)){
  nLV=3
}


# x=x-matrix(apply(x,2,mean),nrow=dim(x)[1],ncol=dim(x)[2],byrow=T) # centrage de x sur les colonnes (par descripteur)
m=dim(x)[2] # nombre de descripteurs
n=dim(x)[1]

Id=1:n
IdSpectra=dimnames(x)[[1]]

levgroup=as.numeric(names(table(y)))
nlevgroup=table(y)
k=length(levgroup)

# ord=order(y)
# x=x[ord,]
# group=y[ord]
# IdSpectra=IdSpectra[ord]

n=c()
id=list()
for (i in 1:k){
  n[i] = nlevgroup[[i]]
  id[[i]] = which(y==levgroup[i])
}


xtrain=as.matrix(x)
ytrain=y

xax=round(as.numeric(names(x[1,])),2)
names=dimnames(x)[[2]]

# Create a dummy response matrix
dummy = I(stats::model.matrix(~y-1, data.frame(y = factor(ytrain))))
dummy
ytrain.mat = matrix(dummy, nrow = length(y) )


#################################################
#################################################
#
#                   II.PLS-DA                   #
#
#################################################
#################################################

# choix du nombre de variables latentes :
# sur le training set par cross-validation


# -------------------
# Cross-validation
# -------------------
# -------------------

pls1 = pls::cppls(ytrain.mat ~ xtrain, scale = FALSE, validation = "CV", model = TRUE)

# ncomp = min(5, pls1$ncomp)
# pls::mvrValstats(pls1, estimate="train")

#% of variance explained
# ------------------------
explvarYcum=100 * drop(pls::R2(pls1, estimate = "train", intercept = FALSE)$val)
explvarY=explvarYcum
for (i in 2:dim(explvarY)[2]) {
  explvarY[,i]=explvarYcum[,i]-explvarYcum[,(i-1)]
}
ExpVarY = list(explvarY = explvarY, explvarYcum = explvarYcum)


explvarX=pls::explvar(pls1) # X variance explained by the components
explvarXcum=cumsum(pls::explvar(pls1))

ExpVarX = list(explvarX = explvarX, explvarXcum = explvarXcum)


if (drawRMSEP==TRUE){
  #R^2: coefficient of multiple determination
  R2 = drop(pls::R2(pls1, estimate="train", intercept = FALSE)$val) #"train": training or calibration data estimate

  # RMSEP
  RMSEP = drop(pls::RMSEP(pls1, intercept = FALSE)$val)
  if (dim(RMSEP)[2]<3) {par(mfrow=c(1,dim(RMSEP)[2]),oma = c(0, 0, 4, 0),mar=c(5, 4, 4, 4) + 0.1)
  }else{ par(mfrow=c(ceiling(dim(RMSEP)[2]/3),3),oma = c(0, 0, 4, 0),mar=c(5, 4, 4, 4) + 0.1)}

  for (i in 1:dim(RMSEP)[2]){
    plot(R2[i,1:10], type="l", ylab="R2", xlab="", main=paste0("Y", i))

    par(new=TRUE)
    plot(RMSEP[1,i,1:10], type="l", lty = 2, col="red", axes=FALSE, xlab="", ylab="")
    mtext("RMSEP",side=4,col="red",line=3)
    axis(4, col="red",col.axis="red",las=1)

  }
  mtext("PLS: Cross-validated RMSEP curves", outer = TRUE, cex = 1.5)
  }


# -----------------------------------------------------
#    Estimation du modele final avec nLV composantes
# -----------------------------------------------------
# -----------------------------------------------------

# Estimation du modele final
#
pls = pls::cppls(ytrain.mat ~ xtrain, ncomp = nLV, scale = FALSE, validation = "CV")


Xscores=cbind(pls$scores)

# Estimation de la capacite predictive
#
#R^2: coefficient of multiple determination
R2 = pls::R2(pls, estimate="train", intercept = FALSE)$val[1,,nLV] #"train": training or calibration data estimate

# RMSEP
RMSEP = pls::RMSEP(pls, intercept = FALSE)$val[1,,nLV]

#Q2 : cross-validated R^2

## cross-validated R^2:
Q2 = pls::R2(pls, estimate = "CV", intercept = FALSE)$val[1,,nLV]

#=====
# PRESS
PRESS = apply(pls$validation$PRESS,2, sum)

# PRESS = apply(pls::mvrValstats(pls, estimate="CV")$SSE[1,,],2,sum) # same PRESS values

# SSE
pls::mvrValstats(pls, estimate="train")$SSE[1,,] #same residual sum of squares: residuals = pls$residuals[,,(nLV-1)]^2
SS = apply(pls::mvrValstats(pls, estimate="train")$SSE[1,,],2,sum) 
SS = SS[-nLV-1]


Q2cum = 1- prod(PRESS/SS)



#     Preparation des sorties
# ----------------------------------------
# ----------------------------------------

#
rpls<-list(original.dataset = xtrain,  Y= ytrain.mat, nLV=nLV,
           RMSEP = RMSEP, R2 = R2, Q2 = Q2, Q2cum = Q2cum,  ExpVarY = ExpVarY, ExpVarX = ExpVarX,
           coefficients =  pls$coefficients, scores = pls$scores, loadings = pls$loadings)

return(rpls)
}
