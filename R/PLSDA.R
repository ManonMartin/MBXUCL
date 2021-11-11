
#' @export PLSDA
#' @title PLS-DA
#'
#' @description
#' PLS-DA model for multiple classes (2 or more).
#'
#' @param x A numeric matrix of spectral intensities.
#' @param y The real groups' memberships of the samples.
#' @param nLV The number of latent variables for PLS-DA.
#' @param validation The type of validation (see `pls::plsr()` for more info).
#' @param ... Additional arguments to be passed to `pls::plsr()` (eg: `segments`).
#'
#' @param drawRMSEP Logical, if \code{'TRUE'}, will draw plot of RMSEP values to decide on the number of latent variables.
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
#'   \item{\code{VIP}}{The variable VIPs}
#'    }
#' @author Manon Martin
#'
#' @examples
#' data('HumanSerum')
#' 
#' PLSDA = PLSDA(x = HumanSerumSpectra,  y = ClassHS, nLV = NULL, drawRMSEP = TRUE, validation = "LOO")
#' 
#' PLSDA = PLSDA(x = HumanSerumSpectra,  y = ClassHS, nLV = NULL, drawRMSEP = TRUE, validation = "CV", segments = 5)
#'
#'@importFrom pls cppls mvrValstats RMSEP R2 MSEP
#'@importFrom stats model.matrix
#'@importFrom graphics axis mtext par plot
#'@importFrom plsVarSel VIP



PLSDA <- function(x, y, nLV = NULL, drawRMSEP = TRUE, 
                  validation = c("CV", "none", "LOO"), ...) {


  checkArg(nLV, "num", can.be.null = TRUE)
  checkArg(drawRMSEP, "bool", can.be.null = FALSE)
  validation <- match.arg(validation)

  ######################### I. Initialisations #

  if (is.null(nLV)) {
    nLV <- 3
  }

  m <- dim(x)[2]  # nombre de descripteurs
  n <- dim(x)[1]

  xtrain <- as.matrix(x)
  ytrain <- y

  options(warn = -1)
  xax <- round(as.numeric(names(x[1, ])), 2)
  if (sum(is.na(xax)) > 0) {
    xax <- c(1:m)
  }
  options(warn = 1)


  # Create a dummy response matrix
  dummy <- I(stats::model.matrix(~y - 1, data.frame(y = factor(ytrain))))
  dummy
  ytrain.mat <- matrix(dummy, nrow = length(y))

  if (nlevels(factor(ytrain))==2){
    ytrain.mat <- ytrain.mat[,1,drop=FALSE]
  }

  ################################################# II.PLS-DA #

  # choix du nombre de variables latentes : sur le training set par
  # cross-validation


  # ------------------- Cross-validation ------------------- -------------------

  pls1 <- pls::plsr(ytrain.mat ~ xtrain, scale = FALSE, validation = validation,
                    model = TRUE, ...)

  # ncomp = min(5, pls1$ncomp) pls::mvrValstats(pls1, estimate='train')

  # % of variance explained ------------------------
  explvarYcum <- 100 * drop(pls::R2(pls1, estimate = "train", intercept = FALSE)$val)
  
  if (nlevels(factor(ytrain))==2){
    explvarYcum <- t(as.matrix(explvarYcum))
    rownames(explvarYcum) <- "Y1"
  } 
  explvarY <- explvarYcum
  
  
  for (i in 2:dim(explvarY)[2]) {
    explvarY[, i] <- explvarYcum[, i] - explvarYcum[, (i - 1)]
  }
  ExpVarY <- list(explvarY = explvarY, explvarYcum = explvarYcum)


  explvarX <- pls::explvar(pls1)  # X variance explained by the components
  explvarXcum <- cumsum(pls::explvar(pls1))

  ExpVarX <- list(explvarX = explvarX, explvarXcum = explvarXcum)


  if (drawRMSEP == TRUE) {
    
    nvals <- min(20,pls1$ncomp)
    
    # R^2: coefficient of multiple determination
    Q2 <- drop(pls::R2(pls1, estimate = "CV", intercept = FALSE)$val) 
    if (nlevels(factor(ytrain))==2){
      Q2 <- t(as.matrix(Q2))
      rownames(Q2) <- "Y1"
    } 
    
    Q2 <- Q2[,1:nvals, drop=FALSE]
    Q2 <- t(Q2)
    
    # RMSEP
    RMSEP <- drop(pls::RMSEP(pls1, intercept = FALSE, estimate = "CV")$val)
    if (nlevels(factor(ytrain))==2){
      RMSEP <- t(as.matrix(RMSEP))
      rownames(RMSEP) <- "Y1"
    } 
    
    RMSEP <- RMSEP[,1:nvals, drop=FALSE]
    RMSEP <- t(RMSEP)
    
    if (dim(RMSEP)[2] < 3) {
      par(mfrow = c(1, dim(RMSEP)[2]), oma = c(0, 0, 4, 0), mar = c(3, 4, 2,
        4) + 0.1)
    } else {
      par(mfrow = c(ceiling(dim(RMSEP)[2]/3), 3), oma = c(0, 0, 4, 0), mar = c(3,
        4, 2, 4) + 0.1)
    }

    for (i in 1:dim(RMSEP)[2]) {
      plot(1:nvals, Q2[1:nvals,i], type = "b", ylab = "Q2", xlab = "", main = paste0("Y",i))
      abline(v=nLV)
      
      par(new = TRUE)
      plot(1:nvals, RMSEP[1:nvals,i], type = "b", lty = 2, col = "red", axes = FALSE,
            xlab = "", ylab = "")
      abline(v=nLV)
      mtext("RMSECV", side = 4, col = "red", line = 3)
      axis(4, col = "red", col.axis = "red", las = 1)

    }
    mtext("PLS: Cross-validated RMSEP curves", outer = TRUE, cex = 1.5)
  }


  # ----------------------------------------------------- Estimation du modele
  # final avec nLV composantes
  # -----------------------------------------------------
  # -----------------------------------------------------

  # Estimation du modele final
  pls <- pls::plsr(ytrain.mat ~ xtrain, ncomp = nLV, scale = FALSE, 
                   validation = validation, ...)


  Xscores <- cbind(pls$scores)

  # Estimation de la capacite predictive R^2: coefficient of multiple determination
  R2 <- pls::R2(pls, estimate = "train", intercept = FALSE)$val[1, , nLV]  #'train': training or calibration data estimate

  # RMSEP
  RMSEP <- pls::RMSEP(pls, estimate = "CV", intercept = FALSE)$val[1, , nLV]

  # Q2 : cross-validated R^2

  ## cross-validated R^2:
  Q2 <- pls::R2(pls, estimate = "CV", intercept = FALSE)$val[1, , nLV]

  

  # PRESS = apply(pls::mvrValstats(pls, estimate='CV')$SSE[1,,],2,sum) # same PRESS
  # values

  # SSE
  pls::mvrValstats(pls, estimate = "train")$SSE[1, , ]  #same residual sum of squares: residuals = pls$residuals[,,(nLV-1)]^2
  
  if (nlevels(factor(ytrain))==2){
    Q2cum <- Q2
  } else{ 
    # ===== PRESS
    PRESS <- apply(pls$validation$PRESS, 2, sum)
    
    SS <- apply(pls::mvrValstats(pls, estimate = "train")$SSE[1, , ], 2, sum) 
    SS <- SS[-nLV - 1]
    Q2cum <- 1 - prod(PRESS/SS)
  }

  

  VIP <- plsVarSel::VIP(pls, opt.comp = nLV)

  # Preparation des sorties ----------------------------------------
  # ----------------------------------------

  #
  rpls <- list(original.dataset = xtrain, Y = ytrain.mat, nLV = nLV, RMSEP = RMSEP,
              R2 = R2, Q2 = Q2, Q2cum = Q2cum, ExpVarY = ExpVarY, ExpVarX = ExpVarX,
              coefficients = pls$coefficients, scores = pls$scores, loadings = pls$loadings,
              VIP = VIP)

  return(rpls)
}
