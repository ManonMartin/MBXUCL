# =============================================================== fonction OPLSDA

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
#'   \item{\code{b}}{Model coefficients}
#'   \item{\code{Tp}}{Predictive scores}
#'   \item{\code{Pp}}{Predictive loadings}
#'   \item{\code{W}}{X-weights}
#'   \item{\code{C}}{y-weights}
#'   \item{\code{Tortho}}{Orthogonal scores}
#'   \item{\code{Portho}}{Orthogonal loadings}
#'   \item{\code{Wortho}}{Orthogonal weights matrix}
#'   \item{\code{Selected.biomarkers}}{Vector of identified biomarkers}
#'   \item{\code{CV}}{Criterion for the number of orthogonal components to keep}
#'   \item{\code{original.dataset}}{Original dataset}
#'   \item{\code{Xopls}}{OPLS-filtered X matrix}
#' }
#'
#' @details
#' The function allows only one predictive component since it is designed for a single dependent variable.
#' It is based on the NIPALS PLS algorithm: after the removal of orthogonal components from the X matrix, a PLS1
#' NIPALS is ran on the filtered X matrix.
#'
#' @examples
#' data('DataSimul')
#' x = DataSimul[['x']]
#' y = DataSimul[['y']]
#' oplsda.res = OPLSDA(x=x, y=y, impT = FALSE,impG = FALSE, no=2, nb = 15, out.path = '.')
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline barplot legend text title
#' @importFrom stats cov sd


OPLSDA <- function(x, y, impT = FALSE, impG = FALSE, no = 2, nb = 15, out.path = ".") {

  # checks

  # if (sum(!y %in% c(0,1)) >0) { y = as.factor(y) if (nlevels(y) >2) { stop('There
  # is more than 2 levels in y') } else { warning('levels of y are automaticaly
  # transformed in 0 and 1s') y[y==levels(y)[1]] = 0 y[y==levels(y)[2]] = 1 } }

  y <- as.numeric(y)

  checkArg(impT, "bool", can.be.null = FALSE)
  checkArg(impG, "bool", can.be.null = FALSE)
  checkArg(no, "int", can.be.null = FALSE)
  checkArg(nb, "int", can.be.null = FALSE)
  checkArg(out.path, "str", can.be.null = FALSE)

  colmeans <- apply(x, 2, mean)
  x <- x - matrix(apply(x, 2, mean), nrow = dim(x)[1], ncol = dim(x)[2], byrow = T)  # centrage de x sur les colonnes (par descripteur)


  xoriginal <- x
  varnames <- dimnames(x)[[2]]
  obsnames <- dimnames(x)[[1]]
  xtrain <- x
  ytrain <- y

  m <- dim(x)[2]  # nombre de descripteurs
  n <- dim(x)[1]

  if (nb > m){
    stop("nb is superior to the total number of descriptors, you need to reduce its value")
  }


  options(warn = -1)
  xax <- as.numeric(colnames(x))
  if (sum(is.na(xax)) > 0) {
    xax <- c(1:m)
  }
  options(warn = 1)


  # pls_nipals function ------------------------------

  pls_nipals <- function(xtrain, ytrain, np = 1) {

    xtrain <- scale(xtrain, center = T, scale = F)

    n <- dim(xtrain)[1]
    m <- dim(xtrain)[2]

    if (np < 1)  {
      stop("np should be at least 1")
    }

    Tp <- matrix(NA, ncol = np, nrow = n, dimnames = list(obsnames, NULL))
    Pp <- matrix(NA, ncol = np, nrow = m, dimnames = list(varnames, NULL))
    C <- matrix(NA, nrow = 1, ncol = np)
    W <- matrix(NA, ncol = np, nrow = m)

    for (j in 1:np)  {
      W[, j] <- as.vector((t(xtrain) %*% ytrain) %*% solve(t(ytrain) %*% ytrain))  # X-weights w
      norm.w <- sqrt(W[, j] %*% W[, j])  # norm of w (vector)
      W[, j] <- W[, j] %*% solve(norm.w)  # normalised w
      Tp[, j] <- (xtrain %*% W[, j]) %*% solve(t(W[, j]) %*% W[, j])  # X-scores
      C[, j] <- t(ytrain) %*% Tp[, j] %*% solve(t(Tp[, j]) %*% Tp[, j])  # y-weights
      Pp[, j] <- as.vector((t(xtrain) %*% Tp[, j]) %*% solve(t(Tp[, j]) %*%
        Tp[, j]))  # 1st PLS loading
      xtrain <- xtrain - Tp[, j] %*% t(Pp[, j])
    }

    b <- W %*% solve(t(Pp) %*% W) %*% as.vector(C)

    return(list(Tp = Tp, Pp = Pp, C = C, W = W, b = b))
  }


  # Initialisation des matrices de sortie de l'algorithme it\'eratif :
  Tortho <- matrix(NA, nrow = n, ncol = no, dimnames = list(obsnames, NULL))
  Portho <- matrix(NA, nrow = m, ncol = no, dimnames = list(varnames, NULL))
  Wortho <- matrix(NA, nrow = m, ncol = no)
  CV <- matrix(NA, nrow = 1, ncol = no)
  VarXortho <- matrix(NA, nrow = 1, ncol = no)
  VarX <- matrix(NA, nrow = 1, ncol = no)


  # Debut de l'algorithme de l'OPLS-DA avec le calcul des poids de X :
  w <- (t(xtrain) %*% ytrain) %*% solve(t(ytrain) %*% ytrain)  # X-weights w
  w <- as.vector(w)
  norm.w <- sqrt(w %*% w)  # norm of w (vector)
  w <- w %*% solve(norm.w)  # normalised w


  # creation de la boucle pour extraire les composantes orthogonales

  for (i in 1:no) {
    t <- (xtrain %*% w) %*% solve(t(w) %*% w)  # X-scores
    c <- t(ytrain) %*% t %*% solve(t(t) %*% t)  # y-weights
    u <- ytrain %*% c %*% solve(t(c) %*% c)  # y-scores
    p <- as.vector((t(xtrain) %*% t) %*% solve(t(t) %*% t))  # 1st PLS loading
    norm.p <- sqrt(p %*% p)
    Wortho[, i] <- as.vector(p - (w %*% (t(w) %*% p) %*% solve(t(w) %*% w)))  # Xortho-weights
    norm.wo <- sqrt(Wortho[, i] %*% Wortho[, i])
    Wortho[, i] <- Wortho[, i] %*% solve(norm.wo)  # normalised wortho
    CV[, i] <- norm.wo %*% solve(norm.p)  # criterion for the number of orthogonal components to keep
    Tortho[, i] <- xtrain %*% Wortho[, i] %*% solve(t(Wortho[, i]) %*% Wortho[,
      i])  # Xortho-scores
    Portho[, i] <- t(xtrain) %*% Tortho[, i] %*% solve(t(Tortho[, i]) %*% Tortho[,
      i])  # ortho. loading

    varortho <- Tortho[, i] %*% t(Portho[, i])
    var <- t %*% t(p)

    xtrain <- xtrain - Tortho[, i] %*% t(Portho[, i])

    # matrice var-cov
    covo <- cov(varortho)  # otho
    covp <- cov(xtrain)  # pred
    covorig <- cov(cbind(xoriginal))

    # total variation X
    ## var orthog retiree, en %
    VarXortho[, i] <- 100 * (sum(diag(covo))/sum(diag(covorig)))
    ## var non orthog, en %
    VarX[, i] <- 100 * (sum(diag(covp))/sum(diag(covorig)))

  }


  # X orthogonal
  xortho <- Tortho %*% t(Portho)


  # Apply nipals pls on filtered X matrix (Xnew)
  xnew <- xtrain  # deflated matrix
  res_nipals <- pls_nipals(xnew, ytrain, np = 1)

  invisible(list2env(res_nipals, envir = environment()))

  Tp = res_nipals[["Tp"]]
  Pp = res_nipals[["Pp"]]
  C = res_nipals[["C"]]
  W = res_nipals[["W"]]
  b = res_nipals[["b"]]

  b <- as.vector(b)

  # In order to apply b directly on a new X, b becomes bcorr
  bcorr <- (diag(1, m, m) - Wortho %*% solve(t(Portho) %*% Wortho) %*% t(Portho)) %*%  b
  names(bcorr) <- varnames
  bcorr <- as.vector(bcorr)


  # Recherche des nb biomarkeurs ayant les b les plus grands
  namindbiom <- order(abs(b))[((m + 1) - (1:nb))]
  indbiom <- b[namindbiom]


  # sortie
  ropls <- list(b = bcorr, Tp = Tp, Pp = Pp, W = W, C = C, Tortho = Tortho, Portho = Portho,
    Wortho = Wortho, Selected.biomarkers = indbiom, CV = CV, original.dataset = xoriginal,
    Xopls = xnew)


  # Sorties graphiques

  if (impG == TRUE) {

    # Validation plot
    # pdf(file.path(out.path, "OPLS_Validation.pdf"), width = 10, height = 6)
    #
    # COL <- rep("gray93", no)
    # # par(mar=c(4,4,4,4))
    # mp <- barplot(t(CV), axes = F, axisnames = F, border = 1, col = COL)
    # axis(1, at = mp, labels = c(1:no))
    # axis(2)
    # title(main = "OPLS: Choice of the n[orthog. Components]", xlab = "Orthogonal OPLS-DA Components",
    #   ylab = "||Wortho|| / ||p||")
    #
    # dev.off()


    # Xvar explique
    pdf(file.path(out.path, "OPLS_Xvar.pdf"), width = 6, height = 6)

    COL <- rep("gray93", no + 1)
    COL[1] <- "cyan1"
    mp1 <- graphics::barplot(c(VarX[no], VarXortho), axes = F, axisnames = F,
                            width = rep(0.7, no), border = 1, col = COL)
    axis(1, at = mp1, tick = F, labels = c("PC", paste0("OC", 1:no)))
    axis(2)
    mtext(paste(round(c(VarX[no], VarXortho), 1), "%", sep = " "), col = "gray40",
          at = mp1, side = 1)
    title(main = "OPLS: Percentage of X-variance \n explained by each component",
          xlab = "OPLS Components", ylab = "% X-variance")

    dev.off()

    # OPLS-DA score scatter plot
    index <- ceiling(no/2)

    pdf(file.path(out.path, "OPLS_Score.pdf"), width = 12, height = 6 * index)
    col <- numeric()
    col[ytrain == 0] <- 3
    col[ytrain == 1] <- 2

    par(mfrow = c(ceiling(no/2), 2))
    for (i in 1:no) {
      max.pc1 <- 1.1 * (max(abs(Tp[, i])))
      max.pc2 <- 1.1 * (max(abs(Tortho[, i])))

      lim <- c()
      if (max.pc1 > max.pc2) {
        lim <- c(-max.pc1, max.pc1)
      } else  {
        lim <- c(-max.pc2, max.pc2)
      }

      plot(Tp[, no], Tortho[, i], col = col, pch = 19, xlim = lim, ylim = lim,
          xlab = paste0("Predictive T score [", 1, "]"),
          ylab = paste0("Orthogonal T score [",i, "]"), main = "OPLS-DA score scatter plot")
      abline(h = 0, v = 0, lty = 2, col = "gray")
      legend(x = "topright", border = T, y.intersp = 0.7, cex = 1, pch = 19,
            legend = list("Group1", "Group0"), col = c(2, 3))
    }

    dev.off()


    # S-plot correlation and covariance matrices covariance
    s <- as.matrix(x, ncol = ncol(x))
    p1 <- c()
    for (i in 1:ncol(s))  {
      scov <- cov(s[, i], Tp[, no])
      p1 <- matrix(c(p1, scov), ncol = 1)  # covariance x-T
    }
    # correlation
    pcorr1 <- c()
    Tno <- as.matrix(Tp[, no], ncol = 1)
    for (i in 1:nrow(p1)) {
      den <- apply(Tno, 2, sd, na.rm = TRUE) * sd(s[, i])
      corr1 <- p1[i, ]/den
      pcorr1 <- matrix(c(pcorr1, corr1), ncol = 1)  # correlation
    }

    # plot
    pdf(file.path(out.path, "OPLS_Splot.pdf"), width = 10, height = 8)

    plot(p1, pcorr1, xlab = "p(cov)[1]", ylab = "p(corr)[1]", main = "S-plot (OPLS-DA)",
        ylim = c(min(pcorr1, na.rm = T) * 1.1, max(pcorr1, na.rm = T) * 1.1),
        xlim = c(min(p1, na.rm = T) * 1.1, max(p1, na.rm = T) * 1.1))
    sel <- p1 * pcorr1
    sel <- order(sel, decreasing = TRUE)[1:nb]
    text(p1[sel], pcorr1[sel], labels = colnames(s)[sel], cex = 0.7, pos = 1)
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)

    dev.off()


    ######################################## plot of the biomarker coefficients ##

    # OPLS pretreated PLS coefficients
    delta <- mean(sort(abs(b))[m - nb + c(0, 1)])

    pdf(file.path(out.path, "OPLS_coef.pdf"), width = 10, height = 6)

    par(mar = c(4, 2, 2, 1))

    plot(b, type = "l", xaxt = "n", yaxt = "n", main = "OPLS: Vector of descriptors' rank", xlab = "ppm")
    abline(h = 0)
    abline(h = delta * c(-1, 1), lty = 2)
    axis(side = 2, cex.axis = 0.7)
    axis(side = 1, at = c(1, seq(50, m, 50)), labels = xax[c(1, seq(50, m, 50))],
      cex.axis = 0.7)

    dev.off()



    # loadings plot
    pdf(file.path(out.path, "OPLS_LoadingPred.pdf"), width = 10, height = 6)

    plot(Pp[, no], xaxt = "n", type = "l", main = "OPLS pretreated PLS coefficients (loadings)")
    axis(side = 2, cex.axis = 0.7)
    axis(side = 1, at = c(1, seq(50, m, 50)), labels = xax[c(1, seq(50, m, 50))],
      cex.axis = 0.7)

    dev.off()

  }

  ########################### text impression
  if (impT == TRUE)
  { # cat('\n Dataset', dataname)
    print((ropls))
  }
  return(ropls)
}



# ===============================================================
# =============================================================== Fonction OPLSDA_pred
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
#' data('DataSimul')
#' x = DataSimul[['x']]
#' y = DataSimul[['y']]
#' oplsda.res = OPLSDA(x=x[-c(1:5),], y=y[-c(1:5)],
#'          impT = FALSE,impG = FALSE, no=2, nb = 15, out.path = '.')
#' OPLSDA_pred(ropls = oplsda.res, x.new = x[c(1:5),])


OPLSDA_pred <- function(ropls, x.new) {
  y.pred <- x.new %*% ropls$b
  names(y.pred) <- rownames(x.new)
  return(y.pred)
}


# ===============================================================

############################## Fonction cvOPLSDA
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
#' @param ImpG If \code{TRUE}, prints a validation plot.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{RMSECV}}{RMSECV from the OPLSDA cross-validation.}
#'   \item{\code{folds_i}}{A vector indicating the group of observations for cross-validation.}
#'   \item{\code{Prop1}}{Table of proportions of 1's in y for each subgroup.}
#'   \item{\code{yres}}{List with true and predicted classes for each \code{NumOrtho}.}
#' }
#'
#' @examples
#' data('DataSimul')
#' x = DataSimul[['x']]
#' y = DataSimul[['y']]
#' cvOPLSDA.res = cvOPLSDA(x = x, y = y, k_fold = 10, NumOrtho = 5)
#' plot(cvOPLSDA.res$RMSECV)
#'
#' @importFrom plyr ddply

cvOPLSDA <- function(x, y, k_fold = 10, NumOrtho = 1, ImpG = FALSE) {

  # checks if (sum(!y %in% c(0,1)) >0) { warning('y is not a vector of 1 and 0s') }

  checkArg(k_fold, "int", can.be.null = FALSE)
  checkArg(NumOrtho, "int", can.be.null = FALSE)

  library("plyr")

  df <- data.frame(rowname = rownames(x), Class = y)
  n <- dim(x)[1]

  createFolds <- function(x, k) {
    n <- nrow(x)
    x$FOLDS <- rep(1:k, length.out = n)[sample(n, n)]
    x
  }


  folds <- plyr::ddply(.data = df, .variables = .(df$Class), .fun = plyr::here(createFolds), k = k_fold)
  folds_i <- folds$FOLDS

  # Prop1 = plyr::ddply(.data = folds,.variables = .(FOLDS) ,.fun =
  # plyr::here(plyr::summarise), prop = sum(.(Class))/length(.(Class)))

  Prop1 <- stats::aggregate(folds$Class, by = list(Category = folds$FOLDS), FUN = mean)

  index <- vector("list", k_fold)
  Xtrain <- vector("list", k_fold)
  Ytrain <- vector("list", k_fold)
  Xnew <- vector("list", k_fold)
  for (k in 1:k_fold) {
    index[[k]] <- which(folds_i == k)
    Xtrain[[k]] <- x[-index[[k]], ]
    Ytrain[[k]] <- y[-index[[k]]]
    Xnew[[k]] <- x[index[[k]], ]
  }

  RMSECV <- c()
  yres <- vector("list")
  for (j in 1:NumOrtho) {

    ropls.train <- mapply(OPLSDA, x = Xtrain, y = Ytrain, MoreArgs = list(impT = FALSE,
                          impG = FALSE, no = j, out.path = "."), SIMPLIFY = FALSE, USE.NAMES = TRUE)

    ypred <- mapply(OPLSDA_pred, ropls = ropls.train, x.new = Xnew, SIMPLIFY = FALSE,
                    USE.NAMES = TRUE)

    Ypred <- unlist(ypred)
    sorted <- order(names(Ypred))
    Ypred <- Ypred[order(names(Ypred))]

    ytrue <- df$Class
    names(ytrue) <- df$rowname
    ytrue <- ytrue[order(names(ytrue))]

    yres[[j]] <- cbind(ytrue = ytrue, Ypred = Ypred)

    RMSECV[(j)] <- sqrt(sum((ytrue - Ypred)^2)/n)

  }

  rcvOPLSDA <- list(RMSECV = RMSECV, folds_i = folds_i, Prop1 = Prop1, yres)


  if (ImpG == TRUE) {
    plot(RMSECV, xaxt = "n", ylab = "RMSECV", main = "Validation plot", type = "b",
         xlab = "Number of orthogonal components")
    axis(side = 1, at = c(1:NumOrtho), labels = c(1:NumOrtho))
  }
  return(rcvOPLSDA)

}



