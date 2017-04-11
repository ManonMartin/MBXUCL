#' @export Inertia
#' @title Measure of inertia
#'
#' @description
#' Calculates the within, between and total inertia of a data matrix
#'
#' @param x A numeric matrix.
#' @param y A numeric vector of data clusters.
#' @param print If \code{TRUE}, prints the inertia values
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{Between_within}}{A matrix with the following elements: Name.group = Group membership,
#'   InertiaT = Total inertia, InertiaB = Between inertia, InertiaW = Within inertia.}
#'   \item{\code{Per_group}}{A matrix with the following elements: Inertia_group = Inertia per group,
#'   Inertia_group100 =  Inertia per group (percentage of ITOT pergroup), Inertia_moy_group = Mean Inertia per group.}
#' }
#'
#' @examples
#' data('HumanSerum')
#' Res = Inertia(x = HumanSerumSpectra, y = ClassHS, print = TRUE)
#'


Inertia=function (x, y, print = FALSE)
{
  if (!is.numeric(x)) {
    stop(deparse(substitute(x)), " is not numeric.")
  }
  if (!is.numeric(y)) {
    stop(deparse(substitute(y)), " is not numeric.")
  }
  if (!is.logical(print)) {
    stop(deparse(substitute(print)), " is not logical.")
  }
  # Initialize matrices
  nGroup <- length(unique(y))
  m <- dim(x)[2]
  n <- dim(x)[1]
  yG <- c()
  nG <- c()
  xG <- vector("list", nGroup)
  xmG <- vector("list", nGroup)
  # Calculate group and ground means
  j <- 1
  for (i in unique(y)) {
    yG[j] <- i
    nG[j] <- length(which(y == i))
    xG[[j]] <- x[which(y == i), ]
    xmG[[j]] <- apply(xG[[j]], 2, mean)
    j <- j + 1
  }
  G <- apply(x, 2, mean)
  # Computes beween groups inertia
  IB <- c()
  for (i in 1:length(unique(y))) {
    ib <- nG[i] * (sum((xmG[[i]] - G)^2))
    IB <- c(IB, ib)
  }
  InertiaB <- sum(IB)
  # Computes total inertia
  MMG <- matrix(G, nrow = n, ncol = m, byrow = TRUE)
  IT <- (x - MMG)^2
  InertiaT <- sum(IT)
  # Compute all global results
  iner.inter <- (InertiaB/InertiaT)
  iner.inter100 <- iner.inter * 100
  iner.intra <- 1 - iner.inter
  iner.intra100 <- iner.intra * 100
  InertiaW <- (iner.intra * InertiaT)
  PInertiaT <- InertiaT * 100/InertiaT
  PInertiaW <- InertiaW * 100/InertiaT
  PInertiaB <- InertiaB * 100/InertiaT
  res1 <- matrix(data = c(InertiaB, InertiaW, InertiaT, PInertiaB,
                          PInertiaW, PInertiaT), ncol = 3, dimnames = list(c("Value",
                                                                             "Percentage"), c("BI", "WI", "TI")), byrow = TRUE)
  # Computes inertia per group
  Inertia_group <- c()
  Inertia_moy_group <- c()
  for (i in 1:nGroup) {
    xjmoyen <- matrix(xmG[[i]], nrow = nG[i], ncol = m, byrow = TRUE)
    In_g <- sum((xG[[i]] - xjmoyen)^2)
    Inertia_group <- c(Inertia_group, In_g)
    In_m_g <- (In_g/nG[i])
    Inertia_moy_group <- c(Inertia_moy_group, In_m_g)
  }
  Inertia_TOT_group <- sum(Inertia_group)
  Inertia_group100 <- 100 * c(Inertia_group, Inertia_TOT_group)/Inertia_TOT_group
  numbers <- c(nG, n)
  res2 <- matrix(data = c(numbers, c(Inertia_group, Inertia_TOT_group),
                          Inertia_group100, c(Inertia_moy_group, NA)), ncol = 4, dimnames = list(c(paste0("Group ",
                                                                                                          unique(y)), "Total"), c("N", "Inertia_group", "Inertia_group100",
                                                                                                                                  "Inertia_moy_group")),byrow=FALSE)
  if (print == TRUE) {
    cat("\n Total Inertia", InertiaT)
    cat("\n Intergroup Inertia", InertiaB)
    cat("\n Intergroup Inertia (% of Itot)", iner.inter100)
    cat("\n Within Inertia ", InertiaW)
    cat("\n Within Inertia (% of Itot)", iner.intra100, "\n")
    cat("\n Inertia per group ", Inertia_group, "\n")
    cat("\n Inertia per group (% of I_pergroup)", Inertia_group100,
        "\n")
    cat("\n Mean Inertia per group ", Inertia_moy_group,
        "\n")
  }
  return(list(Between_within = res1, Per_group = res2))
}
