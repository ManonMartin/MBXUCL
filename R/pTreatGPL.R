### PRE TRAITEMENTS GPL (GLOBAL PEAK LIST) #####################################

# Suppression of negative intensities
# Suppression of the water zone (square)
# Normalization (sum = 1)



#' @export pTreatGPL
#' @title 2D-NMR Global Peak List pre-treatments
#'
#' @description
#' Elimination of negative intensities (=0), Elimination of the water zone,
#' Normalization of the intensities (sum = 1)
#'
#' @param Intensities A numeric matrix of intensities (from individual GPL object).
#' @param Coordinates A numeric matrix of coordinates (C1, C2).
#' @param LowW Number corresponding to the lower value of the water zone to suppress (4.63).
#' @param UpW Number corresponding to the upper value of the water zone to suppress (5.1).
#'
#' @return A pre-treated matrix of intensities (with rows = spectra)
#'
#' @author Baptiste Feraud
#'
#' @examples
#' path <-  system.file("extdata", package = "MBXUCL")
#' GPL = read.table(file.path(path, "UrineGPL.csv"), sep=",",
#'  header=TRUE, row.names = 1, na.strings = "." )
#'
#' Int = GPL[,seq(3,dim(GPL)[2]- 1, 2)]
#' Coord = GPL[,c(1,2)]
#' IntT <- pTreatGPL(Int, Coord, 4.63, 5.1)
#' IntT[,1:2]
#'


pTreatGPL = function(Intensities, Coordinates, LowW, UpW){

# checks

if (missing(Intensities)){
warning("Intensities is missing with no default value")
}

if (missing(Coordinates)){
warning("Coordinates is missing with no default value")
}

if (missing(LowW)){
warning("LowW is missing with no default value")
}

if (missing(UpW)){
warning("UpW is missing with no default value")
}


if (! is.numeric(LowW)) {
  stop(deparse(substitute(LowW)), " is not numeric.")
}

if (! is.numeric(UpW)) {
  stop(deparse(substitute(UpW)), " is not numeric.")
}

if (LowW > UpW) {
  stop(deparse(substitute(UpW)), " the lower value of the water zone can't be greater than the upper value.")
}


#Negative intensities ==> zeros
Intensities[Intensities < 0] = 0

Mat_corr <- cbind(Coordinates, Intensities)

Mat_corr <- subset(Mat_corr, apply(Intensities, 1, sum) != 0)

#Suppression de la zone de l'eau
Mat_corr_w <- subset(Mat_corr, !((Mat_corr[,1] > LowW & Mat_corr[,1] < UpW) &
                                                      (Mat_corr[,2] > LowW & Mat_corr[,2] < UpW) ))

#Normalization, sum = 1
Intensities <- t(Mat_corr_w[-c(1,2)])
IntT <- Intensities / apply(Intensities, 1, sum)

return(IntT)
}

