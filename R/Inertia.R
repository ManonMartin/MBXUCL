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
#'   \item{\code{Name.group}}{Group membership}
#'   \item{\code{InertiaT}}{Total inertia}
#'   \item{\code{InertiaB}}{Between inertia}
#'   \item{\code{InertiaW}}{Within inertia}
#' }
#'
#' @examples
#' data("iris")
#' data = as.matrix((iris[,1:4]))
#' class = as.numeric(iris[,5])
#'
#' Res = Inertia(x = data, y = class, print = TRUE)
#'


Inertia = function(x, y, print=FALSE) {

# checks
if (! is.numeric(x)) {
  stop(deparse(substitute(x)), " is not numeric.")
}

if (! is.numeric(y)) {
  stop(deparse(substitute(y)), " is not numeric.")
}

if (! is.logical(print)) {
    stop(deparse(substitute(print)), " is not logical.")
}


# le calcul du barycentre de chaque groupe :
#############################################

nGroup=length(unique(y))

#Creation of group lists for data, y and barycenters
m=dim(x)[2]
n=dim(x)[1]

yG=c() # number of the class
nG=c() # number of elements for each group
xG=vector("list", nGroup) # data for elements of i's group
xmG=vector("list", nGroup) # barycenter for elements of i's group
j=1
for (i in unique(y)){
  yG[j]=i
  nG[j]=length(which(y==i))
  xG[[j]]=x[which(y==i),]
  xmG[[j]]=apply(xG[[j]],2,mean)
  j=j+1
}


# Raw mean
G=apply(x,2,mean)


#### Calculation of Between Inertia  (sum ng*d2(ug,u))
#############################################
IB=c()
for (i in 1:length(unique(y))){
  ib = nG[i]*(sum((xmG[[i]]-G)^2))
  IB=c(IB,ib)
  InertiaB=sum(IB)
}


#### Calculation of Total Inertia
#############################################
xmG=matrix(G,nrow=n,ncol=m,byrow=TRUE)
it=(x-xmG)^2
IT=apply(it,1,sum)


InertiaT = sum(IT)


iner.inter <- (InertiaB/InertiaT) # Between-group Inertia (% of tot inertia) :
iner.inter100=iner.inter*100


#### Calculation of Within Inertia
#############################################
iner.intra = 1-iner.inter
InertiaW=(iner.intra * InertiaT)



iner.intra100=iner.intra*100 # Within-group Inertia (% of tot inertia) :


PInertiaT = InertiaT*100/InertiaT
PInertiaW = InertiaW*100/InertiaT
PInertiaB = InertiaB*100/InertiaT


res1=matrix(data = c(InertiaB, InertiaW, InertiaT, PInertiaB, PInertiaW, PInertiaT), ncol = 3,
           dimnames = list(c("Value", "Percentage"), c("BI", "WI", "TI")), byrow = TRUE)



#### Calcul de l'inertie de chaque groupe
###############################################

Inertia_group=c() # inertie par groupe
Inertia_moy_group=c() # inertie moyenne par groupe


for (i in 1:nGroup){
  xjmoyen=matrix(xmG[[i]],nrow=nG[i],ncol=m,byrow=TRUE)
  xij = x[which(y==unique(y)[i]),]
  # inertie par groupe
  In_g = sum((xij-xjmoyen)^2) # Somme_i(xij –xjmoyen)2  
  Inertia_group=c(Inertia_group,In_g)
  # inertie moyenne par groupe
  In_m_g = (1/nG[i])*sum((xij-xjmoyen)^2) # (1/taille du groupe nj)*Somme_i(xij –xjmoyen)2
  Inertia_moy_group=c(Inertia_moy_group,In_m_g)
}

Inertia_TOT_group = sum(Inertia_group)

## pourcentage d'inertie
Inertia_group100 = 100*c(Inertia_group,Inertia_TOT_group)/Inertia_TOT_group # Somme_i(xij –xjmoyen)2  /Somme_ij(xij –xjmoyen)2 



res2=matrix(data = c(Inertia_group, Inertia_TOT_group,Inertia_group100, Inertia_moy_group, NA), ncol = 3,
           dimnames = list(c(paste0("Group ", unique(y)), "Total"), c("Inertia_group", "Inertia_group100", "Inertia_moy_group")), 
           byrow = FALSE)




if (print==TRUE) {
  cat("\n Total Inertia",InertiaT)
  cat("\n Intergroup Inertia",InertiaB)
  cat("\n Intergroup Inertia (% of Itot)",iner.inter100)
  cat("\n Within Inertia ",InertiaW)
  cat("\n Within Inertia (% of Itot)",iner.intra100, "\n")
  cat("\n Inertia per group ",Inertia_group, "\n")
  cat("\n Inertia per group (% of I_pergroup)",Inertia_group100, "\n")
  cat("\n Mean Inertia per group ",Inertia_moy_group, "\n")
}




return(list(Between_within = res1, Per_group = res2))


}
