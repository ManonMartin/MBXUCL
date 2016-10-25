#' @export DrawSL
#' @title Scores and Loadings plots
#'
#' @description
#' Draws scores and Loadings plots for PCA analysis derived from the SVDforPCA function.
#'
#' @param obj The objects resulting from a PCA analysis with the SVDforPCA function.
#' @param type.obj The type of object to be plotted.
#' @param drawNames If \code{TRUE}, will show the observations names on the Scores plot.
#' @param createWindow If \code{TRUE}, will create a new window for the plot.
#' @param main Plot title. If \code{NULL}, default title is provided.
#' @param class Optional numeric vector giving the class of the observations.
#' @param axes Numerical vector indicating the PC axes that are drawn. Only the two first values are considered for scores plot.
#' @param type.graph The type of plot, either \code{"scores"} or \code{"loadings"}
#' @param loadingstype The type of Loadings plot, either a line plot (\code{"l"}) or points with histogram-like vertical lines (\code{"p"}).
#' @param num.stacked Number of stacked plots if \code{type} is \code{"loadings"}.
#' @param xlab Label of the x-axis.
#'
#' @return A score or loading plot in the current device.

#'
#' @examples
#'
#' data("HumanSerum")
#' res.PCA = SVDforPCA(HumanSerumSpectra)
#' class = ClassHS
#'
#' DrawSL(res.PCA, drawNames=TRUE, type.obj = "PCA",
#' createWindow=FALSE, main = "PCA score plot for HumanSerum dataset",
#'   class = class, axes =c(1,2), type.graph ="scores")
#'
#' DrawSL(res.PCA, drawNames=TRUE, type.obj = "PCA",
#' createWindow=FALSE, main = "PCA loadings plot for HumanSerum dataset",
#'    axes = 1, type.graph ="loadings", loadingstype="l")
#'
#' @importFrom grDevices dev.new
#' @import ggplot2
#' @import reshape2


DrawSL <- function (obj, type.obj = c("PCA", "PLSDA", "OPLSDA"), drawNames=TRUE,
                           createWindow=FALSE, main = NULL, class = NULL, axes =c(1,2),
                           type.graph =c("scores", "loadings"), loadingstype=c("l", "p"), num.stacked = 4, xlab = NULL) {

  checkArg(main, "str", can.be.null=TRUE)

   if (!is.null(class)){
    class = as.factor(class)
  }

  loadingstype=match.arg(loadingstype)
  type.graph = match.arg(type.graph)
  type.obj = match.arg(type.obj)


  m = dim(obj$original.dataset)[1]
  n = dim(obj$original.dataset)[2]

  # class

  if(!is.null(class) && is.vector(class, mode = "any") && length(class)!=m){
    stop("the length of class is not equal to the nrow of data matrix")
  }

  # axes
  if(!is.vector(axes, mode = "numeric")){
    stop("axes is not a numeric vector")
  }


  Xax = axes[1]
  Yax = axes[2]


  # Eigenvalues
  if (type.obj == "PCA")
  {eig <- obj$eigval
  # Variances in percentage
  variance <- eig*100/sum(eig)
  }


  # scores
  if (type.obj == "OPLSDA") {
    obj$scores = cbind(Tp = obj$Tp, obj$Tortho)
    colnames(obj$scores) = c("Tp", paste0("To", 1:dim(obj$Tortho)[2]))
  }

  class(obj$scores) = "numeric"
  scores = as.data.frame(obj$scores)

  # loadings
  if (type.obj == "OPLSDA") {
    obj$loadings = cbind(obj$Pp, obj$Portho)
    colnames(obj$loadings) = c("Pp", paste0("Po", 1:dim(obj$Portho)[2]))
    }
  class(obj$loadings) = "numeric"
  loadings = obj$loadings

  if (type.obj != "OPLSDA") {
  colnames(loadings) = paste0("Loading", c(1:dim(loadings)[2]))
  }

  loadings = as.data.frame(loadings)

  plots <- list()
  plot = list()
  Var = rowname = value = NULL # only for R CMD check

##########################################

if (type.graph == "scores") {

  if (createWindow) {
    grDevices::dev.new(noRStudioGD = TRUE)
  }
  Xlim=c(min(scores[,Xax])*1.4, max(scores[,Xax])*1.4)
  Ylim=c(min(scores[,Yax])*1.4, max(scores[,Yax])*1.4)

  plots <- ggplot2::ggplot(scores, ggplot2::aes(get(colnames(scores)[Xax]),get(colnames(scores)[Yax]))) +
    ggplot2::xlim(Xlim) +
    ggplot2::ylim(Ylim)

  if(is.null(class)) {
    plots <- plots + ggplot2::geom_jitter()
  } else {plots <- plots +  ggplot2::geom_jitter(ggplot2::aes(colour = class, shape = class))}

  plots <- plots + ggplot2::ggtitle(main) +
    ggplot2::geom_vline(xintercept = 0, size = 0.1) +
    ggplot2::geom_hline(yintercept = 0, size = 0.1) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray60", size = 0.2), panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = "gray98")) +
    if (type.obj == "PCA")
    {ggplot2::labs(x=paste0("PC",Xax," (", round(variance[Xax],2) ,"%)"), y=paste0("PC",Yax," (", round(variance[Yax],2) ,"%)"))
    } else {ggplot2::labs(x=paste0("Tp",Xax), y=paste0("Tp",Yax))}

   if (drawNames) {

    if(is.null(class)) {
      plots = plots + ggplot2::geom_text(ggplot2::aes(x = scores[,Xax], y = scores[,Yax], label = rownames(obj$original.dataset)),
                                         hjust = 0, nudge_x = (Xlim[2]/25),  show.legend = FALSE, size = 2)
    } else {plots = plots + ggplot2::geom_text(ggplot2::aes(x = scores[,Xax], y = scores[,Yax], label = rownames(obj$original.dataset), colour = class, shape = class),
                                               hjust = 0, nudge_x = (Xlim[2]/25), show.legend = F, size = 2)}
  }

  plots

} else {

  loadings = loadings[,axes]

  if (is.vector(loadings)) {
    n = 1
  }else {n = ncol(loadings)}

  j=1
  i = 1
  while (i <= n)
  {
    if (createWindow) {
      grDevices::dev.new(noRStudioGD = TRUE)
    }

    last = min(i + num.stacked-1, n)

    if (n == 1) {
      melted <- reshape2::melt(t(loadings), varnames=c("rowname", "Var"))
    }else {melted <- reshape2::melt(t(loadings[, i:last]), varnames=c("rowname", "Var"))}


    plot <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var, y = value))
    if (loadingstype == "l") {
      plot = plot + ggplot2::geom_line()
    } else {
      plot = plot + ggplot2::geom_point(size = 0.5)
    }

    plot = plot + ggplot2::ggtitle(main) +
      ggplot2::facet_grid(rowname ~ ., scales = "free_y") +
      ggplot2::theme(legend.position="none") +
      ggplot2::labs(x=xlab, y = "Loadings") +
      ggplot2::geom_hline(yintercept = 0, size = 0.5, linetype = "dashed", colour = "gray60")
      if (type.obj == "PCA"){
        plot = plot + ggplot2::annotate("text", x = -Inf, y = Inf, label = paste0("(",round(variance[i:last],2), "%)"), vjust=1, hjust=1)
      }

    if (is.numeric(melted[1,"Var"])){
      if ((melted[1,"Var"] - melted[(dim(melted)[1]),"Var"])>0) {
        plot =  plot + ggplot2::scale_x_reverse()
      }
    }

    #
    #         require("gridExtra")
    # plots
    plots[[j]] = plot
    i = last + 1
    j=j+1
  }

  plots
}

} # END
