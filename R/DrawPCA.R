#' @export DrawPCA
#' @title Scores and Loadings plots
#'
#' @description
#' Draws scores and Loadings plots for PCA analysis derived from the SVDforPCA function.
#'
#' @param PCAobj The objects resulting from a PCA analysis with the SVDforPCA function.
#' @param drawNames If \code{TRUE}, will show the observations names on the Scores plot.
#' @param createWindow If \code{TRUE}, will create a new window for the plot.
#' @param main Plot title. If \code{NULL}, default title is provided.
#' @param class Optional numeric vector giving the class of the observations.
#' @param axes Numerical vector indicating the PC axes that are drawn. Only the two first values are considered for Scores plot and only the first is considered for Loadings plot.
#' @param type The type of plot, either \code{"scores"} or \code{"loadings"}
#' @param loadingstype The type of Loadings plot, either a line plot (\code{"l"}) or points with histogram-like vertical lines (\code{"p"}).
#' @param num.stacked Number of stacked plots if \code{type} is \code{"loadings"}.
#' @param xlab Label of the x-axis.
#'
#' @return A score or loading plot in the current device.

#'
#' @examples
#'
#' data("HumanSerum")
#' PCAobj = SVDforPCA(HumanSerumSpectra)
#' class = ClassHS
#'
#' DrawPCA(PCAobj, drawNames=TRUE,
#' createWindow=FALSE, main = "PCA score plot for HumanSerum dataset",
#'   class = class, axes =c(1,2), type ="scores")
#'
#' DrawPCA(PCAobj, drawNames=TRUE,
#' createWindow=FALSE, main = "PCA loadings plot for HumanSerum dataset",
#'    axes = 1, type ="loadings", loadingstype="l")
#'
#' @importFrom grDevices dev.new
#' @import ggplot2
#' @import reshape2


DrawPCA <- function (PCAobj, drawNames=TRUE,
                           createWindow=FALSE, main = NULL, class = NULL, axes =c(1,2),
                           type =c("scores", "loadings"), loadingstype=c("l", "p"), num.stacked = 4, xlab = NULL) {

  checkArg(main, "str", can.be.null=TRUE)
  class = as.factor(class)
  loadingstype=match.arg(loadingstype)
  type = match.arg(type)

  m = dim(PCAobj$original.dataset)[1]
  n = dim(PCAobj$original.dataset)[2]

  # class

  if(is.vector(class, mode = "any") & length(class)!=m){
    stop("the length of class is not equal to the nrow of data matrix")
  }

  # axes
  if(!is.vector(axes, mode = "numeric")){
    stop("axes is not a numeric vector")
  }


  Xax = axes[1]
  Yax = axes[2]




  # Eigenvalues
  eig <- PCAobj$eigval

  # Variances in percentage
  variance <- eig*100/sum(eig)

  # scores
  scores = as.data.frame(PCAobj$pcs)

  # loadings
  loadings = PCAobj$pcv

  colnames(loadings) = paste0("Loading", c(1:dim(loadings)[2]))
  loadings = as.data.frame(loadings)

  plots <- list()
  Var = rowname = value = NULL # only for R CMD check

##########################################

if (type == "scores") {

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
    ggplot2::labs(x=paste0("PC",Xax," (", round(variance[Xax],2) ,"%)"), y=paste0("PC",Yax," (", round(variance[Yax],2) ,"%)"))

  if (drawNames) {
    if(is.null(class)) {
      plots = plots + ggplot2::geom_text(ggplot2::aes(x = scores[,Xax], y = scores[,Yax], label = rownames(PCAobj$original.dataset)),
                                         hjust = 0, nudge_x = (Xlim[2]/25),  show.legend = FALSE, size = 2)
    } else {plots = plots + ggplot2::geom_text(ggplot2::aes(x = scores[,Xax], y = scores[,Yax], label = rownames(PCAobj$original.dataset), colour = class, shape = class),
                                               hjust = 0, nudge_x = (Xlim[2]/25), show.legend = F, size = 2)}
  }

  print(ggplot2::last_plot())


  # if(is.null(class)) {
  #   graphics::plot(PCAobj$pcs[,Xax], PCAobj$pcs[,Yax],
  #        xlab=paste0("PC",Xax," (", round(variance[Xax],2) ,"%)"), xlim=Xlim,
  #        ylab=paste0("PC",Yax," (", round(variance[Yax],2) ,"%)"), ylim=Ylim,
  #        main=ifelse(is.null(main), "PCA Scores plot", main))
  #   graphics::abline(v=0, h=0, lty = 2)
  #   if (drawNames==TRUE) {
  #     graphics::text(PCAobj$pcs[,Xax], PCAobj$pcs[,Yax], rownames(PCAobj$original.dataset), pos=c(2,3))
  #   }
  # }else{
  #
  #   graphics::plot(PCAobj$pcs[,Xax], PCAobj$pcs[,Yax], col=class,
  #        xlab=paste0("PC",Xax," (", round(variance[Xax],2) ,"%)"), xlim=Xlim,
  #        ylab=paste0("PC",Yax," (", round(variance[Yax],2) ,"%)"), ylim=Ylim,
  #        main=ifelse(is.null(main), "PCA Scores plot", main))
  #   graphics::abline(v=0, h=0, lty = 2)
  #   if (drawNames==TRUE) {
  #     graphics::text(PCAobj$pcs[,Xax],PCAobj$pcs[,Yax],labels=rownames(PCAobj$original.dataset), pos=c(2,3), col=class)
  #   }
  # }



} else if (type == "loadings"){

  loadings = loadings[,axes]

  if (is.vector(loadings)) {
    n = 1
  }else {n = ncol(loadings)}


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


    plots <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var, y = value))
    if (loadingstype == "l") {
      plots = plots + ggplot2::geom_line()

    } else if (loadingstype == "p") {

      plots = plots + ggplot2::geom_point(size = 0.5)
    } else {warning("loadingstype is misspecified")}

    plots = plots + ggplot2::ggtitle(main) +
      ggplot2::facet_grid(rowname ~ ., scales = "free_y") +
      ggplot2::theme(legend.position="none") +
      ggplot2::labs(x=xlab, y = "Loadings") +
      ggplot2::geom_hline(yintercept = 0, size = 0.5, linetype = "dashed", colour = "gray60") +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = paste0("(",round(variance[i:last],2), "%)"), vjust=1, hjust=1)

    if ((melted[1,"Var"] - melted[(dim(melted)[1]),"Var"])>0) {
      plots =  plots + ggplot2::scale_x_reverse()
    }
    #
    #         require("gridExtra")
    i = last + 1
    print(ggplot2::last_plot())
  }


  # if (loadingstype == "l") {
  #   graphics::plot(loadings[,axes[1]], col="blue", xaxt="n", type="l", ylab = "Loading", xlab="Variables",
  #                  main=paste0(ifelse(is.null(main), "PCA Loading plot", main)),
  #                  sub =  paste0("\n PC", axes[1], " (", round(variance[axes[1]],2), "%)"))
  #   xat=round(as.numeric(rownames(loadings)),3)
  #   at=seq(1,length(xat), length(xat)/20)
  #   graphics::abline(h=0, lty=2)
  #   graphics::axis(side=1, at=at, labels=round(xat[at],2))
  # } else if (loadingstype == "p") {
  #   graphics::plot(loadings[,axes[1]], col="blue",  type="p", ylab = "Loading", xlab="Variables",
  #                  main=paste0(ifelse(is.null(main), "PCA Loading plot", main)),
  #                              sub =  paste0("\n PC", axes[1], " (", round(variance[axes[1]],2), "%)"))
  #   graphics::lines(loadings[,axes[1]], col="blue",  type="h")
  #   graphics::abline(h=0, lty=2)
  #   graphics::axis(side=1, at=1:dim(loadings)[1], labels=rownames(loadings))
  # } else warning("loadingstype is misspecified")

}

} # END
