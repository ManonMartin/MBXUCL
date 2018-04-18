#' @export DrawScores
#' @title Scores plots
#'
#' @description
#' Draws scores plots for the SVDforPCA, PLSDA or OPLSDA functions.
#'
#' @param obj The objects resulting from a PCA, PLSDA or OPLSDA analysis.
#' @param type.obj The type of object to be plotted.
#' @param drawNames If \code{TRUE}, will show the observations names on the scores plot.
#' @param createWindow If \code{TRUE}, will create a new window for the plot.
#' @param main Plot title. If \code{NULL}, default title is provided.
#' @param color Optional character, factor or numeric vector giving the color of the observations.
#' @param pch Optional character, factor or numeric vector giving the pch of the observations.
#' @param size The points size.
#' @param cex.lab The size of points labels.
#' @param axes Numerical vector indicating the PC axes that are drawn. Only the two first values are considered for scores plot. See details#' @param num.stacked Number of stacked plots if \code{type} is \code{'loadings'}.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param drawEllipses If \code{TRUE}, will draw ellipses with the \code{ggplot2::stat_ellipse} with groups coresponding to the color vector.
#' @param typeEl The type of ellipse, either "norm" (multivariate normal distribution), "t" (multivariate t-distribution) and "euclid" draws a circle with the radius equal to level, representing the euclidean distance from the center.
#' @param levelEl The confidence level at which to draw an ellipse.
#' @param drawPolygon If \code{TRUE}, will relate the points linked to the same group (color vector)
#' @param noLegend If \code{TRUE}, no legend is drawn.
#' @param legend_color_manual If not \code{NULL}, a named character vector giving manually the colors to be drawn, named with the levels of the color vector.
#' @param legend_shape_manual If not \code{NULL}, a named character vector giving manually the shapes to be drawn, named with the levels of the pch vector.
#' @return A score or loading plot in the current device.

#' @details
#' If \code{type.obj} is \code{'OPLSDA'}, axes = 1 represents the predictive score vector, axes = 2 represents the first orthogonal score vector, etc.
#'
#' @examples
#'
#' data('HumanSerum')
#' res.PCA = SVDforPCA(HumanSerumSpectra)
#' class = ClassHS
#'
#' DrawScores(res.PCA, drawNames=TRUE, type.obj = 'PCA',
#' createWindow=FALSE, main = 'PCA score plot for HumanSerum dataset',
#'   color = class, axes =c(1,2))
#'
#'
#' @importFrom grDevices dev.new
#' @import ggplot2
#' @import reshape2
#' @import gridExtra
#' @importFrom plyr ddply

DrawScores <- function(obj, type.obj = c("PCA", "PLSDA", "OPLSDA"), drawNames = TRUE,
                   createWindow = FALSE, main = NULL, color = NULL, pch = NULL, size = 1,
                   cex.lab = 3, axes = c(1, 2), xlab = NULL, ylab = NULL, drawEllipses = FALSE,
                   typeEl = "norm", levelEl = 0.9, drawPolygon = FALSE, noLegend = FALSE,
                   legend_color_manual = NULL, legend_shape_manual = NULL) {

  checkArg(main, "str", can.be.null = TRUE)

# color = numeric, or factor
# pch = numeric, or factor


  type.obj <- match.arg(type.obj)


  m <- dim(obj$original.dataset)[1]
  nn <- dim(obj$original.dataset)[2]



  # color

  if (!is.null(color) && is.vector(color, mode = "any") ) {
    if (length(color) != m) {
      stop("the length of color is not equal to the nrow of data matrix")
    }else if (!is.null(legend_color_manual) & is.vector(legend_color_manual, mode = "any") &
              length(legend_color_manual) != nlevels(as.factor(color))) {
    stop("The length of legend_color_manual is not equal to the nlevels of color")
    } else if (sum(!names(legend_color_manual) %in% levels(as.factor(color)))!=0) {
    stop("The names of legend_color_manual are not equal to the levels names of color")
    }
  }

  # pch
  if (!is.null(pch) && is.vector(pch, mode = "any")) {
    if (length(pch) != m) {
    stop("the length of pch is not equal to the nrow of data matrix")
    } else if (!is.null(legend_shape_manual) & is.vector(legend_shape_manual, mode = "any") &
               length(legend_shape_manual) != nlevels(as.factor(pch))) {
      stop("The length of legend_shape_manual is not equal to the nlevels of pch")
    } else if (sum(!names(legend_shape_manual) %in% levels(as.factor(pch)))!=0) {
      stop("The names of legend_shape_manual are not equal to the levels names of pch")
    }
  }





  if(drawPolygon & is.null(color)) {
    stop("drawPolygon is true but color is null ...")
  }



  # axes
  if (!is.vector(axes, mode = "numeric")) {
    stop("axes is not a numeric vector")
  }

  if(drawNames & is.null(rownames(obj$original.dataset))) {
    warning("rownames of the dataset is NULL, hence you cannot use DrawNames=TRUE ==> reset to FALSE")
    drawNames = FALSE
  }


  # Eigenvalues
  if (type.obj == "PCA") {
    eig <- obj$eigval
    # Variances in percentage
    variance <- eig * 100/sum(eig)
  }


  # scores
  Xax <- axes[1]
  Yax <- axes[2]

  if (type.obj == "PCA") {
    XaxName <- paste0("PC", Xax, " (", round(variance[Xax], 2),"%)")
    YaxName <- paste0("PC", Yax, " (", round(variance[Yax], 2), "%)")
  } else if (type.obj == "OPLSDA") {
    XaxName <- ifelse(Xax == 1, "Tp", paste0("To", Xax))
    YaxName <- ifelse(Yax == 1, "Tp", paste0("To", Yax - 1))
  } else { # PLS-DA
    XaxName <- paste0("Tp", Xax)
    YaxName <- paste0("Tp", Yax)
  }



  if (type.obj == "OPLSDA") {
    XaxName <- ifelse(Xax == 1, "Tp", paste0("To", Xax))
    YaxName <- ifelse(Yax == 1, "Tp", paste0("To", Yax - 1))

    obj$scores <- cbind(Tp = obj$Tp, obj$Tortho)
    colnames(obj$scores) <- c("Tp", paste0("To", 1:dim(obj$Tortho)[2]))
  }

  class(obj$scores) <- "numeric"
  scores <- as.data.frame(obj$scores)


  # define color and pch

  if (!is.null(color)) {
    color_factor <- as.factor(color)
    namecolor <- deparse(substitute(color))
  }

  if (!is.null(pch)) {
    pch_factor <- as.factor(pch)
    namepch <- deparse(substitute(pch))
  }


  plots <- list()
  plot <- list()
  Var <- rowname <- value <- NULL  # only for R CMD check

  ##########################################


    # labs
    if (is.null(xlab)) {
      xlab <- XaxName
    }
    if (is.null(ylab)) {
      ylab <- YaxName
    }

    if (createWindow)  {
      grDevices::dev.new(noRStudioGD = TRUE)
    }
    Xlim <- c(min(scores[, Xax]) * 1.4, max(scores[, Xax]) * 1.4)
    Ylim <- c(min(scores[, Yax]) * 1.4, max(scores[, Yax]) * 1.4)



    if (!drawEllipses & !drawPolygon) {
      plots <- ggplot2::ggplot(scores , ggplot2::aes(get(colnames(scores)[Xax]),
                 get(colnames(scores)[Yax])))  + ggplot2::xlim(Xlim) + ggplot2::ylim(Ylim)

    }else {
      plots <- ggplot2::ggplot(scores , ggplot2::aes(get(colnames(scores)[Xax]),get(colnames(scores)[Yax])))
    }

    if (is.null(color) & is.null(pch)) {
      # no color & no shape
      plots <- plots + ggplot2::geom_point(size=size)

    } else if (!is.null(color) & is.null(pch)) {
      # color
       plots <- plots  +
        ggplot2::geom_point(ggplot2::aes(colour = color_factor),size=size)
         if (!is.null(legend_color_manual)) {
           plots <- plots  +
             ggplot2::scale_colour_manual(name = namecolor, breaks = color_factor,
                                          values = legend_color_manual,
                                         guide=guide_legend(order=1))+
             guides(colour = guide_legend(override.aes = list(shape = 15)))
         } else {
           plots <- plots  +
           scale_colour_discrete(name = namecolor, breaks = unique(color_factor),
                                 labels = as.character(unique(color)),
                                 guide=guide_legend(order=1))+
             guides(colour = guide_legend(override.aes = list(shape = 15)))
         }


       if (drawPolygon) {
         dataf <- cbind(color_factor, scores)
         x <- colnames(scores)[Xax]
         y <- colnames(scores)[Yax]
         chulls <- plyr::ddply(dataf, .(color_factor), function(df) df[chull(df[,x], df[,y]), ])

         plots <- plots  +
           geom_polygon(data = chulls, aes(get(colnames(scores)[Xax]),
                                           get(colnames(scores)[Yax]), fill = color_factor, color = color_factor),
                        alpha=0.4,show.legend = FALSE, size =0.3)


       }

      if (drawEllipses) {

        plots <- plots + ggplot2::stat_ellipse(mapping = aes(get(colnames(scores)[Xax]),
                                                             get(colnames(scores)[Yax]),
                                               colour = color_factor),
                                               data = scores, type = typeEl,
                                               level = levelEl)
      }

    } else if (is.null(color) & !is.null(pch)) {
      # shape
      plots <- plots + ggplot2::geom_point(ggplot2::aes(shape = pch_factor), size=size)+
        scale_shape_manual(name = namepch, values=seq(0,26),
                           guide=guide_legend(order=1, shape = 1))


      # legend_shape_manual
      if (!is.null(legend_shape_manual)) {
        plots <- plots  +
          ggplot2::scale_shape_manual(name = namepch, breaks = pch_factor,
                                      values = legend_shape_manual,
                                       guide=guide_legend(order=1))
      }

    } else {
      # color + shape
      if (namecolor!=namepch) {
        plots <- plots + ggplot2::geom_point(ggplot2::aes(colour = color_factor, shape = pch_factor), size=size)+
          scale_shape_manual(values=seq(0,26))
      } else
      plots <- plots + ggplot2::geom_point(ggplot2::aes(colour = color_factor, shape = color_factor), size=size)+
          scale_shape_manual(name = namecolor, values=seq(0,26),
                             guide=guide_legend(order=1, shape = 1))

      # legend_color_manual
      if (!is.null(legend_color_manual)) {
        plots <- plots  +
          ggplot2::scale_colour_manual(name = namecolor, breaks = unique(color_factor),
                                       values = legend_color_manual,
                                       guide=guide_legend(order=1, shape = 1))+
        guides(colour = guide_legend(override.aes = list(shape = 15)))

      }else {
        plots <- plots  +
          scale_colour_discrete(name = namecolor,
                                guide=guide_legend(order=1, shape = 1))
         if (namecolor!=namepch) {
           plots <- plots  + guides(colour = guide_legend(override.aes = list(shape = 15)))
         }
      }

      # legend_shape_manual
      if (namecolor!=namepch) {
        if (!is.null(legend_shape_manual)) {
          plots <- plots  +
            ggplot2::scale_shape_manual(name = namepch, breaks = unique(pch_factor),
                                        values = legend_shape_manual,
                                        guide=guide_legend(order=1))
        }


      } else {

        if (!is.null(legend_shape_manual)) {
          plots <- plots  +
            ggplot2::scale_shape_manual(name = namecolor, breaks = unique(color_factor),
                                        values = legend_shape_manual,
                                        guide=guide_legend(order=1))
        }

      }


      if (drawPolygon) {
        dataf <- cbind(color_factor, scores)
        x <- colnames(scores)[Xax]
        y <- colnames(scores)[Yax]
        chulls <- plyr::ddply(dataf, .(color_factor), function(df) df[chull(df[,x], df[,y]), ])

        plots <- plots  +
        geom_polygon(data = chulls, aes(get(colnames(scores)[Xax]),
                                        get(colnames(scores)[Yax]), fill = color_factor, color = color_factor),
                     alpha=0.4,show.legend = FALSE, size =0.3)
      }


      if (drawEllipses) {

        plots <- plots + ggplot2::stat_ellipse(mapping = aes(get(colnames(scores)[Xax]),
                                                             get(colnames(scores)[Yax]),
                                               colour = color_factor),
                                               data = scores, type = typeEl,
                                               level = levelEl)
      }

    }



    plots <- plots + ggplot2::labs(title = main, x = xlab, y = ylab) + ggplot2::geom_vline(xintercept = 0,
                                                                                           size = 0.1) + ggplot2::geom_hline(yintercept = 0, size = 0.1) + ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray60",
                                                              size = 0.2), panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = "gray98"))


    if (drawNames) {

      if (is.null(color)) {
        plots <- plots + ggplot2::geom_text(ggplot2::aes(x = scores[, Xax],
                                                         y = scores[, Yax], label = rownames(obj$original.dataset)), hjust = 0,
                                            nudge_x = (Xlim[2]/25), show.legend = FALSE, size = cex.lab)
      } else {
        plots <- plots + ggplot2::geom_text(ggplot2::aes(x = scores[, Xax],
                                                         y = scores[, Yax], label = rownames(obj$original.dataset), colour = color_factor),
                                            hjust = 0, nudge_x = (Xlim[2]/25), show.legend = F, size = cex.lab)
      }
    }

    if (noLegend){
      plots <- plots + theme(legend.position="none")
    }

    plots



}  # END
