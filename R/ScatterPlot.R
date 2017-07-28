#' @export ScatterPlot
#' @title ScatterPlot
#'
#' @description
#' Draws scatter plots.
#'
#' @param x x-axis values.
#' @param y y-axis values.
#' @param points_labs If not \code{NULL}, will show the labels for each observation on the plot.
#' @param createWindow If \code{TRUE}, will create a new window for the plot.
#' @param main Plot title. If \code{NULL}, default title is provided.
#' @param color Optional character, factor or numeric vector giving the color of the observations.
#' @param pch Optional character, factor or numeric vector giving the pch of the observations.
#' @param size The points size.
#' @param cex.lab The size of points labels.
#' @param xlab If not \code{NULL}, label for the x-axis.
#' @param ylab If not \code{NULL}, label for the y-axis.
#' @param legend_pch If not \code{NULL}, the labels for the pch legend.
#' @param legend_color If not \code{NULL}, the labels for the color legend.
#' @param drawEllipses If \code{TRUE}, will draw ellipses with the \code{ggplot2::stat_ellipse} with groups coresponding to the color vector.
#' @param typeEl The type of ellipse, either "norm" (multivariate normal distribution), "t" (multivariate t-distribution) and "euclid" draws a circle with the radius equal to level, representing the euclidean distance from the center.
#' @param levelEl The confidence level at which to draw an ellipse.
#'
#' @return A score or loading plot in the current device.

#' @details
#' If \code{type.obj} is \code{'OPLSDA'}, axes = 1 represents the predictive score vector, axes = 2 represents the first orthogonal score vector, etc.
#'
#' @examples
#'
#'x <- c(1:20)
#'y <- sample(20)
#'colors <- c(rep(1,10), rep(2,10))
#'points_pch <-  c(rep(17,10), rep(18,10))
#'labels <- as.character(c(1:20))
#'# automatic legend names and values assignment
#'ScatterPlot(x = x, y = y, points_labs = labels, createWindow=FALSE,
#'            main = 'Scatter plot', color = colors, pch = points_pch,
#'            xlab = "x axis", ylab = "y axis")
#'
#'# decide on the color or pch values
# while adding a personalized legend
#'legend_color = c(rep("A", 10), rep("B", 10))
#'legend_pch = c(rep("C", 10), rep("D", 10))
#'ScatterPlot(x = x, y = y, points_labs=labels, createWindow=FALSE,
#'            main = 'Scatter plot', color = colors, pch = points_pch,
#'            xlab = "x axis", ylab = "y axis",
#'            legend_pch = legend_pch,legend_color=legend_color)
#'
#' @importFrom grDevices dev.new
#' @import ggplot2
#' @import reshape2
#' @import gridExtra

ScatterPlot <- function(x, y, points_labs = NULL, createWindow = FALSE, main = NULL,
                        color = NULL, pch = NULL, size = 1, cex.lab = 3, xlab = NULL, ylab = NULL, legend_pch = NULL,
                        legend_color = NULL, drawEllipses = FALSE, typeEl = "norm", levelEl = 0.9) {


 checkArg(main, "str", can.be.null = TRUE)

 checkArg(x, "num", can.be.null = FALSE)
 checkArg(y, "num", can.be.null = FALSE)

 data <- cbind(x,y)
 data <- as.data.frame(data)

  # color = numeric, or factor
  # pch = numeric, or factor

  m <- length(x)
  # nn <- dim(obj$original.dataset)[2]

  # checks points_labs
  if (!is.null(points_labs) && is.vector(points_labs, mode = "any") && length(points_labs) != m) {
    stop("the length of color is not equal to the nrow of data matrix")
  }

  # checks color
  if (!is.null(color) && is.vector(color, mode = "any") && length(color) != m) {
    stop("the length of color is not equal to the nrow of data matrix")
  }


  # checks pch
  if (!is.null(pch) && is.vector(pch, mode = "any") && length(pch) != m) {
    stop("the length of pch is not equal to the nrow of data matrix")
  }


  # define color and pch

  if (!is.null(color)) {
    color_factor <- as.factor(color)
    namecolor <- deparse(substitute(color))
    if (is.null(legend_color)) {
      valuescolor <- color_factor
    } else {
        if (!sum(diag(table(color_factor, legend_color)))==m & !sum(diag(table(color_factor, legend_color)))==0) {
          warning("pch_factor and legend_pch do not have concordant classes")
          print(table(color_factor, legend_color))
        }
      valuescolor <- legend_color
      }
  }

  if (!is.null(pch)) {
    pch_factor <- as.factor(pch)
    namepch <- deparse(substitute(pch))
    if (is.null(legend_pch)) {
      valuespch <- pch_factor
    } else {
        if (!sum(diag(table(pch_factor, legend_pch)))==m & !sum(diag(table(pch_factor, legend_pch)))==0) {
          warning("pch_factor and legend_pch do not have concordant classes")
          print(table(pch_factor, legend_pch))
        }
      valuespch <- legend_pch
      }
  }


  plots <- list()
  plot <- list()
  Var <- rowname <- value <- NULL  # only for R CMD check

  ##########################################


  # labs
  if (is.null(xlab)) {
    xlab <- deparse(substitute(x))
  }
  if (is.null(ylab)) {
    ylab <- deparse(substitute(y))
  }

  if (createWindow)  {
    grDevices::dev.new(noRStudioGD = TRUE)
  }



  plots <- ggplot2::ggplot(data=data, ggplot2::aes(x=x,y=y))

  if (is.null(color) & is.null(pch)) {
    # no color & no shape
    plots <- plots + ggplot2::geom_point(size=size)

  } else if (!is.null(color) & is.null(pch)) {
    # color
    plots <- plots + ggplot2::geom_point(ggplot2::aes(colour = color_factor), size=size) +
      scale_colour_discrete(name = namecolor, breaks = color_factor,
                            labels = valuescolor,
                            guide=guide_legend(order=1))

    if (drawEllipses) {
      plots <- plots + ggplot2::stat_ellipse(mapping = aes(x=x,y=y,
                                             colour = color_factor),
                                             data = data, type = typeEl,
                                             level = levelEl)
    }


  } else if (is.null(color) & !is.null(pch)) {
    # shape
    plots <- plots + ggplot2::geom_point(ggplot2::aes(shape = pch_factor), size=size) +
      scale_shape_discrete(name = namepch, breaks = pch_factor,
                           labels = valuespch,
                           guide=guide_legend(order=1))
  } else {
    # color + shape
    plots <- plots + ggplot2::geom_point(ggplot2::aes(colour = color_factor, shape = pch_factor), size=size) +
      scale_colour_discrete(name = namecolor, breaks = color_factor,
                            labels = valuescolor,
                            guide=guide_legend(order=1)) +
      scale_shape_discrete(name = namepch, breaks = pch_factor,
                           labels = valuespch,
                           guide=guide_legend(order=2))

    if (drawEllipses) {
      plots <- plots + ggplot2::stat_ellipse(mapping = aes(x=x,y=y,
                                             colour = color_factor),
                                             data = data, type = typeEl,
                                             level = levelEl)
    }
  }


# title, theme
  plots <- plots + ggplot2::labs(title = main, x = xlab, y = ylab) + ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray60",
                                                            size = 0.2), panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = "gray98"))

# points_labs
  if (!is.null(points_labs)) {

    if (is.null(color)) {
      plots <- plots + ggplot2::geom_text(ggplot2::aes(x = x,
                                                       y = y, label = points_labs),
                                          hjust = 0, vjust = 1, show.legend = F, size = cex.lab)
    } else {
      plots <- plots + ggplot2::geom_text(ggplot2::aes(x = x,
                                                       y = y, label = points_labs, colour = color_factor),
                                          hjust = 0, vjust = 1,  show.legend = F, size = cex.lab)
    }
  }

  plots



}  # END
