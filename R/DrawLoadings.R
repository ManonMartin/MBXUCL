#' @export DrawLoadings
#' @title Loadings plots
#'
#' @description
#' Draws Loadings plots for the SVDforPCA, PLSDA or OPLSDA functions.
#'
#' @param obj The objects resulting from a PCA, PLSDA or OPLSDA analysis.
#' @param type.obj The type of object to be plotted.
#' @param createWindow If \code{TRUE}, will create a new window for the plot.
#' @param main Plot title. If \code{NULL}, default title is provided.
#' @param axes Numerical vector indicating the PC axes that are drawn.
#' @param loadingstype The type of Loadings plot, either a line plot (\code{'l'}), points (\code{'p'}) or segments (\code{'s'}).
#' @param num.stacked Number of stacked plots.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param ang Angle to rotate the x axis labels for a better visualisation, either 0, 45 or 90 degrees.
#' @param xaxis_type Specify if the xaxis is numerical or character
#' @param nxaxis Number of thick marks on the xaxis for a character x variable
#' @param hline Numerical scalar. If not \code{NULL}, an horizontal dashed line is drawn at the \code{hline} value
#' @param xaxis_size Size of the x axis.
#'
#' @return A loading plot in the current device.
#'
#' @details
#' Better results printping is obtained by saving the plot as an object
#'
#' @examples
#'
#' data('HumanSerum')
#' res.PCA = SVDforPCA(HumanSerumSpectra)
#'
#' DrawLoadings(res.PCA, type.obj = 'PCA',
#' createWindow = FALSE, main = 'PCA loadings plot for HumanSerum dataset',
#'    axes = 1,  loadingstype='l', xlab = "L1")
#'
#' @importFrom grDevices dev.new
#' @import ggplot2
#' @import reshape2
#' @import gridExtra
#' @importFrom rlang .data



DrawLoadings <- function(obj, type.obj = c("PCA", "PLSDA", "OPLSDA"),
                         createWindow = FALSE, main = NULL,  axes = c(1, 2),
                         loadingstype = c("l", "p", "s"), num.stacked = 4, xlab = NULL, ylab = NULL,
                         ang = c("0", "45", "90"), xaxis_type = c("numerical", "character"),
                         nxaxis = 10, hline = NULL, xaxis_size = 11) {

  checkArg(main, "str", can.be.null = TRUE)
  checkArg(nxaxis, "num", can.be.null = FALSE)
  checkArg(hline, "num", can.be.null = TRUE)




  loadingstype <- match.arg(loadingstype)

  type.obj <- match.arg(type.obj)
  xaxis_type <- match.arg(xaxis_type)

  ang <- match.arg(ang)

  if (ang %in% c("0", "45")) {
    vjust <- 1
    hjust <- 0.5
  } else {
    vjust <- 0.5
    hjust <- 1
  }


  m <- dim(obj$original.dataset)[1]
  nn <- dim(obj$original.dataset)[2]



  # axes
  if (!is.vector(axes, mode = "numeric")) {
    stop("axes is not a numeric vector")
  }


  # Eigenvalues
  if (type.obj == "PCA") {
    eig <- obj$eigval
    # Variances in percentage
    variance <- eig * 100/sum(eig)
  }

  # loadings
  if (type.obj == "OPLSDA") {
    obj$loadings <- cbind(obj$Pp, obj$Portho)
    colnames(obj$loadings) <- c("Pp", paste0("Po", 1:dim(obj$Portho)[2]))
  }
  class(obj$loadings) <- "numeric"
  loadings <- obj$loadings

  if (type.obj != "OPLSDA"){
    colnames(loadings) <- paste0("Loading", c(1:dim(loadings)[2]))
  }

  loadings <- as.data.frame(loadings)

  plots <- list()
  plot <- list()
  Var <- rowname <- value <- NULL  # only for R CMD check

  ##########################################

  if (!is.null(ylab) && length(ylab) != 1) {
    warning("Since ylab is longer than 1, the first element is selected")
    ylab <- ylab[1]
  }

  if (!is.null(xlab) && length(xlab) != 1) {
    warning("Since xlab is longer than 1, the first element is selected")
    xlab <- xlab[1]
  }

  # loadings <- loadings[, axes]

  # labs
  if (is.null(xlab)) {
    xlab <- "Index"
  }

  if (is.null(ylab)) {
    ylab <- "Value"
  }


  n <- length(axes)


  j <- 1
  i <- 1
  while (i <= n) {

    last <- min(i + num.stacked - 1, n)

    if(i<last){
      load <- t(loadings[, i:last])
    }else{
      load <- as.matrix(t(loadings[, i]))
      names(load) <- rownames(loadings)
    }
      melted <- reshape2::melt(load, varnames = c("rowname", "Var"))
      if (n==1) {
        melted[,"rowname"] <-  rep(colnames(loadings)[axes],nn)
        if (xaxis_type == "numerical"){
          melted[,"Var"] <-  round(as.numeric(rownames(loadings)),3)
        } else{melted[,"Var"] <-  rownames(loadings)}
      }

      # to silence CMD check
      # in the future you will be able to use ggplot2::vars(.data$facet)
      # but this is currently not supported (tidyverse/ggplot2#2963)
      facet <- NULL; rm(facet)

      #######

  if (xaxis_type == "numerical") {
    if (type.obj == "PCA") {
      ylabs <- paste0("L", axes[i:last]," ",
                      paste0("(", round(variance[i:last],2), "%)"))
    }else {
      ylabs <- paste0("L", i:last)
    }

    plot_data <- data.frame(
      facet = ylabs,
      x = melted$Var,
      y = melted$value
    )
    plot <- ggplot2::ggplot(data = plot_data, ggplot2::aes(x = .data$x, y = .data$y))

      # plot <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var, y = value))
    } else  {
      melted$index <- as.numeric(as.factor(melted$Var))
      melted <- as.data.frame(melted)

      if (type.obj == "PCA") {
        ylabs <- paste0("L", i:last," ",
                        paste0("(", round(variance[i:last],2), "%)"))
      }else {
        ylabs <- paste0("L", i:last)
      }

      plot_data <- data.frame(
        facet = ylabs,
        x = melted$index,
        y = melted$value
      )

      plot <- ggplot2::ggplot(data = plot_data, ggplot2::aes(x = .data$x, y = .data$y))

    }

      plot <- plot + ggplot2::theme_bw()

    if (loadingstype == "p") {
      plot <- plot + ggplot2::geom_point(size=0.5)
      if (xaxis_type == "character"){
        plot <- plot + ggplot2::scale_x_continuous(breaks = seq(1, nn, floor(nn/nxaxis)),
                                                   labels = rownames(loadings)[seq(1, nn, floor(nn/nxaxis))])
      }
    } else if (loadingstype == "l")  {
      plot <- plot + ggplot2::geom_line()
      if (xaxis_type == "character"){
        plot <- plot + ggplot2::scale_x_continuous(breaks = seq(1, nn, floor(nn/nxaxis)),
                                                   labels = rownames(loadings)[seq(1, nn, floor(nn/nxaxis))])
      }
    } else  {

      plot <- plot + ggplot2::geom_hline(yintercept = 0, size = 0.5, linetype = "dashed", colour = "gray60")

        if (xaxis_type == "numerical"){
          plot <- plot + ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0),
                                               size = 0.5, lineend = "round")
        } else {
          plot <- plot + ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0),
                                               size = 0.5, lineend = "round")

          plot <- plot + ggplot2::scale_x_continuous(breaks = seq(1, nn, floor(nn/nxaxis)),
                                        labels = rownames(loadings)[seq(1, nn, floor(nn/nxaxis))])

        }

      }


    plot <- plot + ggplot2::labs(title = main, x = xlab, y = ylab) +
       ggplot2::facet_grid(facet ~  ., scales = "free_y") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = as.numeric(ang), vjust = vjust, hjust = hjust,
                                                         size = xaxis_size)) +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 90)) +
      ggplot2::theme(legend.position = "none")
      if (!is.null(hline)){
        plot <- plot + ggplot2::geom_hline(yintercept = hline, size = 0.5, linetype = "dashed", colour = "gray60")
      }

#
# #######
#     if (type.obj == "PCA") {
#       if (xaxis_type == "numerical")  {
#         plot <- plot + ggplot2::annotate("text", x = -Inf, y = Inf,
#                                          label = paste0("(",round(variance[i:last], 2), "%)"), vjust = 2, hjust = 1.5)
#       } else {
#         plot <- plot + ggplot2::annotate("text", x = 1, y = max(loadings),
#                                          label = paste0("(",round(variance[i:last], 2), "%)"), hjust = 0)
#       }
#
#       }
#
# #########


    if (xaxis_type == "numerical")  {
      if ((melted[1, "Var"] - melted[(dim(melted)[1]), "Var"]) > 0) {
        plot <- plot + ggplot2::scale_x_reverse()
      }
    }

    plots[[j]] <- plot
    i <- last + 1

    if (createWindow)  {
      grDevices::dev.new(noRStudioGD = TRUE)
    }

    # gridExtra::grid.arrange(plot)
    j <- j + 1

  }

  return(plots)

}  # END
