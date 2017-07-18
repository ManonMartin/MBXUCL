#' @export LinePlot
#' @title Line plots
#'
#' @description
#' Draws Line plots.
#'
#' @param X A numerical matrix containing the rows to be drawn.
#' @param createWindow If \code{TRUE}, will create a new window for the plot.
#' @param main Plot title. If \code{NULL}, default title is provided.
#' @param rows Numerical vector indicating the X matrix rows that are drawn.
#' @param type The type of plot, either a line plot (\code{'l'}), points (\code{'p'}) or segments (\code{'s'}).
#' @param num.stacked Number of stacked plots.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param ang Angle to rotate the x axis labels for a better visualisation, either 0, 45 or 90 degrees.
#' @param xaxis_type Specify if the xaxis is numerical or character (corresponds to the colnames of X).
#' @param nxaxis Number of thick marks on the xaxis for a character x variable.
#'
#' @return A loading plot in the current device.
#'
#' @details
#' Better results printping is obtained by saving the plot as an object
#'
#' @examples
#'
#' X <- rbind(sample(1:20,20), sample(21:40,20), sample(41:60,20))
#' colnames(X) <- as.character(1:20)
#'
#' a <- LinePlot(X = X, createWindow = FALSE, main = "line plot",  rows = c(1, 2),
#'              type = "l", num.stacked = 4, xlab = "x-axis", ylab = "y-axis",
#'              ang = "0", xaxis_type = "numerical", nxaxis = 10)
#'
#' @importFrom grDevices dev.new
#' @import ggplot2
#' @import reshape2
#' @import gridExtra



LinePlot <- function(X, createWindow = FALSE, main = NULL,  rows=NULL,
                         type = c("l", "p", "s"), num.stacked = 4, xlab = NULL, ylab = NULL,
                         ang = c("0", "45", "90"), xaxis_type = c("numerical", "character"), nxaxis = 10) {

  checkArg(main, "str", can.be.null = TRUE)
  checkArg(nxaxis, "num", can.be.null = FALSE)

  if (is.vector(X)){
    if (xaxis_type=="numerical"){
      X <- matrix(X, nrow = 1, dimnames = list(deparse(substitute(X)), 1:nn))
    }else{
      X <- matrix(X, nrow = 1, dimnames = list(deparse(substitute(X)), paste0("V",1:nn)))
    }
  }

  type <- match.arg(type)

  xaxis_type <- match.arg(xaxis_type)

  ang <- match.arg(ang)

  if (ang %in% c("0", "45")) {
    vjust <- 1
    hjust <- 0.5
  } else {
    vjust <- 0.5
    hjust <- 1
  }


  m <- dim(X)[1]
  nn <- dim(X)[2]

  if (is.null(rownames(X))){
    if (m==1){
      rownames(X) = deparse(substitute(X))
      }else{rownames(X) <- 1:m}
  }

  if (is.null(colnames(X))) {
    if (xaxis_type=="numerical"){
    colnames(X) <- 1:nn
    }else{
      colnames(X) <- paste0("V",1:nn)
    }
  }

    X <- as.data.frame(X)

  plots <- list()
  plot <- list()
  Var <- colname <- value <- NULL  # only for R CMD check

  ##########################################


  # labs
  if (is.null(rows)){
    n <- m
  } else {n <- length(rows)} # number of line plots to draw


  j <- 1
  i <- 1
  while (i <= n) {


    last <- min(i + num.stacked - 1, n)

    melted <- reshape2::melt(t(X[i:last, ]), varnames = c("Var", "colname"))
    if (n==1){
      if (!is.null(rows)) {
        melted[,"colname"] <-  rep(row.names(X)[rows],nn)
      } else { melted[,"colname"] <-  rep(row.names(X),nn)}
    }



    if (xaxis_type == "numerical") {
      plot <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var, y = value))
    }else{
      melted$index <- as.numeric(as.factor(melted$Var))
      melted <- as.data.frame(melted)
      plot <- ggplot2::ggplot(data = melted, ggplot2::aes(x = index, y = value))
    }

    plot <- plot + ggplot2::theme_bw()
    if (type == "p") {
      plot <- plot + ggplot2::geom_point(size=0.5)
      if (xaxis_type == "character"){
        plot <- plot + ggplot2::scale_x_continuous(breaks = seq(1, nn, floor(nn/nxaxis)),
                                                   labels = colnames(X)[seq(1, nn, floor(nn/nxaxis))])
      }
    } else if (type == "l")  {
      plot <- plot + ggplot2::geom_line()
      if (xaxis_type == "character"){
        plot <- plot + ggplot2::scale_x_continuous(breaks = seq(1, nn, floor(nn/nxaxis)),
                                                   labels = colnames(X)[seq(1, nn, floor(nn/nxaxis))])
      }
    } else  {
      if (xaxis_type == "numerical"){
        plot <- plot + ggplot2::geom_segment(ggplot2::aes(xend = Var, yend = 0),
                                             size = 0.5, lineend = "round")
      } else {
        plot <- plot + ggplot2::geom_segment(ggplot2::aes(xend = index, yend = 0),
                                             size = 0.5, lineend = "round")

        plot <- plot + ggplot2::scale_x_continuous(breaks = seq(1, nn, floor(nn/nxaxis)),
                                                   labels = colnames(X)[seq(1, nn, floor(nn/nxaxis))])

      }

    }


    plot <- plot + ggplot2::labs(title = main, x = xlab, y = ylab) + ggplot2::facet_grid(colname ~., scales = "free_y") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = as.numeric(ang), vjust = vjust, hjust = hjust)) +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 90)) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_hline(yintercept = 0, size = 0.5, linetype = "dashed", colour = "gray60")


    if (is.numeric(melted[1, "Var"]))  {
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
