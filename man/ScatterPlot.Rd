% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ScatterPlot.R
\name{ScatterPlot}
\alias{ScatterPlot}
\title{ScatterPlot}
\usage{
ScatterPlot(
  x,
  y,
  points_labs = NULL,
  createWindow = FALSE,
  main = NULL,
  color = NULL,
  pch = NULL,
  size = 1,
  cex.lab = 3,
  xlab = NULL,
  ylab = NULL,
  legend_pch = NULL,
  legend_color = NULL,
  drawEllipses = FALSE,
  typeEl = "norm",
  levelEl = 0.9
)
}
\arguments{
\item{x}{x-axis values.}

\item{y}{y-axis values.}

\item{points_labs}{If not \code{NULL}, will show the labels for each observation on the plot.}

\item{createWindow}{If \code{TRUE}, will create a new window for the plot.}

\item{main}{Plot title. If \code{NULL}, default title is provided.}

\item{color}{Optional character, factor or numeric vector giving the color of the observations. If \code{length(color)} = 1, the unique color is kept for all the points.}

\item{pch}{Optional character, factor or numeric vector giving the pch of the observations.}

\item{size}{The points size.}

\item{cex.lab}{The size of points labels.}

\item{xlab}{If not \code{NULL}, label for the x-axis.}

\item{ylab}{If not \code{NULL}, label for the y-axis.}

\item{legend_pch}{If not \code{NULL}, the labels for the pch legend.}

\item{legend_color}{If not \code{NULL}, the labels for the color legend.}

\item{drawEllipses}{If \code{TRUE}, will draw ellipses with the \code{ggplot2::stat_ellipse} with groups coresponding to the color vector.}

\item{typeEl}{The type of ellipse, either "norm" (multivariate normal distribution), "t" (multivariate t-distribution) and "euclid" draws a circle with the radius equal to level, representing the euclidean distance from the center.}

\item{levelEl}{The confidence level at which to draw an ellipse.}
}
\value{
A score or loading plot in the current device.
}
\description{
Draws scatter plots.
}
\details{
If \code{type.obj} is \code{'OPLSDA'}, axes = 1 represents the predictive score vector, axes = 2 represents the first orthogonal score vector, etc.
}
\examples{

x <- c(1:20)
y <- sample(20)
colors <- c(rep(1,10), rep(2,10))
points_pch <-  c(rep(17,10), rep(18,10))
labels <- as.character(c(1:20))
# automatic legend names and values assignment
ScatterPlot(x = x, y = y, points_labs = labels, createWindow=FALSE,
           main = 'Scatter plot', color = colors, pch = points_pch,
           xlab = "x axis", ylab = "y axis")

# decide on the color or pch values
legend_color = c(rep("A", 10), rep("B", 10))
legend_pch = c(rep("C", 10), rep("D", 10))
ScatterPlot(x = x, y = y, points_labs=labels, createWindow=FALSE,
           main = 'Scatter plot', color = colors, pch = points_pch,
           xlab = "x axis", ylab = "y axis",
           legend_pch = legend_pch,legend_color=legend_color)

}
