\docType{methods}
\name{plotGlobalFit}
\alias{plotGlobalFit}
\alias{plotGlobalFit,Appac,numeric,list-method}
\title{
Plot the global (APPAC) fit function
}
\description{
S4 method to plot an APPAC global fit.
\loadmathjax{}
}
\usage{
\S4method{plotGlobalFit}{Appac,numeric,list}(object, size, colors)
%plotGlobalFit(object, colors)
}
\arguments{
  \item{object}{
an object of class \code{\link{Appac}}
}
  \item{size}{
a numeric value which controls the size of the axis labels and titles
}
  \item{colors}{
    a list containing the colors of the color scheme to be used in the plot:
      \describe{
%      \itemize{
        \item{highlight_color: }{the color in which the datapoints are represented}
        \item{lowlight_color: }{the color in which the errorbars are represented}
        \item{line_color: }{the color of the fitted function line}
        \item{fill_color: }{the color in which the standard uncertainties of the fitted function line are plotted}
%      }
    }
  }
}
\details{
APPAC fits a function of the form
\mjeqn{\hat{Y} = Y_{ref} \\times{(1 + \\kappa (P - P_{ref}))}}{area = area.ref * (1 + kappa * (P - P.ref))}
to a given data set. \mjeqn{P_{ref}}{P.ref} is set by the user; the fitted
area at the reference pressure \mjeqn{Y_{ref}}{Y.ref} is the expected value of
the area of a component in a sample at \mjeqn{P_{ref}}{P.ref}.
\mjeqn{Y_{ref}}{Y.ref} is obtained by local fits of area vs. P.
\code{appac.plot} plots the slope vs the area.ref of the local fits;
\mjeqn{\\kappa}{kappa} is given by the slope of the fitted line.
}
\value{
a ggplot object
}

\seealso{
\code{\link{appac}}
}
\examples{
## Please be patient, the calculation takes a while.

## Not run:
library(appac)

data("PGC_TCD")

## define P.ref as the median of the measured ambient pressures at injection time
## data set contains some NA, because of missing peaks
P.ref <- median(PGC_TCD$air.pressure, na.rm = TRUE)

## define the column names of the dataset
appac.colnames <- list(
  sample.col    = "sample.name",
  peak.col      = "peak.name",
  date.col      = "injection.date",
  pressure.col  = "air.pressure",
  area.col      = "raw.area"
)

appac.control <- list(min.data.points = 5, drift.model = "linear")

## appac will throw warnings that sample and peak names have been changed
x <- suppressWarnings(
       appac(data = PGC_TCD, 
       P.ref = P.ref, 
       appac.control = appac.control, 
       appac.colnames = appac.colnames)
     )

## View results
colors <- list(highlight_color="black", lowlight_color="darkgrey", 
               line_color="darkblue", fill_color="lightblue")
plotGlobalFit(x, size = 4, colors)
}

\keyword{plots}
