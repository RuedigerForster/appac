\docType{methods}
\name{plotGlobalFit}
\alias{plotGlobalFit}
\alias{plotGlobalFit,Appac-method}
\title{
Plot the global (APPAC) fit function
}
\description{
S4 method to plot an APPAC global fit.
\loadmathjax{}
}
\usage{
\S4method{plotGlobalFit}{Appac}(object, size = 4, ...)
%plotGlobalFit(object, colors)
}
\arguments{
  \item{object}{
an object of class \code{\link{Appac}}
}
  \item{size}{
a numeric value which controls the size of the axis labels and titles
}
  \item{\dots}{
    additional arguments:

    \code{colors}: a named list containing the color scheme to be used in the plot:
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
\mjeqn{\hat{Y} = Y_{ref} \\times{(1 + \\kappa (P - P_{ref}))}}{area = area.ref * (1 + kappa * (P - P_ref))}
to a given data set. \mjeqn{P_{ref}}{P_ref} is set by the user; the fitted
area at the reference pressure \mjeqn{Y_{ref}}{Y.ref} is the expected value of
the area of a component in a sample at \mjeqn{P_{ref}}{P_ref}.
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
## Not run:

library(dplyr)

## PGC_TCD is a short dataset of weekly measurements of 1 control sample taken on
## process analyzer. The data frame must be filtered as it contains several
## channels.
data("PGC_TCD")
PGC_TCD <- PGC_TCD \%>\% filter(channel == "WLD_RL2")

## define the column names of interest
appac_colnames <- list(
  sample_col    = "sample.name",
  peak_col      = "peak.name",
  date_col      = "injection.date",
  pressure_col  = "air.pressure",
  area_col      = "raw.area"
)

## check column, sample and peak names
PGC_TCD <- check_cols(PGC_TCD, appac_colnames)

## identify the column names of PGC_TCD
str(PGC_TCD)

## define P_ref as the median of the measured ambient pressures at injection time
P_ref <- median(PGC_TCD$Air_Pressure)

appac_control <- list(min_data_points = 5, drift_model = "linear")

## appac will throw warnings that sample and peak names have been changed
x <- suppressWarnings(
       appac(data = PGC_TCD,
       P_ref = P_ref,
       appac_control = appac_control,
       appac_colnames = appac_colnames)
     )

## View results
colors <- list(highlight_color="black", lowlight_color="darkgrey",
               line_color="darkblue", fill_color="lightblue")
plotGlobalFit(x, size = 4, colors)
}

\keyword{plots}
