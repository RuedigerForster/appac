\docType{methods}
\name{plotControlChart}
\alias{plotControlChart}
\alias{plotControlChart,Appac,character,character-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the raw and fitted area vs. time of a peak of a sample
}
\description{
The control chart is a plot of the raw and the fitted (pressure corrected and
drift compensated) areas vs. time. It is a variant of a Shewart control chart.
}
\usage{
\S4method{plotControlChart}{Appac,character,character}(object, sample, peak, size = 4, ...)
}
\arguments{
  \item{object}{
an object of class \code{\link{Appac}}
}
  \item{sample}{
a character value representing the name of the sample
}
  \item{peak}{
a character value representing the name of the peak
}
  \item{size}{
a numeric value which controls the size of the axis labels and titles
}
  \item{\dots}{
    additional arguments:

    \code{colors}: a named list containing the color scheme to be used in the plot:
    \describe{
%      \itemize{
        \item{highlight_color:}{the color in which the datapoints of the
        corrected areas are represented}
        \item{lowlight_color:}{the color in which the datapoints of the raw
        areas are represented}
        \item{line_color:}{the color of the fitted function line}
        \item{fill_color:}{the color in which the standard uncertainties of the
        fitted function line are plotted}
%      }
    }
  }
}
\value{
a ggplot object
}

\seealso{
\code{\link{plotLocalFit}}
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
x <- appac(data = PGC_TCD,
           P_ref = P_ref,
           appac_control = appac_control,
           appac_colnames = appac_colnames)

## View results
s <- "X17k"
p <- "iso.Butan"
colors <- list(highlight_color="black", lowlight_color="darkgrey",
               line_color="darkblue", fill_color="lightblue")

plotControlChart ( object = x, sample = s, peak = p, size = 4, colors = colors)
}

\keyword{plots}
