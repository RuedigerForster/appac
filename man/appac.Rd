\name{appac}
\alias{appac}
\docType{methods}
\title{
Atmospheric Pressure Peak Area Correction (APPAC)
}
\description{
APPAC is a method for posterior correction of peak areas obtained from GC/FID or
GC/TCD chromatograms. The corrections are obtained from daily repetitive
measurements on dedicated control samples. The algorithm takes into account the
influences of atmospheric pressure and instrument drift. Change points, such as
abrupt changes, in the data are detected and any resulting bias is compensated.
For this purpose, the data are divided into episodes at the breakpoints. The
episodes are analysed separately and then recombined. Additionally, the detector
response is linearized in the output.
}
\usage{
appac ( data,
        P.ref = 1013.25,
        appac.control  = list(min.data.points = 50),
        appac.colnames = list(sample.col = NA,
                              peak.col = NA,
                              date.col = NA,
                              pressure.col = NA,
                              area.col = NA)
      )
}
\arguments{
  \item{data}{
    a data frame containing daily repetitive measurements of one or more samples
    from a single GC channel.
  }
  \item{P.ref}{
    numeric. The atmospheric pressure (in hPa) to which the instrument response
    is corrected. Any arbitrary pressure within the üressure range of the data 
    set may be chosen; however, in this case the means of the input and the 
    output areas will differ. It is recommended to choose the median of the 
    pressure measurements of the dataset.
  }
  \item{appac.control}{
    a list of the control elements which determine how the data are processed.
	    \itemize{
	      \item \code{min.data.points:} integer(>0). The minimum number of data
	      points belonging to a peak from which an episode is analysed.
	    }
    }
  \item{appac.colnames}{
  a list of column names of \code{df}, which contain the following information:
	    \itemize{
	      \item \code{sample.col:} character. The name of the column containing 
	      the names of the samples
	      \item \code{peak.col:} character. The name of the column containing the 
	      names of the peaks
	      \item \code{pressure.col:} character. The name of the column containing 
	      the measured atmospheric pressures
	      \item \code{area.col:} character. The name of the column containing the
	      measured peak areas
      }
  }
}
\details{
Sample and peak names in the input data frame are used as variables in the code
and thus should comply to \bold{R} naming conventions. A valid variable name must
start with a letter and may contain letters, digits, period(.) and
underscore(_). If this is not the case, non-complying names will be converted
and a warning will be thrown.

Dates in the input dataframe must comply to the format \code{\%Y-\%m-\%d}.
See \code{\link{strptime}} for further information.

appac makes extensive use of the correlation, which can be found in the peak
areas of the chromatograms of identical samples. It is reuired that the
chromatograms to be evaluated contain at least two evaluated peaks, however
more than two evaluated peaks are highly recommended.
}
\value{
an object of class \code{\link{Appac}}
}

\section{Warning}{....}

\examples{
## Please be patient, the calculation takes a while.

## Not run:
library(appac)

## PLOT_FID is a long dataset of daily measurements of 5 control samples taken on
## a dedicated natural gas analyzer.
data("PLOT_FID")

## truncate the data set in favor of a shorter calculation time
PLOT_FID <- PLOT_FID[1:30000, ]


## define P.ref as the median of the measured ambient pressures at injection time
## data set contains some NA, because of missing peaks
P.ref <- median(PLOT_FID$air.pressure, na.rm = TRUE)

## define the column names of the dataset
appac.colnames <- list(
  sample.col    = "sample.name",
  peak.col      = "peak.name",
  date.col      = "injection.date",
  pressure.col  = "air.pressure",
  area.col      = "raw.area"
)

appac.control <- list(min.data.points = 50)

## appac will throw warnings that sample and peak names have been changed
x <- appac(data = PLOT_FID, P.ref = P.ref, 
           appac.control = appac.control, 
           appac.colnames = appac.colnames)

## view results
print(paste("Air pressure correction coefficient (kappa):", 
      sprintf("\%.3e", coefficients(x@correction)), "[1/hPa]"))

## plot the fit
colors <- list(highlight_color="black", lowlight_color="darkgrey", 
               line_color="darkblue", fill_color="lightblue")
plotGlobalFit(x, colors)

}
