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
        P_ref = 1013.25,
        appac_colnames = list(sample_col = NA,
                              peak_col = NA,
                              date_col = NA,
                              pressure_col = NA,
                              area_col = NA),
        appac_control  = list(min_data_points = 50,
                              drift_model = c("linear", "quadratic"))
      )
}
\arguments{
  \item{data}{
    a data frame containing daily repetitive measurements of one or more samples
    from a single GC channel.
  }
  \item{P_ref}{
    numeric. The atmospheric pressure (in hPa) to which the instrument response
    is corrected. Any arbitrary pressure within the üressure range of the data 
    set may be chosen; however, in this case the means of the input and the 
    output areas will differ. It is recommended to choose the median of the 
    pressure measurements of the dataset.
  }
  \item{appac_colnames}{
  a list of column names of \code{data}, which contain the following information:
	    \itemize{
	      \item \code{sample_col}: character. The name of the column containing 
	      the names of the samples
	      \item \code{peak_col}: character. The name of the column containing the 
	      names of the peaks
	      \item \code{pressure_col}: character. The name of the column containing 
	      the measured atmospheric pressures
	      \item \code{area_col}: character. The name of the column containing the
	      measured peak areas
      }
  }
  \item{appac_control}{
    a list of the control elements which determine how the data are processed.
	    \itemize{
	      \item \code{min_data_points}: integer(>0). The minimum number of data
	      points belonging to a peak from which an episode is analysed.
	    	\item \code{drift_model}: character. Experimental feature, determines how drift is modeled.
	   }
  }
%%  \item{\dots}{additional parameters. Currently not in use.}
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
\seealso{
\code{\link{check_cols}}
}
\section{Warning}{....}

\examples{
## Not run:

library(dplyr)

## PGC_TCD is a short dataset of weekly measurements of 1 control sample taken on
## a dedicated process gas chromatograph. It contains more than one channel, so 
## the data.frame must be filtered.
data("PGC_TCD")
PGC_TCD <- PGC_TCD \%>\% filter(channel == "WLD_RL2")

## define the column names of the dataset
appac_colnames <- list(
  sample_col    = "sample.name",
  peak_col      = "peak.name",
  date_col      = "injection.date",
  pressure_col  = "air.pressure",
  area_col      = "raw.area"
)

## check column, sample and peak names
PGC_TCD <- check_cols(PGC_TCD, appac_colnames)

## define P_ref as the median of the measured ambient pressures at injection time
## data set contains some NA, because of missing peaks
P_ref <- median(PGC_TCD$Air_Pressure, na.rm = TRUE)


appac_control <- list(min_data_points = 5, drift_model = "linear")

## appac will throw warnings that sample and peak names have been changed
x <- appac(data = PGC_TCD, P_ref = P_ref, 
           appac_colnames = appac_colnames,
           appac_control = appac_control)

## view results
print(paste("Air pressure correction coefficient (kappa):", 
      sprintf("\%.3e", coefficients(x@correction)), "[1/hPa]"))

## plot the fit
colors <- list(highlight_color="black", lowlight_color="darkgrey", 
               line_color="darkblue", fill_color="lightblue")
plotGlobalFit(object = x, size = 2, colors = colors)

}
