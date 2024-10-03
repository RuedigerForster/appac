#
# library(stats) # rstudent, lm, cooks.distance, cor.test
# library(dplyr) # pipe operator, filter, mutate
# library(data.table) # definition of IDate, as.IDate
# library(Rbeast) # beast
#
# source("./R/helper.R")
# source("./R/S4_definitions.R")
# source("./R/appac_step_1.R")
# source("./R/appac_step_2.R")
# source("./R/appac_step_3.R")

#
# appac.control$fit.drift = c("linear","quadratic")
#

appac <- function(data,
                  P.ref = 1013.25,
                  appac.colnames = list(
                    sample.col = NA,
                    peak.col = NA,
                    date.col = NA,
                    pressure.col = NA,
                    area.col = NA
                  ),
                  appac.control = list(
                    min.data.points = 50,
                    drift.model = c("linear", "quadratic")
                  )) {
  if (missing(data)) {
    stop("Input 'data' is missing")
  }

  data <- check.cols(data, appac.colnames)

  appac.control$drift.model <- match.arg(appac.control$drift.model, c("linear", "quadratic"))

  #----------------------------------------------------------
  # check if min.data.points is a positive integer
  #----------------------------------------------------------
  appac.control$min.data.points <- as.integer(appac.control$min.data.points)
  if (appac.control$min.data.points < 0L) {
    stop("'appac.control$min.data.points' must have a positive value")
  }

  #----------------------------------------------------------
  # check if P.ref is within the range of the input pressures
  #----------------------------------------------------------
  if (!is.numeric(P.ref) | P.ref < min(data$air.pressure, na.rm = TRUE) | P.ref > max(data$air.pressure, na.rm = TRUE)) {
    stop("P.ref: ", P.ref, " is out of range: ", min(data$air.pressure, na.rm = TRUE), " <  P.ref < ", max(data$air.pressure, na.rm = TRUE))
  }

  #----------------------------------------------------------
  # check if data contains all peaks for every file
  # nrow(data) must be a multiple of length(peaks)
  #----------------------------------------------------------
  if (nrow(data) %% length(unique(data$peak.name)) != 0) {
    stop("Missing data points")
  }




  step_1 <- appac_step_1(data, P.ref, appac.control)
  step_2 <- appac_step_2(step_1, P.ref, appac.control)
  step_3 <- appac_step_3(step_2, P.ref, appac.control)
  return(step_3)
}
