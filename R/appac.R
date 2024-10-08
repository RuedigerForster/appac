
appac <- function(
  data,
  P_ref = 1013.25,
  appac_control = list(
    min_data_points = 50,
    drift_model = c("linear", "quadratic")
  )
) {
  
  # appac_colnames = list(
  #   sample_col = NA,
  #   peak_col = NA,
  #   date_col = NA,
  #   pressure_col = NA,
  #   area_col = NA
  # ),
  

  if (missing(data)) {
    stop("Input 'data' is missing")
  }

  if (!all(colnames(data) %in% column_names)) {
    stop("Unknown column names. Please call 'check_cols(data, appac_colnames)' before calling 'appac()'.")
  }

  appac_control$drift_model <- match.arg(
    appac_control$drift_model, c("linear", "quadratic")
  )

  #----------------------------------------------------------
  # check if min_data_points is a positive integer
  #----------------------------------------------------------
  appac_control$min_data_points <- as.integer(appac_control$min_data_points)
  if (appac_control$min_data_points < 0L) {
    stop("'appac_control$min_data_points' must have a positive value")
  }

  #----------------------------------------------------------
  # check if P_ref is within the range of the input pressures
  #----------------------------------------------------------
  if (!is.numeric(P_ref) || P_ref < min(data$Air_Pressure, na.rm = TRUE) ||
        P_ref > max(data$Air_Pressure, na.rm = TRUE)) {
    stop("P_ref: ", P_ref, " is out of range: ",
      min(data$Air_Pressure, na.rm = TRUE), " <  P_ref < ",
      max(data$Air_Pressure, na.rm = TRUE)
    )
  }

  #----------------------------------------------------------
  # check if data contains all peaks for every file
  # nrow(data) must be a multiple of length(peaks)
  #----------------------------------------------------------
  if (nrow(data) %% length(unique(data$Peak_Name)) != 0) {
    stop("Missing data points")
  }




  step_1 <- appac_step_1(data, P_ref, appac_control)
  step_2 <- appac_step_2(step_1, P_ref, appac_control)
  step_3 <- appac_step_3(step_2, P_ref, appac_control)
  return(step_3)
}
