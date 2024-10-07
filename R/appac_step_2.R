appac_step_2 <- function(appac, P_ref, appac_control) {

  Correction <- appac@correction
  Drift <- appac@drift
  spls <- names(Correction@samples)

  #----------------------------------------------------------
  # local fits:
  #
  # fit area vs P using compensated areas
  # fit_results: list >> fit_results_table: data.frame
  #--------------------------------------------------------

  fit_data <- lapply(seq_along(spls), function(s) {
    data.frame(
      date = Drift@samples[[s]]$date,
      P = Drift@samples[[s]]$pressure,
      Drift@samples[[s]]$compensated.raw.area
    )
  })
  fit_results <- lapply(fit_data, function(x) {
    .fit_level(
      as.matrix(x[, -c(1, 2)]),
      x[, 2] - P_ref,
      colnames(x)[-c(1, 2)],
      family = VGAM::uninormal()
    )
  })
  names(fit_results) <- spls
  for (i in seq_along(fit_results)) {
    fit_results[[i]]$n <- length(fit_data[[i]]$date)
  }
  fit_results_table <- data.frame(do.call(rbind, fit_results))
  rownames(fit_results_table) <- NULL
  Correction@local.fits <- fit_results

  #----------------------------------------------------------
  # global fit:
  #
  #----------------------------------------------------------
  fit_data <- data.frame(
    slope = fit_results_table$slope,
    area_ref = fit_results_table$area.ref,
    area_ref2 = fit_results_table$area.ref^2,
    weight = fit_results_table$n
  )

  global_fit <- stats::lm(
    slope ~ 0 + area_ref + area_ref2,
    data = fit_data,
    weights = weight
  )
  coefficients_global_fit <- coefficients(global_fit)
  names(coefficients_global_fit) <- c("kappa", "lambda")
  # p-values
  p_values <- summary(global_fit)$coefficients[, 4]
  names(p_values) <- c("kappa", "lambda")
  # std errors
  std_err <- summary(global_fit)$coefficients[, 2]
  names(std_err) <- c("kappa", "lambda")
  Correction@global.fit <- list(
    coefficients = coefficients_global_fit,
    p.values = p_values,
    standard.errors = std_err,
    P_ref = P_ref
  )

  # #----------------------------------------------------------
  # # set area_ref to the expected values obtained from global_fit
  # # and refit slope
  # #----------------------------------------------------------
  slope <- stats::predict(global_fit)
  a <- coefficients_global_fit[1] / coefficients_global_fit[2] / 2
  if (a < 0) {
    area_ref <- sqrt(abs(slope / coefficients_global_fit[2] - a^2)) + a
  } else {
    area_ref <- sqrt(abs(slope / coefficients_global_fit[2] + a^2)) - a
  }

  #----------------------------------------------------------
  # recalculate corrected + expected area
  #----------------------------------------------------------
  for (i in spls) {
    P <- P(object = Correction, sample = i)
    area <- rawAreas(object = Correction, sample = i)
    aref <- Correction@local.fits[[i]]$area.ref
    expected_area <- .get_expected_area(
      P = P,
      aref = aref,
      pref = P_ref,
      coefficients = coefficients_global_fit
    )
    corrected_area <- .get_corrected_area(
      area = area,
      P = P,
      pref = P_ref,
      coefficients_global_fit
    )
    Correction@samples[[i]]$corrected.area <- corrected_area
    Correction@samples[[i]]$expected.area <- expected_area
    colnames(Correction@samples[[i]]$corrected.area) <-
      colnames(Correction@samples[[i]]$raw.area)
    colnames(Correction@samples[[i]]$expected.area) <-
      colnames(Correction@samples[[i]]$raw.area)
  }

  #----------------------------------------------------------
  # recalculate drift
  #----------------------------------------------------------
  area_ref <- lapply(Correction@local.fits, function(x) x$area.ref)
  drift_data <- .get_drift(
    samples = Correction@samples,
    area_ref = area_ref,
    type = "corrected.area",
    drift_model = appac_control$drift_model
  )
  Drift@drift.factors <- drift_data$drift

  compensated_raw_areas <- lapply(spls, function(x) {
    .get_compensated_area(
      x,
      P_ref = P_ref,
      drift = Drift,
      correction = Correction,
      type = "raw.area",
      debias = FALSE
    )
  })
  compensated_corrected_areas <- lapply(spls, function(x) {
    .get_compensated_area(
      x,
      P_ref = P_ref,
      drift = Drift,
      correction = Correction,
      type = "corrected.area",
      debias = FALSE
    )
  })
  for (i in seq_along(compensated_raw_areas)) {
    Drift@samples[[i]] <- list(
      date = compensated_raw_areas[[i]]$date,
      pressure = compensated_raw_areas[[i]]$P,
      compensated.raw.area = as.matrix(compensated_raw_areas[[i]][, -c(1, 2)])
    )
    Correction@samples[[i]]$compensated.corrected.area <-
      as.matrix(compensated_corrected_areas[[i]][, -c(1, 2)])
  }
  names(Drift@samples) <- names(compensated_raw_areas) <- spls

  return(methods::new("Appac", drift = Drift, correction = Correction))
}
