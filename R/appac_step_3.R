appac_step_3 <- function(appac, P_ref, appac_control) {

  n <- NULL
  Correction <- appac@correction
  Drift <- appac@drift
  spls <- names(Correction@samples)

  #----------------------------------------------------------
  # Balancing of local fits:
  # compensate fit.results for bias in residuals
  #
  # we will find some overcompensation of the pressure fit
  # which can easily be seen in a plot of residuals vs. P
  #----------------------------------------------------------
  expected_areas <- lapply(seq_along(Drift@samples), function(x) {
    .get_expected_area(
      P = Drift@samples[[x]]$pressure,
      aref = Correction@local.fits[[x]]$area.ref,
      pref = P_ref,
      coefficients = Correction@global.fit$coefficients
    )
  })
  residual_areas <- lapply(seq_along(Drift@samples), function(x) {
    Drift@samples[[x]]$compensated.raw.area - expected_areas[[x]]
  })
  P <- lapply(Drift@samples, function(x) x$pressure - P_ref)
  fit_bias <- lapply(
    seq_along(P),
    function(x) {
      stats::lm(formula = residual_areas[[x]] ~ P[[x]],
                na.action = stats::na.omit)
    }
  )
  cor_intercept <- sapply(fit_bias, function(x) stats::coefficients(x)[1, ])
  cor_slope <- sapply(fit_bias, function(x) coefficients(x)[2, ])
  std_errors <- sapply(fit_bias, function(x) sqrt(diag(stats::vcov(x))))
  idx <- rep(c(TRUE, FALSE), length(std_errors) / 2)
  if (length(Drift@samples) == 1) {
    cor_intercept <- as.numeric(cor_intercept)
    cor_slope <- as.numeric(cor_slope)
    cor_se_intercept <- std_errors[idx]
    cor_se_slope <- std_errors[!idx]
    Correction@local.fits[[1]]$area.ref <-
      Correction@local.fits[[1]]$area.ref + cor_intercept
    Correction@local.fits[[1]]$slope <-
      Correction@local.fits[[1]]$slope + cor_slope
    Correction@local.fits[[1]]$se.area.ref <-
      sqrt(Correction@local.fits[[1]]$se.area.ref^2 + cor_se_intercept^2)
    Correction@local.fits[[1]]$se.slope <-
      sqrt(Correction@local.fits[[1]]$se.slope^2 + cor_se_slope^2)
  } else {
    cor_intercept <- sapply(fit_bias, function(x) coefficients(x)[1, ])
    cor_slope <- sapply(fit_bias, function(x) coefficients(x)[2, ])
    std_errors <- lapply(fit_bias, function(x) sqrt(diag(vcov(x))))
    idx <- rep(c(TRUE, FALSE), length(std_errors) / 2)
    cor_se_intercept <- sapply(std_errors, function(x) unname(x[idx]))
    cor_se_slope <- sapply(std_errors, function(x) unname(x[!idx]))
    for (i in seq_along(Correction@samples)) {
      Correction@local.fits[[i]]$area_ref <-
        Correction@local.fits[[i]]$area.ref + cor_intercept[[i]]
      Correction@local.fits[[i]]$slope <-
        Correction@local.fits[[i]]$slope + cor_slope[[i]] # "-"
      Correction@local.fits[[i]]$se.area_ref <-
        Correction@local.fits[[i]]$se.area.ref + cor_se_intercept[[i]]
      Correction@local.fits[[i]]$se.slope <-
        Correction@local.fits[[i]]$se.slope + cor_se_slope[[i]]
    }
  }

  #----------------------------------------------------------
  # global fit:
  #----------------------------------------------------------
  fit_data <- data.frame(do.call(rbind, Correction@local.fits))
  global_fit <- stats::lm(slope ~ 0 + area.ref, data = fit_data, weights = n)
  coefficients_global_fit <- coefficients(global_fit)
  names(coefficients_global_fit) <- "kappa"
  # p-values
  p_values <- summary(global_fit)$coefficients[, 4]
  names(p_values) <- "p.value.kappa"
  # std errors
  std_err <- summary(global_fit)$coefficients[, 2]
  names(std_err) <- "std.error.kappa"
  Correction@global.fit <- list(
    coefficients = coefficients_global_fit,
    p.values = p_values,
    standard.errors = std_err,
    P_ref = P_ref
  )

  #----------------------------------------------------------
  # recalculate corrected + expected area
  #----------------------------------------------------------
  for (i in spls) {
    area_ref <- Correction@local.fits[[i]]$area_ref
    P <- P(object = Correction, sample = i)
    area <- rawAreas(object = Correction, sample = i)
    aref <- Correction@local.fits[[i]]$area.ref
    expected.area <- .get_expected_area(
      P = P,
      aref = aref,
      pref = P_ref,
      coefficients = coefficients_global_fit
    )
    corrected.area <- .get_corrected_area(
      area = area,
      P = P,
      pref = P_ref,
      coefficients = coefficients_global_fit
    )
    Correction@samples[[i]]$corrected.area <- corrected.area
    Correction@samples[[i]]$expected.area <- expected.area
    colnames(Correction@samples[[i]]$corrected.area) <-
      colnames(Correction@samples[[i]]$raw.area)
    colnames(Correction@samples[[i]]$expected.area) <-
      colnames(Correction@samples[[i]]$raw.area)
  }

  #----------------------------------------------------------
  # recalculate drift
  #----------------------------------------------------------
  area_ref <- lapply(Correction@local.fits, function(x) x$area.ref)
  drift.data <- .get_drift(Correction@samples,
    area_ref = area_ref,
    type = "corrected.area",
    drift_model = appac_control$drift_model
  ) # "linear")  # "corrected.area"
  Drift@drift.factors <- drift.data$drift

  compensated_raw_areas <- lapply(spls, function(x) {
    .get_compensated_area(x,
                          P_ref = P_ref,
                          drift = Drift,
                          correction = Correction,
                          type = "raw.area")
  })
  compensated_corrected_areas <- lapply(spls, function(x) {
    .get_compensated_area(x,
                          P_ref = P_ref,
                          drift = Drift,
                          correction = Correction,
                          type = "corrected.area")
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
