appac_step_2 <- function(appac, P.ref, appac.control) {
  Correction <- appac@correction
  Drift <- appac@drift
  spls <- names(Correction@samples)

  #----------------------------------------------------------
  # local fits:
  #
  # fit area vs P using compensated areas
  # fit.results: list >> fit.results.table: data.frame
  #--------------------------------------------------------

  fit.data <- lapply(seq_along(spls), function(s) {
    data.frame(
      date = Drift@samples[[s]]$date,
      P = Drift@samples[[s]]$pressure,
      Drift@samples[[s]]$compensated.raw.area
    )
  })
  fit.results <- lapply(fit.data, function(x) {
    fit.level(as.matrix(x[, -c(1, 2)]), x[, 2] - P.ref, colnames(x)[-c(1, 2)], family = VGAM::uninormal()) # VGAM::cauchy()
  })
  names(fit.results) <- spls
  for (i in seq_along(fit.results)) {
    fit.results[[i]]$n <- length(fit.data[[i]]$date)
  }
  fit.results.table <- data.frame(do.call(rbind, fit.results))
  rownames(fit.results.table) <- NULL
  Correction@local.fits <- fit.results

  #----------------------------------------------------------
  # global fit:
  #
  #----------------------------------------------------------
  fit.data <- data.frame(
    slope = fit.results.table$slope,
    area.ref = fit.results.table$area.ref,
    area.ref2 = fit.results.table$area.ref^2,
    weight = fit.results.table$n
  )

  global.fit <- stats::lm(slope ~ 0 + area.ref + area.ref2, data = fit.data, weights = weight)
  coefficients.global.fit <- coefficients(global.fit)
  names(coefficients.global.fit) <- c("kappa", "lambda")
  # p-values
  p.values.global.fit <- summary(global.fit)$coefficients[, 4]
  names(p.values.global.fit) <- c("kappa", "lambda")
  # std errors
  std.err.global.fit <- summary(global.fit)$coefficients[, 2]
  names(std.err.global.fit) <- c("kappa", "lambda")
  Correction@global.fit <- list(
    coefficients = coefficients.global.fit,
    p.values = p.values.global.fit,
    standard.errors = std.err.global.fit,
    P.ref = P.ref
  )
  # rm(list = c("fit.data", "fit.results", "fit.results.table", "global.fit"))

  # #----------------------------------------------------------
  # # set area.ref to the expected values obtained from global.fit
  # # and refit slope
  # #----------------------------------------------------------
  slope <- stats::predict(global.fit)
  a <- coefficients.global.fit[1] / coefficients.global.fit[2] / 2
  if (a < 0) {
    area.ref <- sqrt(abs(slope / coefficients.global.fit[2] - a^2)) + a
  } else {
    area.ref <- sqrt(abs(slope / coefficients.global.fit[2] + a^2)) - a
  }

  #----------------------------------------------------------
  # recalculate corrected + expected area
  #----------------------------------------------------------
  for (i in spls) {
    peak <- colnames(Correction@samples[[i]]$raw.area)
    area.ref <- Correction@local.fits[[i]]$area.ref
    P <- P(object = Correction, sample = i)
    date <- dates(object = Correction, sample = i)
    area <- rawAreas(object = Correction, sample = i)
    aref <- Correction@local.fits[[i]]$area.ref
    expected.area <- get.expected.area(P = P, aref = aref, P_ref = P.ref, coefficients = coefficients.global.fit)
    corrected.area <- get.corrected.area(area = area, P = P, P_ref = P.ref, coefficients.global.fit)
    Correction@samples[[i]]$corrected.area <- corrected.area
    Correction@samples[[i]]$expected.area <- expected.area
    colnames(Correction@samples[[i]]$corrected.area) <- colnames(Correction@samples[[i]]$raw.area)
    colnames(Correction@samples[[i]]$expected.area) <- colnames(Correction@samples[[i]]$raw.area)
  }

  #----------------------------------------------------------
  # recalculate drift
  #----------------------------------------------------------
  area.ref <- lapply(Correction@local.fits, function(x) x$area.ref)
  drift.data <- get.drift(
    samples = Correction@samples,
    area.ref = area.ref,
    type = "corrected.area",
    drift.model = appac.control$drift.model
  ) # "linear") # , area.ref = area.ref "corrected.area"
  Drift@drift.factors <- drift.data$drift

  compensated.raw.areas <- lapply(spls, function(x) {
    get.compensated.area(x, P_ref = P.ref, drift = Drift, correction = Correction, type = "raw.area", debias = FALSE)
  })
  compensated.corrected.areas <- lapply(spls, function(x) {
    get.compensated.area(x, P_ref = P.ref, drift = Drift, correction = Correction, type = "corrected.area", debias = FALSE)
  })
  for (i in seq_along(compensated.raw.areas)) {
    Drift@samples[[i]] <- list(
      date = compensated.raw.areas[[i]]$date,
      pressure = compensated.raw.areas[[i]]$P,
      compensated.raw.area = as.matrix(compensated.raw.areas[[i]][, -c(1, 2)])
    )
    Correction@samples[[i]]$compensated.corrected.area <- as.matrix(compensated.corrected.areas[[i]][, -c(1, 2)])
  }
  names(Drift@samples) <- names(compensated.raw.areas) <- spls

  return(methods::new("Appac", drift = Drift, correction = Correction))
}
