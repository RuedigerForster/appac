appac_step_3 <- function(appac, P.ref, appac.control) {
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
  expected.areas <- lapply(seq_along(Drift@samples), function(x) {
    get.expected.area(
      Drift@samples[[x]]$pressure,
      Correction@local.fits[[x]]$area.ref,
      P.ref,
      Correction@global.fit$coefficients
    )
  })
  residual.areas <- lapply(seq_along(Drift@samples), function(x) {
    Drift@samples[[x]]$compensated.raw.area - expected.areas[[x]]
  })
  .P <- lapply(Drift@samples, function(x) x$pressure - P.ref)
  fit.bias <- lapply(seq_along(.P), function(x) stats::lm(formula = residual.areas[[x]] ~ .P[[x]], na.action = stats::na.omit)) # aref[[x]] +
  cor.intercept <- sapply(fit.bias, function(x) stats::coefficients(x)[1, ])
  cor.slope <- sapply(fit.bias, function(x) coefficients(x)[2, ])
  std.errors <- sapply(fit.bias, function(x) sqrt(diag(stats::vcov(x))))
  idx <- rep(c(T, F), length(std.errors) / 2)
  if (length(Drift@samples) == 1) {
    cor.intercept <- as.numeric(cor.intercept)
    cor.slope <- as.numeric(cor.slope)
    cor.se_intercept <- std.errors[idx]
    cor.se_slope <- std.errors[!idx]
    Correction@local.fits[[1]]$area.ref <- Correction@local.fits[[1]]$area.ref + cor.intercept # "-"
    Correction@local.fits[[1]]$slope <- Correction@local.fits[[1]]$slope + cor.slope # "-"
    Correction@local.fits[[1]]$se.area.ref <- Correction@local.fits[[1]]$se.area.ref + cor.se_intercept
    Correction@local.fits[[1]]$se.slope <- Correction@local.fits[[1]]$se.slope + cor.se_slope
  } else {
    cor.intercept <- sapply(fit.bias, function(x) coefficients(x)[1, ])
    cor.slope <- sapply(fit.bias, function(x) coefficients(x)[2, ])
    std.errors <- lapply(fit.bias, function(x) sqrt(diag(vcov(x))))
    idx <- rep(c(T, F), length(std.errors) / 2)
    cor.se_intercept <- sapply(std.errors, function(x) unname(x[idx]))
    cor.se_slope <- sapply(std.errors, function(x) unname(x[!idx]))
    for (i in seq_along(Correction@samples)) {
      Correction@local.fits[[i]]$area.ref <- Correction@local.fits[[i]]$area.ref + cor.intercept[[i]] # "-"
      Correction@local.fits[[i]]$slope <- Correction@local.fits[[i]]$slope + cor.slope[[i]] # "-"
      Correction@local.fits[[i]]$se.area.ref <- Correction@local.fits[[i]]$se.area.ref + cor.se_intercept[[i]]
      Correction@local.fits[[i]]$se.slope <- Correction@local.fits[[i]]$se.slope + cor.se_slope[[i]]
    }
  }

  #----------------------------------------------------------
  # global fit:
  #----------------------------------------------------------
  fit.data <- data.frame(do.call(rbind, Correction@local.fits))
  global.fit <- stats::lm(slope ~ 0 + area.ref, data = fit.data, weights = n) # + area.ref2
  coefficients.global.fit <- coefficients(global.fit)
  names(coefficients.global.fit) <- "kappa" # c("kappa", "lambda")
  # p-values
  p.values.global.fit <- summary(global.fit)$coefficients[, 4]
  names(p.values.global.fit) <- "p.value.kappa"
  # names(p.values.global.fit) <- "kappa"
  # std errors
  std.err.global.fit <- summary(global.fit)$coefficients[, 2]
  names(std.err.global.fit) <- "std.error.kappa"
  Correction@global.fit <- list(
    coefficients = coefficients.global.fit,
    p.values = p.values.global.fit,
    standard.errors = std.err.global.fit,
    P.ref = P.ref
  )

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
    expected.area <- get.expected.area(P = P, aref = aref, P_ref = P.ref, coefficients.global.fit)
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
  drift.data <- get.drift(Correction@samples,
    area.ref = area.ref,
    type = "corrected.area",
    drift.model = appac.control$drift.model
  ) # "linear")  # "corrected.area"
  Drift@drift.factors <- drift.data$drift

  compensated.raw.areas <- lapply(spls, function(x) {
    get.compensated.area(x, P_ref = P.ref, drift = Drift, correction = Correction, type = "raw.area")
  })
  compensated.corrected.areas <- lapply(spls, function(x) {
    get.compensated.area(x, P_ref = P.ref, drift = Drift, correction = Correction, type = "corrected.area")
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
