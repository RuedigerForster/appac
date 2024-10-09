## ERROR: check the warning and intervention limits in residuals

.plot_global_fit <-  function(
    data = NA,
    fit_coefs = NA,
    colors = list(
      highlight_color = NA,
      lowlight_color  = NA,
      line_color      = NA,
      fill_color      = NA
    ),
    size = 4) {

  if (!all(c("x", "y", "u_x", "u_y", "n") %in% colnames(data)))
    stop("Must provide a data frame with columns: 'x', 'y', 'u_x', 'u_y', 'n'")
  if (any(is.null(colors)))
    stop("Colors missing in arguments")
  title <- "Atmospheric Pressure Peak Area Correction (APPAC) function"
  subtitle <- "Global Fit"
  xlab <- "peak area at reference pressure"
  ylab <- "slope of peak area vs. pressure"
  min <- function(x, ux) ifelse(x >= 0, x - abs(ux), x + abs(ux))
  max <- function(x, ux) ifelse(x >= 0, x + abs(ux), x - abs(ux))
  # plot uncertainties with k=2
  xmin <- min(data$x, 2 * data$u_x)
  xmax <- max(data$x, 2 * data$u_x)
  ymin <- min(data$y, 2 * data$u_y)
  ymax <- max(data$y, 2 * data$u_y)
  X <- data.frame(x = data$x, y = data$y, xmin, xmax, ymin, ymax, n = data$n)
  if (is.na(fit_coefs[2])) fit_coefs[2] <- 0
  draw_function <- function(x) return(fit_coefs[1] * x + fit_coefs[2] * x * x)
  plot <- ggplot(data = X, aes(x = x, y = y)) +
    stat_smooth(
      method = "lm",
      data = X,
      formula = y ~ poly(x, 1),
      level = 0.95,
      color = "NA",
      fill = colors$fill_color,
      alpha = 0.5
    ) +
    geom_function(
      fun = draw_function,
      linewidth = 1,
      color = colors$line_color
    ) +
    geom_point(
      size = 3,
      shape = 20,
      color = colors$highlight_color
    )  +
    geom_errorbar(
      X,
      mapping = aes(x = x, ymin = ymin, ymax = ymax),
      linewidth = 0.7,
      color = colors$lowlight_color
    ) +
    geom_errorbarh(
      X,
      mapping = aes(y = y, xmin = xmin, xmax = xmax),
      linewidth = 0.7,
      color = colors$lowlight_color
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    ) +
    main_plot_theme(
      size = size
    )
  return(plot)
}


.plot_local_fit <-  function(
    data = NA,
    sample = NA,
    peak = NA,
    colors = list(
      highlight_color = NA,
      lowlight_color  = NA,
      line_color      = NA,
      fill_color      = NA
    ),
    size = 4,
    show.compensated.areas = TRUE,
    plot.residuals = TRUE,
    bins = 50) {
  # check colors; this error should never occur
  if (any(is.na(colors)))
    stop("Colors missing in arguments")
  title <- "Local Fit"
  subtitle <- paste(sample, "-", peak)
  xlab <- "atmospheric pressure (hPa)"
  ylab <- "peak area"
  # plot data
  x <- P(object = data@correction, sample = sample)
  y1 <- rawAreas(object = data@correction, sample = sample)[, peak]
  y2 <- compensatedRawAreas(object = data@drift, sample = sample)[, peak]
  # outliers
  outliers <- .has_outliers(data@correction, sample, conf.int = 0.95)
  if (!is.null(outliers)) { #  && outliers$pass[peak]
    s <- ifelse(outliers$outliers[, peak], 2, 1)
  } else {
    s <- rep(1, length(dates))
  }
  X <- data.frame(x = x, y1 = y1, y2 = y2, s = s)
  std_dev <- sd(X[, 3])
  aref <- X(data@correction)[[sample]]
  slop <- Y(data@correction)[[sample]]
  idx <- names(aref) == peak
  aref <- aref[idx]
  slop <- slop[idx]
  pref <- data@correction@global.fit$P_ref
  draw_function <- function(x) return(aref + slop * (x - pref))
  # plot margins
  span <- 0.025 * aref
  l.y.limit <- aref - span
  h.y.limit <- aref + span
  p1 <- ggplot(na.omit(X)) +
    # raw areas
    geom_point(
      mapping = aes(x = x, y = y1),
      color = colors$lowlight_color,
      size = 1,
      alpha = 1,
      shape = 20
    )
  # drift + bias compensated areas
  if (show.compensated.areas) {
    p1 <- p1 +
      geom_point(
        data = X,
        mapping = aes(x = x, y = y2, alpha = s),
        color = colors$highlight_color,
        size = 1,
        shape = 20
      ) +
      scale_alpha(
        range = c(1, 0.01)
      )
  }
  p1 <- p1 +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab, y = ylab
    ) +
    ylim(
      l.y.limit,
      h.y.limit
    ) +
    geom_function(
      fun = draw_function,
      linewidth = 1,
      color = colors$line_color
    ) +
    main_plot_theme(
      size = size
    )
  # residuals
  ylab <- "residuals"
  # plot data
  x <- P(object = data@correction, sample = sample)
  y1 <- residualAreas(
    object = data@correction,
    sample = sample,
    type = "raw"
  )[, peak]
  # y2 <- residualAreas(
  #   object = data@correction,
  #   sample = sample,
  #   type = "corrected"
  # )[, peak]
  
  y2 <- compensatedRawAreas(
    object = data@drift,
    sample = sample
  )[, peak] - expectedAreas(
    object = data@correction,
    sample = sample
  )[, peak]
  l.y.limit <- -span / 2.5
  h.y.limit <- span / 2.5
  Z <- data.frame(cbind(x = x, y1 = y1, y2 = y2, s = s))
  std_dev <- sd(Z[, 3])
  binwidth <- (max(Z[, 2], na.rm = TRUE) - min(Z[, 2], na.rm = TRUE)) / bins
  p2 <- ggplot(Z) +
    # raw areas
    geom_point(
      data = Z,
      mapping = aes(x = x, y = y1),
      color = colors$lowlight_color,
      size = 1,
      alpha = 0.5,
      shape = 20
    )
  # drift + bias compensated areas
  if (show.compensated.areas) {
    p2 <- p2 +
      geom_point(
        data = Z,
        mapping = aes(x = x, y = y2, alpha = s),
        color = colors$highlight_color,
        size = 1,
        shape = 20
      ) +
      scale_alpha(
        range = c(1, 0.01)
      )
  }
  p2 <- p2 +
    geom_hline(
      aes(yintercept = 0),
      color = colors$line_color,
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_hline(
      yintercept = 2 * std_dev,
      linetype = "dotted",
      color = colors$line_color,
      linewidth = 1
    ) +
    geom_hline(
      yintercept = -2 * std_dev,
      linetype = "dotted",
      color = colors$line_color,
      linewidth = 1
    ) +
    geom_hline(
      yintercept = 3 * std_dev,
      linetype = "dashed",
      color = colors$line_color,
      linewidth = 1
    ) +
    geom_hline(
      yintercept = -3 * std_dev,
      linetype = "dashed",
      color = colors$line_color,
      linewidth = 1
    ) +
    labs(
      title = element_blank(),
      subtitle = element_blank(),
      x = xlab,
      y = ylab
    ) +
    ylim(
      l.y.limit,
      h.y.limit
    ) +
    residuals_plot_theme(
      size = size
    )
  # histogram of residuals
  p3 <- ggplot(
    data = na.omit(Z),
    mapping = aes(x = y1)
  ) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = binwidth,
      na.rm = TRUE,
      position = "identity",
      fill = colors$lowlight_color,
      color = colors$lowlight_color
    )
  if (show.compensated.areas) {
    p3 <- p3 +
      geom_histogram(
        aes(x = y2, y = after_stat(density)),
        binwidth = binwidth,
        na.rm = TRUE,
        position = "identity",
        fill = colors$highlight_color,
        alpha = 0.4
      )
  }
  p3 <- p3 +
    geom_vline(
      aes(xintercept = 0),
      color = colors$line_color,
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_function(
      color = colors$line_color,
      linewidth = 1,
      fun = dnorm,
      args = list(mean = mean(y1, na.rm = TRUE), sd = sd(y1, na.rm = TRUE))
    ) +
    histogram_plot_theme(
      size = size
    ) +
    xlim(
      l.y.limit,
      h.y.limit
    ) +
    coord_flip()

  p  <- ggpubr::ggarrange(
    p1,
    p2 + p3 + patchwork::plot_layout(widths = c(3, 1)),
    ncol = 1, nrow = 2, heights = c(6, 3)
  )
  return(p)
}



.plot_control_chart <- function(
    data = NA,
    sample = NA,
    peak = NA,
    colors = list(
      highlight_color = NA,
      lowlight_color  = NA,
      line_color      = NA,
      fill_color      = NA
    ),
    size = 4,
    show.compensated.areas = TRUE,
    show.breakpoints = TRUE,
    plot.residuals = TRUE,
    bins = 50) {

  if (any(is.na(colors)))
    stop("Colors missing in arguments")

  title <- "Control Chart"
  subtitle <- paste(sample, "-", peak)
  xlab <- "date"
  ylab <- "peak area"
  # plot data
  x <- dates(object = data@correction, sample = sample)
  y1 <- rawAreas(object = data@correction, sample = sample)[, peak]
  y2 <- compensatedCorrectedAreas(
    object = data@correction,
    sample = sample
  )[, peak]
  # outliers
  outliers <- .has_outliers(data@correction, sample, conf.int = 0.95)
  if (!is.null(outliers)) { #  && outliers$pass[peak]
    s <- ifelse(outliers$outliers[, peak], 2, 1)
  } else {
    s <- rep(1, length(dates))
  }
  # plot data
  X <- data.frame(x = x, y1 = y1, y2 = y2, s = s)
  X[, 1] <- as.Date(X[, 1], origin = "1970-01-01")
  std_dev <- sd(X[, 3])
  # plot margins
  aref <- X(data@correction)[[sample]]
  idx <- names(aref) == peak
  aref <- aref[idx]
  span <- 0.025 * aref
  l.y.limit <- aref - span
  h.y.limit <- aref + span
  brk <- c(1, 2, 3, 4, 5, 6, 9, 12)
  nb <- as.numeric((max(x) - min(x)) / 300)
  y.breaks <- brk[which(abs(brk - nb) == min(abs(brk - nb)))]
  # control chart
  p1 <- ggplot(data = na.omit(X)) +
    geom_point(
      mapping = aes(x = x, y = y1),
      color = colors$lowlight_color, size = 2,
      alpha = 1,
      shape = 20
    ) +
    scale_x_date(
      date_breaks = paste0(y.breaks, " months"),
      date_labels = "%m/%y"
    )
  if (show.compensated.areas) {
    p1 <- p1 +
      geom_point(
        data = X,
        mapping = aes(x = x, y = y2, alpha = s),
        color = colors$highlight_color,
        size = 1,
        shape = 20
      ) +
      scale_alpha(
        range = c(1, 0.01)
      )
  }
  p1 <- p1 +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    ) +
    ylim(
      l.y.limit,
      h.y.limit
    ) +
    geom_hline(
      yintercept = aref,
      linetype = "solid",
      color = colors$line_color,
      linewidth = 1
    ) +
    geom_hline(
      yintercept = aref + 2 * std_dev,
      linetype = "dotted",
      color = colors$line_color,
      linewidth = 1
    ) +
    geom_hline(
      yintercept = aref - 2 * std_dev,
      linetype = "dotted",
      color = colors$line_color,
      linewidth = 1
    ) +
    geom_hline(
      yintercept = aref + 3 * std_dev,
      linetype = "dashed",
      color = colors$line_color,
      linewidth = 1
    ) +
    geom_hline(
      yintercept = aref - 3 * std_dev,
      linetype = "dashed",
      color = colors$line_color,
      linewidth = 1
    ) +
    main_plot_theme(
      size = size
    )
  # breakpoints
  if (show.breakpoints) {
    bpts <- data@correction@samples[[sample]]$breakpoints
    p1 <- p1 +
      sapply(bpts, function(a) geom_vline(xintercept=a, linetype = "dotted", color = colors$line_color, linewidth = 0.5))
  }
  #
  # residuals
  ylab <- "residuals"
  # plot data
  x <- dates(object = data@correction, sample = sample)
  y1 <- residualAreas(
    object = data@correction,
    sample = sample,
    type = "raw.area"
  )[, peak]
  y2 <- residualAreas(
    object = data@correction,
    sample = sample,
    type = "compensated.corrected.area"
  )[, peak]
    
  #   compensatedCorrectedAreas(
  #   object = data@correction,
  #   sample = sample
  # )[, peak] - aref
  l.y.limit <- -span / 2.5
  h.y.limit <- span / 2.5
  Z <- as.data.frame(cbind(x = x, y1 = y1, y2 = y2, s = s))
  Z[, 1] <- as.Date(Z[, 1], origin = "1970-01-01")
  binwidth <- (max(Z[, 2], na.rm = TRUE) - min(Z[, 2], na.rm = TRUE)) / bins
  p2 <- ggplot(Z) +
    # raw areas
    geom_point(
      data = Z,
      mapping = aes(x = x, y = y1),
      color = colors$lowlight_color,
      size = 1,
      alpha = 0.5,
      shape = 20
    ) +
    labs(
      title = element_blank(),
      subtitle = element_blank(),
      x = xlab,
      y = ylab
    ) +
    ylim(
      l.y.limit,
      h.y.limit
    ) +
    scale_x_date(
      date_breaks = paste0(y.breaks, " months"),
      date_labels = "%m/%y"
    )
  # drift + bias compensated areas
  if (show.compensated.areas) {
    p2 <- p2 +
      geom_point(
        data = Z,
        mapping = aes(x = x, y = y2, alpha = s),
        color = colors$highlight_color,
        size = 1,
        shape = 20
      ) +
      scale_alpha(
        range = c(1, 0.01)
      )
  }
  p2 <- p2 +
    geom_hline(
      aes(yintercept = 0),
      color = colors$line_color,
      linetype = "dashed",
      linewidth = 1
    ) +
    residuals_plot_theme(
      size = size
    )
  # breakpoints
  if (show.breakpoints) {
    bpts <- data@correction@samples[[sample]]$breakpoints
    p2 <- p2 +
      sapply(bpts, function(a) geom_vline(xintercept=a, linetype = "dotted", color = colors$line_color, linewidth = 0.5))
  }
  # histogram of residuals
  p3 <- ggplot(data = na.omit(Z), mapping = aes(x = y1)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = binwidth,
      na.rm = TRUE,
      position = "identity",
      fill = colors$lowlight_color,
      color = colors$lowlight_color
    )
  if (show.compensated.areas) {
    p3 <- p3 +
      geom_histogram(
        aes(x = y2, y = after_stat(density)),
        binwidth = binwidth,
        na.rm = TRUE,
        position = "identity",
        fill = colors$highlight_color,
        alpha = 0.4
      )
  }
  p3 <- p3 +
    geom_vline(
      aes(xintercept = 0),
      color = colors$line_color,
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_function(
      color = colors$line_color,
      linewidth = 1,
      fun = dnorm,
      args = list(mean = mean(y1, na.rm = TRUE), sd = sd(y1, na.rm = TRUE))
    ) +
    histogram_plot_theme(
      size = size
    ) +
    xlim(
      l.y.limit,
      h.y.limit
    ) +
    coord_flip()

  p  <- ggpubr::ggarrange(
    p1, p2 + p3 + patchwork::plot_layout(widths = c(3, 1)),
    ncol = 1, nrow = 2, heights = c(6, 3)
  )
  return(p)
}
