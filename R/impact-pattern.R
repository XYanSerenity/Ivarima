#' Add an unevaluated candidate model
#'
#' This internal helper appends one candidate transfer-function specification to
#' the unevaluated candidate list if the same key has not already been used.
#'
#' @param key A character string used as the candidate identifier.
#' @param xdata A numeric vector or one-column data frame representing the intervention variable.
#' @param xdatar Integer. Numerator order of the transfer function.
#' @param xdatab Integer. Delay parameter of the transfer function.
#' @param xdatah Integer. Denominator order of the transfer function.
#' @param unevaluated_list A list of unevaluated candidate models.
#'
#' @return An updated list of unevaluated candidate models.
#' @keywords internal
add_unevaluated <- function(key, xdata, xdatar, xdatab, xdatah, unevaluated_list) {
  if (!key %in% names(unevaluated_list)) {
    unevaluated_list[[key]] <- list(
      xdata = xdata,
      xdatar = xdatar,
      xdatab = xdatab,
      xdatah = xdatah
    )
  }

  unevaluated_list
}

#' Safely add an unevaluated candidate model
#'
#' @inheritParams add_unevaluated
#' @keywords internal
add_unevaluated_safe <- function(key, xdata, xdatar, xdatab, xdatah, unevaluated_list) {
  if (!key %in% names(unevaluated_list)) {
    unevaluated_list <- add_unevaluated(
      key = key,
      xdata = xdata,
      xdatar = xdatar,
      xdatab = xdatab,
      xdatah = xdatah,
      unevaluated_list = unevaluated_list
    )
  }

  unevaluated_list
}

#' Lag a vector and fill leading positions with zero
#'
#' @param x A numeric vector.
#' @param lag_n A non-negative integer specifying the lag length.
#'
#' @return A numeric vector with the same length as `x`.
#' @export
#'
#' @examples
#' lag0_fill(1:5, 2)
lag0_fill <- function(x, lag_n = 0L) {
  lag_n <- as.integer(lag_n)
  x <- as.numeric(x)
  n <- length(x)

  if (lag_n <= 0L) {
    return(x)
  }

  if (lag_n >= n) {
    return(rep(0, n))
  }

  c(rep(0, lag_n), x[seq_len(n - lag_n)])
}

#' Find the middle position between the strongest positive and negative coefficients
#'
#' @param x A numeric vector.
#'
#' @return An integer position, or 0 if no valid middle position is found.
#' @keywords internal
find_middle_position <- function(x) {
  pos_idx <- which(x > 0)
  neg_idx <- which(x < 0)

  if (length(pos_idx) == 0L || length(neg_idx) == 0L) {
    return(0L)
  }

  max_pos_idx <- pos_idx[which.max(x[pos_idx])]
  min_neg_idx <- neg_idx[which.min(x[neg_idx])]

  if (abs(max_pos_idx - min_neg_idx) == 2L) {
    return(seq(min(max_pos_idx, min_neg_idx) + 1L, max(max_pos_idx, min_neg_idx) - 1L))
  }

  0L
}

#' Extract ARIMA orders from a fitted forecast::Arima object
#'
#' @param model_ar A fitted ARIMA model returned by `forecast::auto.arima()`.
#'
#' @return A numeric vector in the form c(p, d, q).
#' @keywords internal
extract_arima_order <- function(model_ar) {
  model_ar$arma[c(1, 6, 2)]
}

#' Find the first sign-change position
#'
#' @param x A numeric vector.
#'
#' @return The first position before a sign change, or `NULL` if no sign change is detected.
#' @keywords internal
find_sign_change_position <- function(x) {
  s <- sign(x)
  changes <- which(diff(s) != 0)

  if (length(changes) > 0L) {
    changes[1]
  } else {
    NULL
  }
}

#' Find the strongest adjacent sign difference
#'
#' @param coefs A numeric vector of coefficients.
#'
#' @return The index before the strongest sign reversal, or 0 if none is found.
#' @keywords internal
find_max_sign_diff <- function(coefs) {
  n <- length(coefs)

  if (n < 2L) {
    return(0L)
  }

  max_diff <- 0
  max_diff_index <- 0L

  for (i in seq_len(n - 1L)) {
    if (coefs[i] * coefs[i + 1L] < 0) {
      diff_value <- abs(coefs[i] - coefs[i + 1L])

      if (diff_value > max_diff) {
        max_diff <- diff_value
        max_diff_index <- i
      }
    }
  }

  max_diff_index
}

#' Fit the linear transfer-function model for initial v-weight identification
#'
#' @param y A numeric outcome time series.
#' @param x_p A numeric intervention input series.
#' @param order A numeric vector specifying ARIMA order c(p, d, q).
#' @param k Integer. Number of lagged intervention terms used to estimate v-weights.
#'
#' @return A fitted ARIMA model, or `NULL` if model fitting fails.
#' @keywords internal
fit_ltf_model <- function(y, x_p, order, k) {
  x_p <- as.numeric(x_p)
  y <- as.numeric(y)
  k <- as.integer(k)

  if (length(x_p) < k || length(y) < k || k < 1L) {
    return(NULL)
  }

  lagged_xp <- stats::embed(x_p, k)
  colnames(lagged_xp) <- paste0("xpv", seq_len(k) - 1L)
  y_clean <- utils::tail(y, nrow(lagged_xp))

  fit <- try(
    stats::arima(y_clean, order = order, xreg = lagged_xp, method = "ML"),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    NULL
  } else {
    fit
  }
}

#' Identify candidate delay parameters from the estimated v-weight pattern
#'
#' @param ltf_model A fitted LTF model returned by `fit_ltf_model()`.
#' @param k Integer. Number of estimated v-weights.
#' @param threshold Numeric. Relative threshold used for relaxed matching.
#'
#' @return A list containing candidate delay values in element `b`.
#' @keywords internal
identify_candidate_b <- function(ltf_model, k, threshold = 0.3) {
  if (is.null(ltf_model)) {
    return(list(b = 0L))
  }

  coef_table <- lmtest::coeftest(ltf_model)
  v_coef <- utils::tail(coef_table, k)
  v_coef <- v_coef[!is.na(v_coef[, 1]), , drop = FALSE]

  if (nrow(v_coef) == 0L) {
    return(list(b = 0L))
  }

  df_v <- data.frame(v_coef)
  colnames(df_v) <- c("Coefficient", "Std_Error", "Z_Value", "P_Value")

  peak_pos <- which.max(abs(df_v$Coefficient))
  peak_value <- df_v$Coefficient[peak_pos]

  if (is.na(peak_value) || peak_value == 0) {
    return(list(b = 0L))
  }

  if (peak_pos == 1L) {
    sign_pos <- find_max_sign_diff(df_v$Coefficient)

    if (is.null(sign_pos) || length(sign_pos) == 0L || sign_pos == 0L) {
      return(list(b = 0L))
    }

    if (abs(df_v$Coefficient[sign_pos] / peak_value) < threshold &&
        abs(df_v$Coefficient[sign_pos + 1L] / peak_value) < threshold) {
      return(list(b = 0L))
    }

    b_values <- unique(c(0L, sign_pos - 1L))
    b_values <- b_values[!is.na(b_values) & b_values >= 0L]
    return(list(b = b_values))
  }

  silent_count <- 0L

  for (idx in seq_len(peak_pos - 1L)) {
    if (abs(df_v$Coefficient[idx] / peak_value) < threshold) {
      silent_count <- silent_count + 1L
    } else {
      break
    }
  }

  last_silent_pos <- NA_integer_

  for (idx in seq(from = peak_pos - 1L, to = 1L)) {
    if (abs(df_v$Coefficient[idx] / peak_value) <= threshold) {
      last_silent_pos <- idx
      break
    }
  }

  sign_change_pos <- find_sign_change_position(df_v$Coefficient)

  if (is.null(sign_change_pos)) {
    sign_change_pos <- NA_integer_
  }

  middle_pos <- find_middle_position(df_v$Coefficient)

  if (is.null(middle_pos) || length(middle_pos) == 0L) {
    middle_pos <- NA_integer_
  }

  result <- c(
    peak_pos - 1L,
    peak_pos - 2L,
    silent_count,
    last_silent_pos - 1L,
    sign_change_pos - 1L,
    sign_change_pos - 2L,
    sign_change_pos - 3L,
    middle_pos - 1L
  )

  result <- unique(result[!is.na(result) & result >= 0L])

  if (length(result) == 0L) {
    result <- 0L
  }

  list(b = result)
}

#' Extract non-ARMA transfer-function coefficients
#'
#' @param model A fitted ARIMAX model.
#'
#' @return A numeric vector of transfer-function coefficients.
#' @keywords internal
extract_transfer_coef <- function(model) {
  cf <- stats::coef(model)
  cf[!grepl("ar|ma", names(cf), ignore.case = TRUE)]
}

#' Fit an ARIMAX model with transfer-coefficient initialization
#'
#' @param y A numeric outcome time series.
#' @param order A numeric vector specifying ARIMA order c(p, d, q).
#' @param xtransf A numeric vector or data frame used by `TSA::arimax()`.
#' @param transfer A transfer-function specification used by `TSA::arimax()`.
#' @param init_transfer Initial values for transfer-function parameters.
#' @param max_refit Maximum number of refitting iterations.
#'
#' @return A fitted ARIMAX model, or `NULL` if fitting fails.
#' @keywords internal
fit_optimized_arimax <- function(y, order, xtransf, transfer, init_transfer, max_refit = 5L) {
  pq_number <- order[1] + order[3]
  init_coef <- c(rep(NA, pq_number), init_transfer)

  current_model <- try(
    TSA::arimax(
      y,
      order = order,
      xtransf = xtransf,
      transfer = transfer,
      method = "ML",
      init = init_coef
    ),
    silent = TRUE
  )

  if (inherits(current_model, "try-error")) {
    return(NULL)
  }

  for (iter in seq_len(max_refit)) {
    current_init <- c(rep(NA, pq_number), extract_transfer_coef(current_model))

    refit_model <- try(
      TSA::arimax(
        y,
        order = order,
        xtransf = xtransf,
        transfer = transfer,
        method = "ML",
        init = current_init
      ),
      silent = TRUE
    )

    if (inherits(refit_model, "try-error")) {
      break
    }

    if (abs(refit_model$loglik - current_model$loglik) < 0.1 ||
        refit_model$loglik < current_model$loglik) {
      break
    }

    current_model <- refit_model
  }

  current_model
}

#' Fit an interventional ARIMA model with optimized initialization
#'
#' @param y A numeric outcome time series.
#' @param order A numeric vector specifying ARIMA order c(p, d, q).
#' @param xtransf A numeric vector or data frame used by `TSA::arimax()`.
#' @param transfer A transfer-function specification used by `TSA::arimax()`.
#' @param max_expand Maximum number of order-expansion attempts if initial fitting fails.
#' @param max_refit Maximum number of refitting iterations.
#'
#' @return A fitted ARIMAX model, or `NULL` if fitting fails.
#' @keywords internal
optim_intervention_arimax <- function(y, order, xtransf, transfer, max_expand = 2L, max_refit = 5L) {
  target_order <- order
  init_transfer <- 0

  base_model <- try(
    TSA::arimax(y, order = order, xtransf = xtransf, transfer = transfer, method = "ML"),
    silent = TRUE
  )

  if (inherits(base_model, "try-error")) {
    current_order <- order
    found <- FALSE

    for (step in seq_len(max_expand)) {
      if (current_order[1] > 0 && current_order[1] < 3) {
        current_order <- c(current_order[1] + 1L, current_order[2], current_order[3])
      } else if (current_order[3] > 0 && current_order[3] < 3) {
        current_order <- c(current_order[1], current_order[2], current_order[3] + 1L)
      } else {
        break
      }

      expanded_model <- try(
        TSA::arimax(
          y,
          order = current_order,
          xtransf = xtransf,
          transfer = transfer,
          method = "ML"
        ),
        silent = TRUE
      )

      if (!inherits(expanded_model, "try-error")) {
        found <- TRUE
        init_transfer <- extract_transfer_coef(expanded_model)
        break
      }
    }

    if (!found) {
      return(NULL)
    }
  } else {
    init_transfer <- extract_transfer_coef(base_model)
  }

  optimized_model <- fit_optimized_arimax(
    y = y,
    order = target_order,
    xtransf = xtransf,
    transfer = transfer,
    init_transfer = init_transfer,
    max_refit = max_refit
  )

  if (!inherits(base_model, "try-error") && is.null(optimized_model)) {
    return(base_model)
  }

  if (!inherits(base_model, "try-error") && !is.null(optimized_model)) {
    return(if (optimized_model$aic < base_model$aic) optimized_model else base_model)
  }

  optimized_model
}

#' Determine whether differencing is required based on residual stationarity
#'
#' @param fit A fitted model.
#' @param alpha Significance level for the KPSS test.
#' @param trace Logical. If `TRUE`, intermediate messages are printed.
#'
#' @return Logical. `TRUE` indicates that differencing is suggested.
#' @keywords internal
need_diff_by_resid <- function(fit, alpha = 0.05, trace = FALSE) {
  if (is.null(fit)) {
    return(TRUE)
  }

  residual_values <- tryCatch(stats::residuals(fit), error = function(e) NULL)

  if (is.null(residual_values)) {
    return(TRUE)
  }

  residual_values <- residual_values[is.finite(residual_values)]
  n <- length(residual_values)

  if (n < 10L) {
    return(FALSE)
  }

  kpss_p <- tryCatch(
    tseries::kpss.test(residual_values, null = "Level")$p.value,
    error = function(e) NA_real_
  )

  if (trace) {
    message("need_diff_by_resid: KPSS p-value = ", kpss_p)
  }

  !is.na(kpss_p) && kpss_p < alpha
}

#' Calculate model information criterion
#'
#' @param fit A fitted model.
#' @param ic A character string specifying the criterion: `"aicc"`, `"aic"`, or `"bic"`.
#' @param n Optional sample size.
#'
#' @return A numeric information criterion value.
#' @keywords internal
get_model_ic <- function(fit, ic = c("aicc", "aic", "bic"), n = NULL) {
  ic <- match.arg(ic)

  if (is.null(fit)) {
    return(Inf)
  }

  if (is.null(n)) {
    n <- length(stats::residuals(fit))
  }

  k <- tryCatch(length(stats::coef(fit)), error = function(e) NA_integer_)
  aic_value <- tryCatch(stats::AIC(fit), error = function(e) Inf)
  bic_value <- tryCatch(stats::BIC(fit), error = function(e) Inf)

  if (ic == "aic") {
    return(aic_value)
  }

  if (ic == "bic") {
    return(bic_value)
  }

  if (is.na(k) || is.infinite(aic_value) || (n - k - 1L) <= 0L) {
    return(Inf)
  }

  aic_value + 2 * k * (k + 1L) / (n - k - 1L)
}

#' Search ARIMA orders under a given transfer function
#'
#' @param y A numeric outcome time series.
#' @param xtransf A numeric vector or data frame used by `TSA::arimax()`.
#' @param transfer A transfer-function specification used by `TSA::arimax()`.
#' @param max.p Maximum AR order.
#' @param max.d Maximum differencing order.
#' @param max.q Maximum MA order.
#' @param ic Information criterion used for comparison.
#' @param diff_alpha Significance level used in residual stationarity testing.
#' @param trace Logical. If `TRUE`, intermediate messages are printed.
#'
#' @return A fitted ARIMAX model, or `NULL` if no valid model is found.
#' @keywords internal
search_arimax_order <- function(
    y,
    xtransf,
    transfer,
    max.p = 5L,
    max.d = 1L,
    max.q = 5L,
    ic = "aicc",
    diff_alpha = 0.05,
    trace = FALSE) {
  y <- as.numeric(y)
  n <- length(y)

  init_fit <- try(
    optim_intervention_arimax(y, c(1L, 0L, 0L), xtransf, transfer),
    silent = TRUE
  )

  if (inherits(init_fit, "try-error") || is.null(init_fit)) {
    d_use <- 0L

    if (trace) {
      message("Initial ARIMA(1,0,0) fit failed; d = 0 is used.")
    }
  } else {
    need_diff <- need_diff_by_resid(init_fit, alpha = diff_alpha, trace = trace)
    d_use <- if (need_diff && max.d >= 1L) 1L else 0L

    if (trace) {
      message("Initial order = (1,0,0); need_diff = ", need_diff, "; selected d = ", d_use)
    }
  }

  tried_env <- new.env(parent = emptyenv())

  fit_order <- function(order) {
    order <- as.integer(order)
    key <- paste(order, collapse = ",")

    if (exists(key, envir = tried_env, inherits = FALSE)) {
      return(get(key, envir = tried_env, inherits = FALSE))
    }

    if (length(order) != 3L || any(order < 0L) ||
        order[1] > max.p || order[2] > max.d || order[3] > max.q) {
      assign(key, NULL, envir = tried_env)
      return(NULL)
    }

    fit <- try(
      optim_intervention_arimax(y, order, xtransf, transfer),
      silent = TRUE
    )

    if (inherits(fit, "try-error") || is.null(fit)) {
      assign(key, NULL, envir = tried_env)
      return(NULL)
    }

    assign(key, fit, envir = tried_env)
    fit
  }

  candidate_orders <- list(
    c(1L, d_use, 0L),
    c(0L, d_use, 1L),
    c(1L, d_use, 1L),
    c(2L, d_use, 0L),
    c(0L, d_use, 2L),
    c(2L, d_use, 1L),
    c(1L, d_use, 2L)
  )

  current_fit <- NULL
  current_order <- NULL
  current_score <- Inf

  for (order in candidate_orders) {
    fit <- fit_order(order)

    if (is.null(fit)) {
      next
    }

    score <- get_model_ic(fit, ic = ic, n = n)

    if (trace) {
      message("Initial candidate: ", paste(order, collapse = ","), "; score = ", score)
    }

    if (score < current_score) {
      current_fit <- fit
      current_order <- order
      current_score <- score
    }
  }

  if (is.null(current_fit)) {
    if (trace) {
      message("No valid model was found among the initial candidates.")
    }

    return(NULL)
  }

  repeat {
    p <- current_order[1]
    d <- current_order[2]
    q <- current_order[3]

    neighbor_orders <- list(
      c(p + 1L, d, q),
      c(p - 1L, d, q),
      c(p, d, q + 1L),
      c(p, d, q - 1L),
      c(p + 1L, d, q + 1L),
      c(p + 1L, d, q - 1L),
      c(p - 1L, d, q + 1L),
      c(p - 1L, d, q - 1L)
    )

    best_local_fit <- current_fit
    best_local_order <- current_order
    best_local_score <- current_score

    for (order in neighbor_orders) {
      order <- as.integer(order)

      if (order[1] < 0L || order[3] < 0L || order[1] > max.p || order[3] > max.q) {
        next
      }

      if (order[2] != d_use) {
        next
      }

      fit <- fit_order(order)

      if (is.null(fit)) {
        next
      }

      score <- get_model_ic(fit, ic = ic, n = n)

      if (trace) {
        message("Neighbor: ", paste(order, collapse = ","), "; score = ", score)
      }

      if (score < best_local_score) {
        best_local_fit <- fit
        best_local_order <- order
        best_local_score <- score
      }
    }

    if (all(best_local_order == current_order)) {
      if (trace) {
        message("No better neighboring order was found. The search is stopped.")
      }

      break
    }

    current_fit <- best_local_fit
    current_order <- best_local_order
    current_score <- best_local_score
  }

  if (!is.null(current_order) && current_order[2] == 1L) {
    alternative_order <- c(current_order[1], 0L, current_order[3])
    alternative_fit <- fit_order(alternative_order)

    if (!is.null(alternative_fit)) {
      alternative_score <- get_model_ic(alternative_fit, ic = ic, n = n)

      if (alternative_score <= current_score) {
        current_fit <- alternative_fit
        current_order <- alternative_order
        current_score <- alternative_score
      }
    }
  }

  current_fit
}

#' Fit one candidate interventional ARIMA model
#'
#' @param candidate A candidate specification.
#' @param candidate_key Candidate identifier.
#' @param y A numeric outcome time series.
#' @param lb_lag Lag used in the Ljung-Box test.
#'
#' @return A list containing the fitted result and diagnostics.
#' @keywords internal
fit_one_candidate <- function(candidate, candidate_key, y, lb_lag = 12L) {
  xdata_use <- if (is.data.frame(candidate$xdata)) {
    candidate$xdata
  } else {
    data.frame(x = as.numeric(candidate$xdata))
  }

  if (candidate$xdatab > 0L) {
    x_name <- names(xdata_use)[1]
    x_vector <- lag0_fill(xdata_use[[1]], candidate$xdatab)
    xdata_use <- data.frame(x = x_vector)
    names(xdata_use) <- x_name
  }

  best_model <- try(
    search_arimax_order(y, xdata_use, list(c(candidate$xdatar, candidate$xdatah))),
    silent = TRUE
  )

  if (inherits(best_model, "try-error") || is.null(best_model)) {
    return(list(key = candidate_key, success = FALSE))
  }

  coef_table <- try(lmtest::coeftest(best_model), silent = TRUE)

  if (inherits(coef_table, "try-error")) {
    return(list(
      key = candidate_key,
      success = FALSE,
      best_model = best_model,
      aic = best_model$aic
    ))
  }

  coef_table <- coef_table[
    !rownames(coef_table) %in% "intercept" &
      !grepl("^(ar|ma)", rownames(coef_table), ignore.case = TRUE),
    ,
    drop = FALSE
  ]

  lb_pvalue <- tryCatch(
    stats::Box.test(stats::residuals(best_model), lag = lb_lag, type = "Ljung-Box")$p.value,
    error = function(e) NA_real_
  )

  resid_indep <- !is.na(lb_pvalue) && lb_pvalue > 0.05

  list(
    key = candidate_key,
    candidate = candidate,
    success = TRUE,
    best_model = best_model,
    coef_table = coef_table,
    aic = best_model$aic,
    lb_pvalue = lb_pvalue,
    resid_indep = resid_indep
  )
}

#' Evaluate candidate models serially
#'
#' @param y A numeric outcome time series.
#' @param unevaluated_list List of unevaluated candidates.
#' @param evaluated_list List of evaluated candidates.
#' @param x_p Pulse intervention variable.
#' @param x_s Step intervention variable.
#' @param x_r Ramp intervention variable.
#'
#' @return A list containing evaluated and unevaluated candidate lists.
#' @keywords internal
evaluate_models_serial <- function(y, unevaluated_list, evaluated_list, x_p, x_s, x_r) {
  while (length(unevaluated_list) > 0L) {
    batch_list <- unevaluated_list
    batch_keys <- names(batch_list)
    unevaluated_list <- list()

    for (ii in seq_along(batch_list)) {
      result_i <- fit_one_candidate(batch_list[[ii]], batch_keys[ii], y)

      if (!result_i$success) {
        next
      }

      candidate_key <- result_i$key
      candidate <- result_i$candidate
      coef_table <- result_i$coef_table

      evaluated_list[[candidate_key]] <- list(
        aic = result_i$aic,
        result = data.frame(coef_table),
        lb_pvalue = result_i$lb_pvalue,
        resid_indep = result_i$resid_indep,
        best_model = result_i$best_model,
        candidate = candidate
      )

      if ("Pt-AR1" %in% rownames(coef_table)) {
        pt_ar1_est <- coef_table["Pt-AR1", "Estimate"]
        pt_ar1_p <- coef_table["Pt-AR1", "Pr(>|z|)"]

        if (pt_ar1_est > 0.9 || pt_ar1_est < 0) {
          unevaluated_list <- add_unevaluated_safe(
            paste0("0", candidate$xdatab, candidate$xdatah, "Step"),
            data.frame(St = x_s), 0L, candidate$xdatab, candidate$xdatah, unevaluated_list
          )
          unevaluated_list <- add_unevaluated_safe(
            paste0("1", candidate$xdatab, candidate$xdatah, "StepD"),
            data.frame(St = x_s), 1L, candidate$xdatab, candidate$xdatah, unevaluated_list
          )
        }

        if (is.nan(pt_ar1_p) || pt_ar1_p >= 0.05 || pt_ar1_est < 0) {
          unevaluated_list <- add_unevaluated_safe(
            paste0("0", candidate$xdatab, candidate$xdatah, "Pulse"),
            data.frame(Pt = x_p), 0L, candidate$xdatab, candidate$xdatah, unevaluated_list
          )
        }
      }

      if ("Pt-MA1" %in% rownames(coef_table)) {
        pt_ma1_p <- coef_table["Pt-MA1", "Pr(>|z|)"]

        if (is.nan(pt_ma1_p) || pt_ma1_p >= 0.05) {
          suffix <- sub("^[0-9]+", "", candidate_key)
          unevaluated_list <- add_unevaluated_safe(
            paste0(candidate$xdatar, candidate$xdatab, 0L, suffix),
            data.frame(Pt = x_p), candidate$xdatar, candidate$xdatab, 0L, unevaluated_list
          )
        }
      }

      if ("St-MA1" %in% rownames(coef_table)) {
        st_ma1_p <- coef_table["St-MA1", "Pr(>|z|)"]

        if (is.nan(st_ma1_p) || st_ma1_p >= 0.05) {
          suffix <- sub("^[0-9]+", "", candidate_key)
          unevaluated_list <- add_unevaluated_safe(
            paste0(candidate$xdatar, candidate$xdatab, 0L, suffix),
            data.frame(St = x_s), candidate$xdatar, candidate$xdatab, 0L, unevaluated_list
          )
        }
      }

      if ("St-AR1" %in% rownames(coef_table)) {
        st_ar1_est <- coef_table["St-AR1", "Estimate"]
        st_ar1_p <- coef_table["St-AR1", "Pr(>|z|)"]

        if (st_ar1_est > 0.9 || st_ar1_est < 0) {
          unevaluated_list <- add_unevaluated_safe(
            paste0("0", candidate$xdatab, candidate$xdatah, "Ramp"),
            data.frame(Rt = x_r), 0L, candidate$xdatab, candidate$xdatah, unevaluated_list
          )
        }

        if (is.nan(st_ar1_p) || st_ar1_p >= 0.05 || st_ar1_est < 0) {
          unevaluated_list <- add_unevaluated_safe(
            paste0("0", candidate$xdatab, candidate$xdatah, "Step"),
            data.frame(St = x_s), 0L, candidate$xdatab, candidate$xdatah, unevaluated_list
          )
        }
      }
    }
  }

  list(
    evaluated_list = evaluated_list,
    unevaluated_list = unevaluated_list
  )
}

#' Select the best candidate model by AIC under adequacy constraints
#'
#' @param evaluated_list A list of evaluated candidate models.
#'
#' @return A list containing the selected model name and AIC value.
#' @keywords internal
select_min_aic <- function(evaluated_list) {
  min_aic <- Inf
  min_name <- NULL

  for (candidate_name in names(evaluated_list)) {
    candidate_result <- evaluated_list[[candidate_name]]

    if (is.null(candidate_result$resid_indep) || !isTRUE(candidate_result$resid_indep)) {
      next
    }

    if (is.null(candidate_result$aic)) {
      next
    }

    current_aic <- candidate_result$aic

    if (current_aic >= min_aic) {
      next
    }

    coef_table <- candidate_result$result
    colnames(coef_table) <- c("Coefficient", "Std_Error", "Z_Value", "P_Value")

    pass_shape_rule <- TRUE

    if ("Pt-AR1" %in% rownames(coef_table)) {
      pt_ar1_p <- coef_table["Pt-AR1", "P_Value"]
      pt_ar1_est <- coef_table["Pt-AR1", "Coefficient"]
      pass_shape_rule <- pass_shape_rule && !is.nan(pt_ar1_p) && pt_ar1_p < 0.05 &&
        pt_ar1_est > 0 && pt_ar1_est < 1
    }

    if ("St-AR1" %in% rownames(coef_table)) {
      st_ar1_p <- coef_table["St-AR1", "P_Value"]
      st_ar1_est <- coef_table["St-AR1", "Coefficient"]
      pass_shape_rule <- pass_shape_rule && !is.nan(st_ar1_p) && st_ar1_p < 0.05 &&
        st_ar1_est > 0 && st_ar1_est < 1
    }

    if ("Pt-MA1" %in% rownames(coef_table) && "Pt-MA0" %in% rownames(coef_table)) {
      pt_ma1_p <- coef_table["Pt-MA1", "P_Value"]
      pt_ma1_est <- coef_table["Pt-MA1", "Coefficient"]
      pt_ma0_est <- coef_table["Pt-MA0", "Coefficient"]
      pass_shape_rule <- pass_shape_rule && !is.nan(pt_ma1_p) && pt_ma1_p < 0.05 &&
        pt_ma1_est * pt_ma0_est > 0
    }

    if ("St-MA1" %in% rownames(coef_table) && "St-MA0" %in% rownames(coef_table)) {
      st_ma1_p <- coef_table["St-MA1", "P_Value"]
      st_ma1_est <- coef_table["St-MA1", "Coefficient"]
      st_ma0_est <- coef_table["St-MA0", "Coefficient"]
      pass_shape_rule <- pass_shape_rule && !is.nan(st_ma1_p) && st_ma1_p < 0.05 &&
        st_ma1_est * st_ma0_est > 0
    }

    if ("Rt-MA1" %in% rownames(coef_table)) {
      pass_shape_rule <- FALSE
    }

    if (pass_shape_rule) {
      min_aic <- current_aic
      min_name <- candidate_name
    }
  }

  list(model_name = min_name, aic = min_aic)
}

#' Calculate the automatic impact-pattern identification result
#'
#' This function constructs candidate transfer functions, estimates candidate
#' interventional ARIMA models, evaluates residual independence and impact-shape
#' adequacy, and selects the optimal model by AIC among adequate candidates.
#'
#' @param y A numeric outcome time series.
#' @param x_p A pulse intervention variable.
#' @param x_s A step intervention variable.
#' @param x_r A ramp intervention variable.
#' @param k Integer. Number of lagged intervention terms used in the initial LTF model.
#'
#' @return A list containing the selected model, selected model name, AIC value,
#'   all evaluated candidate models, and candidate delay values.
#' @export
calculate_result <- function(y, x_p, x_s, x_r, k = 11L) {
  y <- as.numeric(y)
  x_p <- as.numeric(x_p)
  x_s <- as.numeric(x_s)
  x_r <- as.numeric(x_r)

  if (length(y) != length(x_p) || length(y) != length(x_s) || length(y) != length(x_r)) {
    stop("`y`, `x_p`, `x_s`, and `x_r` must have the same length.", call. = FALSE)
  }

  evaluated_list <- list()
  unevaluated_list <- list()

  initial_order <- c(1L, 0L, 0L)
  ltf_model <- fit_ltf_model(y, x_p, initial_order, k)

  continue_loop <- TRUE
  loop_index <- 1L

  if (!is.null(ltf_model)) {
    while (continue_loop && loop_index < 6L) {
      ltf_test <- stats::Box.test(ltf_model$residuals, lag = 20L, type = "Ljung-Box")

      if (ltf_test$p.value < 0.05) {
        auto_result <- forecast::auto.arima(
          ltf_model$residuals,
          max.p = 2L,
          max.q = 2L,
          max.d = 1L,
          stepwise = TRUE,
          approximation = TRUE
        )

        residual_order <- extract_arima_order(auto_result)

        if (identical(as.numeric(residual_order), c(0, 0, 0))) {
          continue_loop <- FALSE
        } else {
          ltf_model <- fit_ltf_model(y, x_p, residual_order, k)
        }

        loop_index <- loop_index + 1L
      } else {
        continue_loop <- FALSE
      }
    }
  }

  candidate_b <- identify_candidate_b(ltf_model, k = k)
  candidate_grid <- expand.grid(b = candidate_b$b, h = c(1L, 0L))
  candidate_fits <- list()

  for (i in seq_len(nrow(candidate_grid))) {
    candidate_b_i <- candidate_grid[i, "b"]
    candidate_h_i <- candidate_grid[i, "h"]

    x_p_tmp <- data.frame(Pt = x_p)

    if (candidate_b_i > 0L) {
      x_p_tmp[[1]] <- lag0_fill(x_p_tmp[[1]], candidate_b_i)
    }

    candidate_fit <- try(
      search_arimax_order(
        y = y,
        xtransf = x_p_tmp,
        transfer = list(c(1L, candidate_h_i))
      ),
      silent = TRUE
    )

    if (!inherits(candidate_fit, "try-error") && !is.null(candidate_fit)) {
      candidate_fits[[length(candidate_fits) + 1L]] <- list(
        fit = candidate_fit,
        b = candidate_b_i,
        h = candidate_h_i,
        aic = candidate_fit$aic
      )

      candidate_name <- paste0("1", candidate_b_i, candidate_h_i, "PulseD")
      unevaluated_list <- add_unevaluated(
        key = candidate_name,
        xdata = data.frame(Pt = x_p),
        xdatar = 1L,
        xdatab = candidate_b_i,
        xdatah = candidate_h_i,
        unevaluated_list = unevaluated_list
      )
    }
  }

  if (length(candidate_fits) > 0L) {
    aic_values <- vapply(candidate_fits, function(x) x$aic, numeric(1))
    best_idx <- which.min(aic_values)
    best_b <- candidate_fits[[best_idx]]$b
    best_h <- candidate_fits[[best_idx]]$h
    best_fit <- candidate_fits[[best_idx]]$fit
    candidate_name <- paste0("1", best_b, best_h, "PulseD")

    coef_best <- try(lmtest::coeftest(best_fit), silent = TRUE)

    if (!inherits(coef_best, "try-error")) {
      coef_best <- coef_best[
        !rownames(coef_best) %in% "intercept" &
          !grepl("^(ar|ma)", rownames(coef_best), ignore.case = TRUE),
        ,
        drop = FALSE
      ]

      lb_pvalue_best <- tryCatch(
        stats::Box.test(stats::residuals(best_fit), lag = 12L, type = "Ljung-Box")$p.value,
        error = function(e) NA_real_
      )

      resid_indep_best <- !is.na(lb_pvalue_best) && lb_pvalue_best > 0.05

      evaluated_list[[candidate_name]] <- list(
        aic = best_fit$aic,
        result = data.frame(coef_best),
        lb_pvalue = lb_pvalue_best,
        resid_indep = resid_indep_best,
        best_model = best_fit,
        candidate = list(xdata = data.frame(Pt = x_p), xdatar = 1L, xdatab = best_b, xdatah = best_h)
      )
    }
  } else {
    unevaluated_list <- add_unevaluated(
      key = "1_0_1_PulseD",
      xdata = data.frame(Pt = x_p),
      xdatar = 1L,
      xdatab = 0L,
      xdatah = 1L,
      unevaluated_list = unevaluated_list
    )
  }

  evaluated <- evaluate_models_serial(
    y = y,
    unevaluated_list = unevaluated_list,
    evaluated_list = evaluated_list,
    x_p = x_p,
    x_s = x_s,
    x_r = x_r
  )

  evaluated_list <- evaluated$evaluated_list

  if (length(evaluated_list) == 0L) {
    return(list(
      selected_model_name = NA_character_,
      selected_aic = NA_real_,
      selected_model = NULL,
      evaluated_models = evaluated_list,
      candidate_b = candidate_b$b
    ))
  }

  selected <- select_min_aic(evaluated_list)

  if (is.null(selected$model_name) || !is.finite(selected$aic)) {
    selected_model <- NULL
    selected_model_name <- NA_character_
    selected_aic <- NA_real_
  } else {
    selected_model <- evaluated_list[[selected$model_name]]
    selected_model_name <- selected$model_name
    selected_aic <- selected$aic
  }

  list(
    selected_model_name = selected_model_name,
    selected_aic = selected_aic,
    selected_model = selected_model,
    evaluated_models = evaluated_list,
    candidate_b = candidate_b$b
  )
}

#' Automatically identify the intervention impact pattern
#'
#' This is the main user-facing function. It constructs pulse, step, and ramp
#' intervention variables from the intervention time point, then calls
#' `calculate_result()` to select the optimal interventional ARIMA model.
#'
#' @param y A numeric outcome time series.
#' @param intervention_point Integer. The intervention time point using 1-based indexing.
#' @param k Integer. Number of lagged intervention terms used in the initial LTF model.
#'
#' @return A list containing the selected model, selected model name, AIC value,
#'   all evaluated models, and candidate delay values.
#' @export
#'
#' @examples


auto_identify_impact_pattern <- function(y, intervention_point, k = 11L) {
  y <- as.numeric(y)
  n <- length(y)
  intervention_point <- as.integer(intervention_point)

  if (intervention_point < 1L || intervention_point > n) {
    stop("`intervention_point` must be between 1 and length(y).", call. = FALSE)
  }

  pre_n <- intervention_point - 1L
  post_n <- n - intervention_point + 1L

  x_p <- c(rep(0, pre_n), 1, rep(0, post_n - 1L))
  x_s <- c(rep(0, pre_n), rep(1, post_n))
  x_r <- c(rep(0, pre_n), seq_len(post_n))

  calculate_result(
    y = y,
    x_p = x_p,
    x_s = x_s,
    x_r = x_r,
    k = k
  )
}

#' y <- c(
#'   4.9, 5.2, 5.7, 5.7, 6.2, 6.7, 6.9, 7.1, 6.6, 7.0,
#'   6.9, 6.4, 6.6, 6.4, 7.0, 7.3, 6.0, 6.3, 4.8, 5.3,
#'   5.4, 4.7, 4.9, 4.4, 5.1, 5.3, 6.0, 5.9, 5.9, 5.6,
#'   5.3, 4.5, 4.7, 4.6, 4.3, 5.0, 5.2, 6.2, 5.8, 6.7,
#'   5.7, 6.1, 7.2, 6.5, 6.1, 6.3, 6.4, 7.0, 7.6, 7.2,
#'   7.5, 7.8, 7.2, 7.5, 5.6, 5.7, 4.9, 5.1, 6.2, 6.0,
#'   6.1, 7.5, 7.8, 8.0, 8.0, 8.1, 7.6, 7.1, 6.6, 5.6,
#'   5.9, 6.6, 6.8, 7.8, 7.9, 8.7, 7.7, 7.3, 6.7, 7.5,
#'   6.4, 9.7, 7.5, 7.1, 6.4, 6.0, 5.7, 5.0, 4.2, 5.1,
#'   5.4, 5.1, 5.3, 5.0, 4.8, 4.7, 5.0, 5.4, 4.3, 3.5
#' )
#'
#' result <- auto_identify_impact_pattern(y,  intervention_point = 82, k = 11)
#' result$selected_model$best_model 
#' 
