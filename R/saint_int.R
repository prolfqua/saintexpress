.saint_int_read_tables <- function(si) {
  tables <- .saint_int_prepare_tables(si)
  inter <- tables$inter
  prey <- tables$prey
  bait <- .saint_int_index_bait_replicates(tables$bait)

  inter <- .saint_int_index_interactions(inter, prey, bait)
  ubait <- sort(unique(bait$baitId[!bait$is_ctrl]))
  ip_idx_to_bait_no <- match(bait$baitId[!bait$is_ctrl], ubait)

  matrices <- .saint_int_build_observation_matrices(inter, prey, bait)
  ui <- .saint_int_unique_test_interactions(inter, prey, ip_idx_to_bait_no)
  matrices$test_mat <- .saint_int_mask_self_interactions(matrices$test_mat, ui)

  list(
    inter = inter,
    prey = prey,
    bait = bait,
    ubait = ubait,
    ip_idx_to_bait_no = ip_idx_to_bait_no,
    test_mat = matrices$test_mat,
    ctrl_mat = matrices$ctrl_mat,
    ui = ui
  )
}

.saint_int_prepare_tables <- function(si) {
  stopifnot(names(si) == c("inter", "prey", "bait"))

  inter <- as.data.frame(si$inter, stringsAsFactors = FALSE)
  prey <- as.data.frame(si$prey, stringsAsFactors = FALSE)
  bait <- as.data.frame(si$bait, stringsAsFactors = FALSE)

  names(inter)[seq_len(min(4, ncol(inter)))] <- c(
    "ipId",
    "baitId",
    "preyId",
    "quant"
  )[seq_len(min(4, ncol(inter)))]
  names(bait)[seq_len(min(3, ncol(bait)))] <- c(
    "ipId",
    "baitId",
    "CorT"
  )[seq_len(min(3, ncol(bait)))]
  prey <- .saint_int_normalize_prey_table(prey)

  if (!all(c("ipId", "baitId", "preyId", "quant") %in% names(inter))) {
    stop("inter must have columns ipId, baitId, preyId, quant")
  }
  if (!all(c("preyId", "preyGeneId") %in% names(prey))) {
    stop("prey must have columns preyId and preyGeneId")
  }
  if (!all(c("ipId", "baitId", "CorT") %in% names(bait))) {
    stop("bait must have columns ipId, baitId, CorT")
  }
  if (anyDuplicated(prey$preyId)) {
    stop("duplicate preys in prey file")
  }
  if (!all(bait$CorT %in% c("T", "C"))) {
    stop("3rd column of bait file must be 'T' or 'C'")
  }

  prey$row <- seq_len(nrow(prey))
  list(inter = inter, prey = prey, bait = bait)
}

.saint_int_normalize_prey_table <- function(prey) {
  if (all(c("preyId", "preyGeneId") %in% names(prey))) {
    return(prey)
  }

  prey_raw <- prey
  prey <- data.frame(
    preyId = as.character(prey_raw[[1]]),
    preyLength = 0,
    preyGeneId = as.character(prey_raw[[1]]),
    stringsAsFactors = FALSE
  )
  if (ncol(prey_raw) == 2) {
    prey_length <- .saint_int_parse_length_column(prey_raw[[2]])
    if (prey_length$is_length) {
      prey$preyLength <- prey_length$values
    } else {
      prey$preyGeneId <- as.character(prey_raw[[2]])
    }
  } else if (ncol(prey_raw) > 2) {
    prey$preyLength <- as.numeric(prey_raw[[2]])
    prey$preyGeneId <- as.character(prey_raw[[3]])
  }
  prey
}

.saint_int_parse_length_column <- function(x) {
  if (is.numeric(x)) {
    values <- x
    return(list(is_length = !anyNA(values), values = values))
  }

  numeric_pattern <- paste0(
    "^[-+]?(?:\\d+\\.?\\d*|\\.\\d+)",
    "(?:[eE][-+]?\\d+)?$"
  )
  x <- trimws(as.character(x))
  is_length <- all(nzchar(x)) && all(grepl(numeric_pattern, x))
  values <- if (is_length) as.numeric(x) else numeric()
  list(is_length = is_length, values = values)
}

.saint_int_index_bait_replicates <- function(bait) {
  n_test_seen <- 0L
  n_ctrl_seen <- 0L
  bait$col <- integer(nrow(bait))
  bait$is_ctrl <- bait$CorT != "T"
  for (bait_row in seq_len(nrow(bait))) {
    if (bait$is_ctrl[[bait_row]]) {
      n_ctrl_seen <- n_ctrl_seen + 1L
      bait$col[[bait_row]] <- n_ctrl_seen
    } else {
      n_test_seen <- n_test_seen + 1L
      bait$col[[bait_row]] <- n_test_seen
    }
  }
  bait
}

.saint_int_index_interactions <- function(inter, prey, bait) {
  inter$quant <- log(as.numeric(inter$quant))
  inter$row <- prey$row[match(inter$preyId, prey$preyId)]
  bait_match <- match(inter$ipId, bait$ipId)
  inter$col <- bait$col[bait_match]
  inter$is_ctrl <- bait$is_ctrl[bait_match]
  missing_prey <- is.na(inter$row)
  missing_bait <- is.na(inter$col)
  if (any(missing_prey)) {
    stop("prey ", inter$preyId[which(missing_prey)[[1]]], " not found")
  }
  if (any(missing_bait)) {
    stop("bait not found")
  }
  inter
}

.saint_int_build_observation_matrices <- function(inter, prey, bait) {
  nprey <- nrow(prey)
  n_test_seen <- max(c(0L, bait$col[!bait$is_ctrl]))
  n_ctrl_seen <- max(c(0L, bait$col[bait$is_ctrl]))
  test_mat <- matrix(NA_real_, nrow = nprey, ncol = n_test_seen)
  ctrl_mat <- matrix(NA_real_, nrow = nprey, ncol = n_ctrl_seen)
  for (interaction_row in seq_len(nrow(inter))) {
    prey_row <- inter$row[[interaction_row]]
    replicate_col <- inter$col[[interaction_row]]
    quant <- inter$quant[[interaction_row]]
    if (isTRUE(inter$is_ctrl[[interaction_row]])) {
      ctrl_mat[prey_row, replicate_col] <- quant
    } else {
      test_mat[prey_row, replicate_col] <- quant
    }
  }
  list(test_mat = test_mat, ctrl_mat = ctrl_mat)
}

.saint_int_unique_test_interactions <- function(
  inter,
  prey,
  ip_idx_to_bait_no
) {
  ui_key <- new.env(parent = emptyenv())
  ui <- vector("list", 0)
  test_inter <- inter[!inter$is_ctrl, , drop = FALSE]
  if (nrow(test_inter) > 0) {
    for (interaction_row in seq_len(nrow(test_inter))) {
      bait_col <- ip_idx_to_bait_no[[test_inter$col[[interaction_row]]]]
      key <- paste(test_inter$row[[interaction_row]], bait_col, sep = "\r")
      idx <- ui_key[[key]]
      if (is.null(idx)) {
        idx <- length(ui) + 1L
        ui_key[[key]] <- idx
        ui[[idx]] <- list(
          baitId = test_inter$baitId[[interaction_row]],
          preyId = test_inter$preyId[[interaction_row]],
          preyGeneId = prey$preyGeneId[[test_inter$row[[interaction_row]]]],
          row = test_inter$row[[interaction_row]],
          baitcol = bait_col,
          colIds = integer()
        )
      }
      ui[[idx]]$colIds <- c(ui[[idx]]$colIds, test_inter$col[[interaction_row]])
    }
  }
  ui
}

.saint_int_mask_self_interactions <- function(test_mat, ui) {
  for (entry in ui) {
    if (identical(entry$preyGeneId, entry$baitId)) {
      test_mat[entry$row, entry$colIds] <- NA_real_
    }
  }
  test_mat
}

.saint_int_normalize <- function(parsed) {
  vals <- c(as.vector(parsed$test_mat), as.vector(parsed$ctrl_mat))
  vals <- vals[is.finite(vals) & vals > -1e5]
  if (length(vals) < 2) {
    stop(
      "not enough finite intensity values for SAINTexpress-int normalization"
    )
  }
  center <- mean(vals)
  scale <- sqrt(stats::var(vals) * (length(vals) - 1) / length(vals))
  parsed$test_mat <- (parsed$test_mat - center) / scale
  parsed$ctrl_mat <- (parsed$ctrl_mat - center) / scale
  parsed
}

.saint_int_var_mle <- function(x) {
  if (length(x) <= 1) {
    return(0)
  }
  mean((x - mean(x))^2)
}

.saint_int_var_sample <- function(x) {
  if (length(x) <= 1) {
    stop("sample variance from n<2")
  }
  stats::var(x)
}

.saint_int_mean_sd_estimate <- function(data, t, sd_NA, lse) {
  detected_values <- data[!is.na(data)]
  detected <- length(detected_values)
  if (detected == 0) {
    return(c(mean = t, sd = sd_NA))
  }
  mu <- mean(detected_values)
  if (detected >= 3) {
    return(c(mean = mu, sd = sqrt(.saint_int_var_mle(detected_values))))
  }
  c(mean = mu, sd = exp(mu * lse[[2]] + lse[[1]]))
}

.saint_int_quant_log_pdf <- function(x, mu = 0, sigma = 1) {
  if (is.na(x)) {
    stop("isnan")
  }
  if (abs(x - mu) / sigma > 5) {
    x <- (if (x < 0) -1 else 1) * 5 * sigma + mu
  }
  stats::dnorm(x, mean = mu, sd = sigma, log = TRUE)
}

.saint_int_build_model <- function(parsed, p2p_mapping = NULL, L = 100L) {
  if (L <= 0) {
    stop("L is must be a positive integer.")
  }

  dimensions <- .saint_int_model_dimensions(parsed)
  replicates <- .saint_int_replicate_layout(parsed, dimensions$nbait)
  controls <- .saint_int_estimate_control_model(parsed$ctrl_mat)
  tests <- .saint_int_estimate_test_model(
    parsed$test_mat,
    replicates$bait_no_to_ip_idxes,
    controls
  )
  p2p_mapping <- .saint_int_default_p2p_mapping(p2p_mapping, dimensions$nprey)

  .saint_int_new_model(
    parsed,
    dimensions,
    replicates,
    controls,
    tests,
    p2p_mapping
  )
}

.saint_int_model_dimensions <- function(parsed) {
  list(
    nprey = nrow(parsed$prey),
    nbait = length(parsed$ubait),
    n_ctrl_ip = ncol(parsed$ctrl_mat)
  )
}

.saint_int_replicate_layout <- function(parsed, nbait) {
  bait_no_to_ip_idxes <- .saint_int_bait_replicate_indices(
    parsed$ip_idx_to_bait_no,
    nbait
  )
  list(
    bait_no_to_ip_idxes = bait_no_to_ip_idxes,
    n_rep_vec = vapply(bait_no_to_ip_idxes, length, integer(1))
  )
}

.saint_int_estimate_control_model <- function(ctrl_mat) {
  t_value <- .saint_int_missing_threshold(ctrl_mat)
  ctrl_lse <- .saint_int_control_lse(ctrl_mat)
  ctrl_priors <- .saint_int_control_priors(
    ctrl_mat,
    t_value,
    ctrl_lse$sd_NA,
    ctrl_lse$lse
  )
  list(
    t = t_value,
    eta = ctrl_priors$eta,
    sd_false = ctrl_priors$sd_ctrl,
    sd_NA = ctrl_lse$sd_NA
  )
}

.saint_int_estimate_test_model <- function(
  test_mat_DATA,
  bait_no_to_ip_idxes,
  controls
) {
  test_mat <- .saint_int_test_replicate_matrix(
    test_mat_DATA,
    bait_no_to_ip_idxes
  )
  d <- .saint_int_signal_shift(test_mat_DATA, controls$eta)
  Z <- .saint_int_initial_z(test_mat, controls$eta)
  list(
    test_mat = test_mat,
    test_mat1 = test_mat,
    d = d,
    Z = Z,
    sd_true = .saint_int_true_sd(test_mat, Z, controls$sd_false, controls$sd_NA)
  )
}

.saint_int_default_p2p_mapping <- function(p2p_mapping, nprey) {
  if (is.null(p2p_mapping)) {
    return(rep(list(integer()), nprey))
  }
  p2p_mapping
}

.saint_int_new_model <- function(
  parsed,
  dimensions,
  replicates,
  controls,
  tests,
  p2p_mapping
) {
  list(
    nprey = dimensions$nprey,
    nbait = dimensions$nbait,
    test_mat_DATA = parsed$test_mat,
    test_mat = tests$test_mat,
    test_mat1 = tests$test_mat1,
    n_rep_vec = replicates$n_rep_vec,
    p2p_mapping = p2p_mapping,
    t = controls$t,
    n_ctrl_ip = dimensions$n_ctrl_ip,
    Z = tests$Z,
    beta0 = 0,
    beta1 = 0,
    gamma = 0,
    eta = controls$eta,
    d = tests$d,
    sd_true = tests$sd_true,
    sd_false = controls$sd_false
  )
}

.saint_int_bait_replicate_indices <- function(ip_idx_to_bait_no, nbait) {
  lapply(seq_len(nbait), function(bait_col) {
    which(ip_idx_to_bait_no == bait_col)
  })
}

.saint_int_missing_threshold <- function(ctrl_mat) {
  controls_with_missing <- numeric()
  for (prey_row in seq_len(nrow(ctrl_mat))) {
    ctrl <- ctrl_mat[prey_row, ]
    if (anyNA(ctrl)) {
      controls_with_missing <- c(controls_with_missing, ctrl[!is.na(ctrl)])
    }
  }

  threshold <- mean(controls_with_missing)
  if (!is.finite(threshold)) {
    threshold <- mean(ctrl_mat[is.finite(ctrl_mat)])
  }
  threshold
}

.saint_int_control_lse <- function(ctrl_mat) {
  complete_ctrl <- ctrl_mat[stats::complete.cases(ctrl_mat), , drop = FALSE]
  log_sd_ctrls <- apply(complete_ctrl, 1, function(ctrl) {
    log(sqrt(.saint_int_var_mle(ctrl)))
  })
  mean_ctrls <- rowMeans(complete_ctrl)

  if (!length(log_sd_ctrls)) {
    return(list(lse = c(0, 0), sd_NA = 1))
  }

  sorted_log_sd_ctrls <- sort(log_sd_ctrls)
  median_sd <- exp(sorted_log_sd_ctrls[[
    length(sorted_log_sd_ctrls) %/% 2 + 1L
  ]])

  if (length(log_sd_ctrls) < 2 || stats::var(mean_ctrls) == 0) {
    return(list(
      lse = c(log(stats::median(exp(log_sd_ctrls), na.rm = TRUE)), 0),
      sd_NA = median_sd
    ))
  }

  fit <- stats::lm(log_sd_ctrls ~ mean_ctrls)
  list(lse = unname(stats::coef(fit)), sd_NA = median_sd)
}

.saint_int_control_priors <- function(ctrl_mat, t_value, sd_NA, lse) {
  eta <- numeric(nrow(ctrl_mat))
  sd_ctrl <- numeric(nrow(ctrl_mat))
  for (prey_row in seq_len(nrow(ctrl_mat))) {
    est <- .saint_int_mean_sd_estimate(
      ctrl_mat[prey_row, ],
      t_value,
      sd_NA,
      lse
    )
    eta[[prey_row]] <- est[["mean"]]
    sd_ctrl[[prey_row]] <- max(sd_NA, est[["sd"]])
  }
  list(eta = eta, sd_ctrl = sd_ctrl)
}

.saint_int_test_replicate_matrix <- function(test_mat, bait_no_to_ip_idxes) {
  nprey <- nrow(test_mat)
  nbait <- length(bait_no_to_ip_idxes)
  .saint_int_list_matrix(nprey, nbait, function(prey_row, bait_col) {
    test_mat[prey_row, bait_no_to_ip_idxes[[bait_col]]]
  })
}

.saint_int_signal_shift <- function(test_mat, eta) {
  d <- numeric(nrow(test_mat))
  for (prey_row in seq_len(nrow(test_mat))) {
    detected <- test_mat[prey_row, ]
    detected <- detected[!is.na(detected) & eta[[prey_row]] < detected]
    mu <- if (length(detected)) mean(detected) else 0
    mu <- max(mu, eta[[prey_row]])
    d[[prey_row]] <- max(mu - eta[[prey_row]], log(4))
  }
  d
}

.saint_int_initial_z <- function(test_mat, eta) {
  .saint_int_list_matrix(
    nrow(test_mat),
    ncol(test_mat),
    function(prey_row, bait_col) {
      y <- test_mat[[prey_row, bait_col]]
      ifelse(is.na(y), FALSE, y > (eta[[prey_row]] + log(5)))
    }
  )
}

.saint_int_true_sd <- function(test_mat, Z, sd_false, sd_NA) {
  sd_true <- numeric(nrow(test_mat))
  for (prey_row in seq_len(nrow(test_mat))) {
    z1 <- unlist(lapply(seq_len(ncol(test_mat)), function(bait_col) {
      test_mat[[prey_row, bait_col]][Z[[prey_row, bait_col]]]
    }))
    sd_true[[prey_row]] <- if (length(z1) >= 2) {
      sqrt(.saint_int_var_sample(z1))
    } else {
      sd_false[[prey_row]]
    }
    sd_true[[prey_row]] <- max(sd_NA, sd_true[[prey_row]])
  }
  sd_true
}

.saint_int_list_matrix <- function(nrow, ncol, fn) {
  out <- vector("list", nrow * ncol)
  dim(out) <- c(nrow, ncol)
  for (row_idx in seq_len(nrow)) {
    for (col_idx in seq_len(ncol)) {
      out[[row_idx, col_idx]] <- fn(row_idx, col_idx)
    }
  }
  out
}

.saint_int_llik_mrf_gamma0 <- function(model, beta1) {
  loglik <- 0
  mrf_true <- exp(beta1)
  mrf_false <- exp(0)
  for (prey_row in seq_len(model$nprey)) {
    for (bait_col in seq_len(model$nbait)) {
      y <- model$test_mat1[[prey_row, bait_col]]
      if (all(is.na(y))) {
        next
      }
      k <- model$Z[[prey_row, bait_col]]
      m <- model$n_rep_vec[[bait_col]]
      prod <- prod(ifelse(k, mrf_true, mrf_false) / (mrf_true + mrf_false))
      loglik <- loglik + log(prod) / m
    }
  }
  loglik
}

.saint_int_gsum_mat <- function(model) {
  gsum <- matrix(0, nrow = model$nprey, ncol = model$nbait)
  for (prey_row in seq_len(model$nprey)) {
    links <- model$p2p_mapping[[prey_row]]
    if (!length(links)) {
      next
    }
    for (bait_col in seq_len(model$nbait)) {
      gsum[prey_row, bait_col] <- sum(vapply(
        links,
        function(linked_prey) mean(model$Z[[linked_prey, bait_col]]),
        numeric(1)
      ))
    }
  }
  gsum
}

.saint_int_llik_mrf <- function(model, beta1, gamma, gsum_mat) {
  loglik <- 0
  mrf_false <- exp(0)
  for (prey_row in seq_len(model$nprey)) {
    for (bait_col in seq_len(model$nbait)) {
      y <- model$test_mat1[[prey_row, bait_col]]
      if (all(is.na(y))) {
        next
      }
      k <- model$Z[[prey_row, bait_col]]
      m <- model$n_rep_vec[[bait_col]]
      mrf_true <- exp(beta1 + gamma * gsum_mat[prey_row, bait_col])
      prod <- prod(ifelse(k, mrf_true, mrf_false) / (mrf_true + mrf_false))
      loglik <- loglik + log(prod) / m
    }
  }
  loglik
}

.saint_int_loglikelihood <- function(model) {
  loglik <- 0
  for (prey_row in seq_len(model$nprey)) {
    for (bait_col in seq_len(model$nbait)) {
      y <- model$test_mat1[[prey_row, bait_col]]
      if (all(is.na(y))) {
        next
      }
      k <- model$Z[[prey_row, bait_col]]
      m <- model$n_rep_vec[[bait_col]]
      links <- model$p2p_mapping[[prey_row]]
      gsum <- if (length(links)) {
        sum(vapply(
          links,
          function(linked_prey) mean(model$Z[[linked_prey, bait_col]]),
          numeric(1)
        ))
      } else {
        0
      }
      mrf_true <- exp(model$beta1 + model$gamma * gsum)
      mrf_false <- exp(model$beta0)
      prod_val <- 1
      for (rep in seq_len(m)) {
        y_rep <- if (is.na(y[[rep]])) model$t else y[[rep]]
        dens_true <- .saint_int_quant_log_pdf(
          min(y_rep, model$eta[[prey_row]] + model$d[[prey_row]]),
          model$eta[[prey_row]] + model$d[[prey_row]],
          model$sd_true[[prey_row]]
        )
        dens_false <- .saint_int_quant_log_pdf(
          max(y_rep, model$eta[[prey_row]]),
          model$eta[[prey_row]],
          model$sd_false[[prey_row]]
        )
        prod_val <- prod_val *
          (if (k[[rep]]) {
            mrf_true * exp(dens_true)
          } else {
            mrf_false * exp(dens_false)
          }) /
          (mrf_true + mrf_false)
      }
      loglik <- loglik + log(prod_val) / m
    }
  }
  loglik
}

.saint_int_optimize_gamma0 <- function(model, optimizer) {
  old_beta1 <- model$beta1
  oldf <- .saint_int_llik_mrf_gamma0(model, old_beta1)
  if (
    identical(optimizer, "nloptr") && requireNamespace("nloptr", quietly = TRUE)
  ) {
    res <- nloptr::nloptr(
      x0 = old_beta1,
      eval_f = function(x) -.saint_int_llik_mrf_gamma0(model, x[[1]]),
      lb = -15,
      ub = 15,
      opts = list(algorithm = "NLOPT_LN_COBYLA", ftol_abs = 1e-4, maxeval = 100)
    )
    new_beta1 <- res$solution[[1]]
    maxf <- -res$objective
  } else {
    res <- stats::optim(
      par = old_beta1,
      fn = function(x) -.saint_int_llik_mrf_gamma0(model, x[[1]]),
      method = "L-BFGS-B",
      lower = -15,
      upper = 15,
      control = list(maxit = 100)
    )
    new_beta1 <- res$par[[1]]
    maxf <- -res$value
  }
  model$beta1 <- new_beta1
  if (maxf < oldf) {
    model$beta1 <- old_beta1
  }
  model
}

.saint_int_optimize_mrf <- function(model, optimizer) {
  gsum_mat <- .saint_int_gsum_mat(model)
  old_beta1 <- model$beta1
  old_gamma <- model$gamma
  oldf <- .saint_int_llik_mrf(model, old_beta1, old_gamma, gsum_mat)
  if (
    identical(optimizer, "nloptr") && requireNamespace("nloptr", quietly = TRUE)
  ) {
    res <- nloptr::nloptr(
      x0 = c(old_beta1, old_gamma),
      eval_f = function(x) {
        -.saint_int_llik_mrf(model, x[[1]], x[[2]], gsum_mat)
      },
      lb = c(-15, 0),
      ub = c(15, 10),
      opts = list(algorithm = "NLOPT_LN_COBYLA", ftol_abs = 1e-4, maxeval = 100)
    )
    x <- res$solution
    maxf <- -res$objective
  } else {
    res <- stats::optim(
      par = c(old_beta1, old_gamma),
      fn = function(x) -.saint_int_llik_mrf(model, x[[1]], x[[2]], gsum_mat),
      method = "L-BFGS-B",
      lower = c(-15, 0),
      upper = c(15, 10),
      control = list(maxit = 100)
    )
    x <- res$par
    maxf <- -res$value
  }
  model$beta1 <- x[[1]]
  model$gamma <- x[[2]]
  if (maxf <= .saint_int_llik_mrf(model, x[[1]], 0, gsum_mat)) {
    model$gamma <- 0
  }
  if (maxf < oldf) {
    model$beta1 <- old_beta1
    model$gamma <- old_gamma
  }
  model
}

.saint_int_icm_z <- function(model) {
  pre_calc <- .saint_int_precalculate_rep_logpdf(model)

  for (prey_row in seq_len(model$nprey)) {
    for (bait_col in seq_len(model$nbait)) {
      m <- model$n_rep_vec[[bait_col]]
      for (rep in seq_len(m)) {
        first <- .saint_int_loglikelihood_z(
          model,
          prey_row,
          bait_col,
          rep,
          pre_calc
        )
        model$Z[[prey_row, bait_col]][[rep]] <- !model$Z[[prey_row, bait_col]][[
          rep
        ]]
        second <- .saint_int_loglikelihood_z(
          model,
          prey_row,
          bait_col,
          rep,
          pre_calc
        )
        if (first > second) {
          model$Z[[prey_row, bait_col]][[rep]] <- !model$Z[[
            prey_row,
            bait_col
          ]][[rep]]
        }
      }
    }
  }
  model
}

.saint_int_precalculate_rep_logpdf <- function(model) {
  .saint_int_list_matrix(
    model$nprey,
    model$nbait,
    function(prey_row, bait_col) {
      y <- model$test_mat1[[prey_row, bait_col]]
      m <- model$n_rep_vec[[bait_col]]
      vals <- matrix(0, nrow = m, ncol = 2)
      for (rep in seq_len(m)) {
        y_rep <- if (is.na(y[[rep]])) model$t else y[[rep]]
        vals[rep, 1] <- .saint_int_quant_log_pdf(
          min(y_rep, model$eta[[prey_row]] + model$d[[prey_row]]),
          model$eta[[prey_row]] + model$d[[prey_row]],
          model$sd_true[[prey_row]]
        )
        vals[rep, 2] <- .saint_int_quant_log_pdf(
          max(y_rep, model$eta[[prey_row]]),
          model$eta[[prey_row]],
          model$sd_false[[prey_row]]
        )
      }
      vals
    }
  )
}

.saint_int_loglikelihood_z <- function(
  model,
  prey_row,
  bait_col,
  rep,
  pre_calc
) {
  gsum <- 0
  if (model$gamma != 0) {
    links <- model$p2p_mapping[[prey_row]]
    if (length(links)) {
      gsum <- sum(vapply(
        links,
        function(linked_prey) mean(model$Z[[linked_prey, bait_col]]),
        numeric(1)
      ))
    }
  }
  log_mrf_true <- model$beta1 + model$gamma * gsum
  log_mrf_false <- model$beta0
  pcl <- pre_calc[[prey_row, bait_col]][rep, ]
  if (model$Z[[prey_row, bait_col]][[rep]]) {
    log_mrf_true + pcl[[1]]
  } else {
    log_mrf_false + pcl[[2]]
  }
}

.saint_int_icms <- function(model, with_gamma, optimizer) {
  newllik <- .saint_int_loglikelihood(model)
  model$gamma <- 0
  for (iter in seq_len(15)) {
    oldllik <- newllik
    model <- .saint_int_icm_z(model)
    model <- if (with_gamma) {
      .saint_int_optimize_mrf(model, optimizer)
    } else {
      .saint_int_optimize_gamma0(model, optimizer)
    }
    newllik <- .saint_int_loglikelihood(model)
    if (newllik >= oldllik && exp(newllik - oldllik) - 1 < 1e-3) {
      break
    }
  }
  model
}

.saint_int_calculate_score <- function(model, R = 100L) {
  average_score <- matrix(0, nrow = model$nprey, ncol = model$nbait)
  max_score <- matrix(0, nrow = model$nprey, ncol = model$nbait)
  min_log_odds_score <- matrix(0, nrow = model$nprey, ncol = model$nbait)

  for (prey_row in seq_len(model$nprey)) {
    for (bait_col in seq_len(model$nbait)) {
      cell_score <- .saint_int_calculate_cell_score(
        model,
        prey_row,
        bait_col,
        R
      )
      average_score[prey_row, bait_col] <- cell_score$average
      min_log_odds_score[prey_row, bait_col] <- cell_score$odds
      max_score[prey_row, bait_col] <- cell_score$maximum
    }
  }

  list(average = average_score, maximum = max_score, odds = min_log_odds_score)
}

.saint_int_calculate_cell_score <- function(model, prey_row, bait_col, R) {
  y <- model$test_mat1[[prey_row, bait_col]]
  m <- model$n_rep_vec[[bait_col]]
  log_mrf <- .saint_int_cell_log_mrf(model, prey_row, bait_col)
  rep_scores <- .saint_int_rep_scores(model, prey_row, y, log_mrf)
  tmp_scores <- sort(rep_scores$score, decreasing = TRUE)
  tmp_odds <- sort(rep_scores$odds, decreasing = TRUE)
  max_rep <- R

  list(
    average = if (m > max_rep) {
      mean(tmp_scores[seq_len(max_rep)])
    } else {
      mean(tmp_scores)
    },
    maximum = tmp_scores[[1]],
    odds = tmp_odds[[length(tmp_odds)]]
  )
}

.saint_int_cell_log_mrf <- function(model, prey_row, bait_col) {
  links <- model$p2p_mapping[[prey_row]]
  gsum <- if (length(links)) {
    sum(vapply(
      links,
      function(linked_prey) mean(model$Z[[linked_prey, bait_col]]),
      numeric(1)
    ))
  } else {
    0
  }
  c(true = model$beta1 + model$gamma * gsum, false = model$beta0)
}

.saint_int_rep_scores <- function(model, prey_row, y, log_mrf) {
  score <- numeric(length(y))
  odds <- numeric(length(y))

  for (rep in seq_along(y)) {
    y_rep <- if (is.na(y[[rep]])) model$t else y[[rep]]
    true_log <- log_mrf[["true"]] +
      .saint_int_quant_log_pdf(
        min(y_rep, model$eta[[prey_row]] + model$d[[prey_row]]),
        model$eta[[prey_row]] + model$d[[prey_row]],
        model$sd_true[[prey_row]]
      )
    false_log <- log_mrf[["false"]] +
      .saint_int_quant_log_pdf(
        max(y_rep, model$eta[[prey_row]]),
        model$eta[[prey_row]],
        model$sd_false[[prey_row]]
      )
    score[[rep]] <- .saint_int_logit_probability(true_log, false_log)
    if (is.na(y[[rep]])) {
      score[[rep]] <- 0
    }
    odds[[rep]] <- true_log - false_log
  }

  list(score = score, odds = odds)
}

.saint_int_logit_probability <- function(true_log, false_log) {
  log_scale <- max(true_log, false_log)
  denom <- log_scale +
    log(exp(true_log - log_scale) + exp(false_log - log_scale))
  exp(true_log - denom)
}

.saint_int_computation <- function(model, with_gamma, optimizer) {
  fitted <- .saint_int_icms(
    model,
    with_gamma = with_gamma,
    optimizer = optimizer
  )
  list(
    model = fitted,
    scores = .saint_int_calculate_score(
      fitted,
      R = attr(model, "R", exact = TRUE) %||% 100L
    )
  )
}

.saint_int_format_num <- function(x) {
  ifelse(is.na(x), ".", sprintf("%.3f", x))
}

.saint_int_list_output <- function(model, parsed, scores, topo_scores) {
  out <- lapply(parsed$ui, function(ui) {
    .saint_int_output_row(ui, model, parsed, scores, topo_scores)
  })
  if (!length(out)) {
    return(data.frame())
  }
  res <- do.call(rbind, out)
  numeric_cols <- c(
    "IntensitySum",
    "AvgIntensity",
    "AvgP",
    "MaxP",
    "TopoAvgP",
    "TopoMaxP",
    "SaintScore",
    "OddsScore",
    "FoldChange",
    "BFDR"
  )
  for (col in numeric_cols) {
    res[[col]] <- round(res[[col]], 3)
  }
  res
}

.saint_int_output_row <- function(ui, model, parsed, scores, topo_scores) {
  prey_row <- ui$row
  bait_col <- ui$baitcol
  quant_orig <- exp(model$test_mat[[prey_row, bait_col]])
  quant <- quant_orig
  quant[is.na(quant)] <- exp(model$t)
  avg_quant <- mean(quant)
  ctrl_counts <- exp(parsed$ctrl_mat[prey_row, ])
  avg_p <- scores$average[prey_row, bait_col]
  topo_avg_p <- topo_scores$average[prey_row, bait_col]

  data.frame(
    Bait = ui$baitId,
    Prey = ui$preyId,
    PreyGene = ui$preyGeneId,
    Intensity = paste(.saint_int_format_num(quant_orig), collapse = "|"),
    IntensitySum = sum(quant),
    AvgIntensity = avg_quant,
    NumReplicates = length(quant),
    ctrlIntensity = paste(.saint_int_format_num(ctrl_counts), collapse = "|"),
    AvgP = avg_p,
    MaxP = scores$maximum[prey_row, bait_col],
    TopoAvgP = topo_avg_p,
    TopoMaxP = topo_scores$maximum[prey_row, bait_col],
    SaintScore = max(topo_avg_p, avg_p),
    OddsScore = topo_scores$odds[prey_row, bait_col],
    FoldChange = avg_quant / exp(model$eta[[prey_row]]),
    BFDR = .saint_int_bfdr(avg_p, scores$average),
    boosted_by = .saint_int_boosted_by(model, parsed, prey_row, bait_col),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.saint_int_bfdr <- function(avg_p, average_score) {
  all_average_values <- as.vector(average_score)
  greater <- all_average_values > avg_p
  denom <- sum(greater)
  if (denom == 0) {
    return(0)
  }
  1 - sum(all_average_values[greater]) / denom
}

.saint_int_boosted_by <- function(model, parsed, prey_row, bait_col) {
  links <- model$p2p_mapping[[prey_row]]
  if (!length(links)) {
    return("")
  }
  active <- vapply(
    links,
    function(linked_prey) mean(model$Z[[linked_prey, bait_col]]) > 0,
    logical(1)
  )
  boosted <- parsed$prey$preyId[links[active]]
  if (length(boosted)) paste0(paste(boosted, collapse = "|"), "|") else ""
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

saint_int_r <- function(
  si,
  optimizer = c("base", "nloptr"),
  R = 100L,
  L = 100L,
  p2p_mapping = NULL
) {
  optimizer <- match.arg(optimizer)
  parsed <- .saint_int_read_tables(si)
  parsed <- .saint_int_normalize(parsed)
  model <- .saint_int_build_model(parsed, p2p_mapping = p2p_mapping, L = L)
  attr(model, "R") <- R

  gamma0 <- .saint_int_computation(
    model,
    with_gamma = FALSE,
    optimizer = optimizer
  )
  if (!is.null(p2p_mapping)) {
    topo <- .saint_int_computation(
      model,
      with_gamma = TRUE,
      optimizer = optimizer
    )
  } else {
    topo <- gamma0
  }
  .saint_int_list_output(topo$model, parsed, gamma0$scores, topo$scores)
}
