.saint_spc_read_tables <- function(si) {
  tables <- .saint_int_prepare_tables(si)
  inter <- tables$inter
  prey <- tables$prey
  bait <- .saint_int_index_bait_replicates(tables$bait)

  inter <- .saint_spc_index_interactions(inter, prey, bait)
  ubait <- sort(unique(bait$baitId[!bait$is_ctrl]))
  ip_idx_to_bait_no <- match(bait$baitId[!bait$is_ctrl], ubait)

  matrices <- .saint_spc_build_observation_matrices(inter, prey, bait)
  ui <- .saint_int_unique_test_interactions(inter, prey, ip_idx_to_bait_no)
  matrices$test_mat <- .saint_spc_mask_self_interactions(matrices$test_mat, ui)

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

.saint_spc_index_interactions <- function(inter, prey, bait) {
  inter$quant <- as.integer(as.numeric(inter$quant))
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

.saint_spc_build_observation_matrices <- function(inter, prey, bait) {
  nprey <- nrow(prey)
  n_test_seen <- max(c(0L, bait$col[!bait$is_ctrl]))
  n_ctrl_seen <- max(c(0L, bait$col[bait$is_ctrl]))
  test_mat <- matrix(0L, nrow = nprey, ncol = n_test_seen)
  ctrl_mat <- matrix(0L, nrow = nprey, ncol = n_ctrl_seen)
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

.saint_spc_mask_self_interactions <- function(test_mat, ui) {
  for (entry in ui) {
    if (identical(entry$preyGeneId, entry$baitId)) {
      test_mat[entry$row, entry$colIds] <- 0L
    }
  }
  test_mat
}

.saint_spc_build_model <- function(
  parsed,
  p2p_mapping = NULL,
  L = 100L,
  frequency = 0.5
) {
  if (L <= 0) {
    stop("L is must be a positive integer.")
  }

  dimensions <- .saint_int_model_dimensions(parsed)
  replicates <- .saint_int_replicate_layout(parsed, dimensions$nbait)
  controls <- .saint_spc_estimate_control_model(parsed$ctrl_mat, L = L)
  tests <- .saint_spc_estimate_test_model(
    parsed$test_mat,
    replicates$bait_no_to_ip_idxes,
    controls
  )
  p2p_mapping <- .saint_int_default_p2p_mapping(p2p_mapping, dimensions$nprey)

  .saint_spc_new_model(
    parsed,
    dimensions,
    replicates,
    controls,
    tests,
    p2p_mapping,
    frequency
  )
}

.saint_spc_estimate_control_model <- function(ctrl_mat, L) {
  L <- min(ncol(ctrl_mat), L)
  eta <- numeric(nrow(ctrl_mat))
  lambda2_false <- rep(0.1, nrow(ctrl_mat))
  ctrl_mean <- as.integer(rowMeans(ctrl_mat))

  for (prey_row in seq_len(nrow(ctrl_mat))) {
    ctrl <- sort(ctrl_mat[prey_row, ], decreasing = TRUE)
    eta[[prey_row]] <- max(mean(ctrl[seq_len(L)]), 0.1)
    variance <- .saint_spc_var_sample(ctrl[seq_len(L)])
    if (variance > eta[[prey_row]]) {
      lambda2_false[[prey_row]] <- 1 - sqrt(eta[[prey_row]] / variance)
    }
  }

  list(
    eta = eta,
    ctrl_mean = ctrl_mean,
    lambda2_false = lambda2_false,
    lambda2_true = rep(0.1, nrow(ctrl_mat))
  )
}

.saint_spc_estimate_test_model <- function(
  test_mat_DATA,
  bait_no_to_ip_idxes,
  controls
) {
  test_mat <- .saint_int_test_replicate_matrix(
    test_mat_DATA,
    bait_no_to_ip_idxes
  )
  list(
    test_mat = test_mat,
    test_mat1 = test_mat,
    test_SS = .saint_spc_test_sum_matrix(test_mat),
    d = .saint_spc_signal_shift(test_mat_DATA, controls$eta),
    Z = .saint_spc_initial_z(test_mat, controls$eta)
  )
}

.saint_spc_new_model <- function(
  parsed,
  dimensions,
  replicates,
  controls,
  tests,
  p2p_mapping,
  frequency
) {
  list(
    nprey = dimensions$nprey,
    nbait = dimensions$nbait,
    test_mat_DATA = parsed$test_mat,
    test_mat = tests$test_mat,
    test_mat1 = tests$test_mat1,
    test_SS = tests$test_SS,
    n_rep_vec = replicates$n_rep_vec,
    p2p_mapping = p2p_mapping,
    ctrl_mean = controls$ctrl_mean,
    n_ctrl_ip = dimensions$n_ctrl_ip,
    apply_MRF = .saint_spc_apply_mrf(tests$Z, parsed$test_mat, frequency),
    Z = tests$Z,
    beta0 = 0,
    beta1 = 0,
    gamma = 0,
    eta = controls$eta,
    d = tests$d,
    lambda2_true = controls$lambda2_true,
    lambda2_false = controls$lambda2_false
  )
}

.saint_spc_var_sample <- function(x) {
  if (length(x) <= 1) {
    return(0)
  }
  stats::var(x)
}

.saint_spc_test_sum_matrix <- function(test_mat) {
  out <- matrix(0, nrow = nrow(test_mat), ncol = ncol(test_mat))
  for (prey_row in seq_len(nrow(test_mat))) {
    for (bait_col in seq_len(ncol(test_mat))) {
      out[prey_row, bait_col] <- sum(test_mat[[prey_row, bait_col]])
    }
  }
  out
}

.saint_spc_signal_shift <- function(test_mat_DATA, eta) {
  d <- numeric(nrow(test_mat_DATA))
  for (prey_row in seq_len(nrow(test_mat_DATA))) {
    d[[prey_row]] <- 4 * eta[[prey_row]]
  }
  d
}

.saint_spc_initial_z <- function(test_mat, eta) {
  .saint_int_list_matrix(
    nrow(test_mat),
    ncol(test_mat),
    function(prey_row, bait_col) {
      test_mat[[prey_row, bait_col]] > (2 * eta[[prey_row]])
    }
  )
}

.saint_spc_apply_mrf <- function(Z, test_mat_DATA, frequency) {
  apply_mrf <- logical(nrow(Z))
  for (prey_row in seq_len(nrow(Z))) {
    z_sum <- 0
    for (bait_col in seq_len(ncol(Z))) {
      z_sum <- z_sum + sum(Z[[prey_row, bait_col]])
    }
    apply_mrf[[prey_row]] <- z_sum < ncol(test_mat_DATA) * frequency
  }
  apply_mrf
}

.saint_spc_gp_log_pmf <- function(x, lambda1, lambda2) {
  if (x == 0) {
    return(-lambda1)
  }
  tmp <- lambda1 + x * lambda2
  log(lambda1) + (x - 1) * log(tmp) - tmp - lgamma(x + 1)
}

.saint_spc_gp_log_pmf_mean <- function(k, mean, lambda2) {
  .saint_spc_gp_log_pmf(as.integer(k), mean * (1 - lambda2), lambda2)
}

.saint_spc_llik_mrf_gamma0 <- function(model, beta1) {
  loglik <- 0
  mrf_true <- exp(beta1)
  mrf_false <- exp(0)
  for (prey_row in seq_len(model$nprey)) {
    for (bait_col in seq_len(model$nbait)) {
      y <- model$test_mat1[[prey_row, bait_col]]
      if (all(y == 0)) {
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

.saint_spc_llik_mrf <- function(model, beta1, gamma, gsum_mat) {
  loglik <- 0
  mrf_false <- exp(0)
  for (prey_row in seq_len(model$nprey)) {
    for (bait_col in seq_len(model$nbait)) {
      y <- model$test_mat1[[prey_row, bait_col]]
      if (all(y == 0)) {
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

.saint_spc_loglikelihood <- function(model) {
  loglik <- 0
  for (prey_row in seq_len(model$nprey)) {
    for (bait_col in seq_len(model$nbait)) {
      y <- model$test_mat1[[prey_row, bait_col]]
      if (all(y == 0)) {
        next
      }
      k <- model$Z[[prey_row, bait_col]]
      m <- model$n_rep_vec[[bait_col]]
      log_mrf <- .saint_int_cell_log_mrf(model, prey_row, bait_col)
      prod_val <- 1
      for (rep in seq_len(m)) {
        true_log <- .saint_spc_gp_log_pmf_mean(
          min(y[[rep]], model$eta[[prey_row]] + model$d[[prey_row]]),
          model$eta[[prey_row]] + model$d[[prey_row]],
          model$lambda2_true[[prey_row]]
        )
        false_log <- .saint_spc_gp_log_pmf_mean(
          max(y[[rep]], model$eta[[prey_row]]),
          model$eta[[prey_row]],
          model$lambda2_false[[prey_row]]
        )
        prod_val <- prod_val *
          (if (k[[rep]]) {
            exp(log_mrf[["true"]] + true_log)
          } else {
            exp(log_mrf[["false"]] + false_log)
          }) /
          (exp(log_mrf[["true"]]) + exp(log_mrf[["false"]]))
      }
      loglik <- loglik + log(prod_val) / m
    }
  }
  loglik
}

.saint_spc_precalculate_rep_logpdf <- function(model) {
  .saint_int_list_matrix(
    model$nprey,
    model$nbait,
    function(prey_row, bait_col) {
      y <- model$test_mat1[[prey_row, bait_col]]
      m <- model$n_rep_vec[[bait_col]]
      vals <- matrix(0, nrow = m, ncol = 2)
      for (rep in seq_len(m)) {
        vals[rep, 1] <- .saint_spc_gp_log_pmf_mean(
          min(y[[rep]], model$eta[[prey_row]] + model$d[[prey_row]]),
          model$eta[[prey_row]] + model$d[[prey_row]],
          model$lambda2_true[[prey_row]]
        )
        vals[rep, 2] <- .saint_spc_gp_log_pmf_mean(
          max(y[[rep]], model$eta[[prey_row]]),
          model$eta[[prey_row]],
          model$lambda2_false[[prey_row]]
        )
      }
      vals
    }
  )
}

.saint_spc_optimize_gamma0 <- function(model, optimizer) {
  old_beta1 <- model$beta1
  oldf <- .saint_spc_llik_mrf_gamma0(model, old_beta1)
  if (
    identical(optimizer, "nloptr") && requireNamespace("nloptr", quietly = TRUE)
  ) {
    res <- nloptr::nloptr(
      x0 = old_beta1,
      eval_f = function(x) -.saint_spc_llik_mrf_gamma0(model, x[[1]]),
      lb = -15,
      ub = 15,
      opts = list(algorithm = "NLOPT_LN_COBYLA", ftol_abs = 1e-4, maxeval = 100)
    )
    new_beta1 <- res$solution[[1]]
    maxf <- -res$objective
  } else {
    res <- stats::optim(
      par = old_beta1,
      fn = function(x) -.saint_spc_llik_mrf_gamma0(model, x[[1]]),
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

.saint_spc_optimize_mrf <- function(model, optimizer) {
  gsum_mat <- .saint_int_gsum_mat(model)
  old_beta1 <- model$beta1
  old_gamma <- model$gamma
  oldf <- .saint_spc_llik_mrf(model, old_beta1, old_gamma, gsum_mat)
  if (
    identical(optimizer, "nloptr") && requireNamespace("nloptr", quietly = TRUE)
  ) {
    res <- nloptr::nloptr(
      x0 = c(old_beta1, old_gamma),
      eval_f = function(x) {
        -.saint_spc_llik_mrf(model, x[[1]], x[[2]], gsum_mat)
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
      fn = function(x) -.saint_spc_llik_mrf(model, x[[1]], x[[2]], gsum_mat),
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
  if (maxf <= .saint_spc_llik_mrf(model, x[[1]], 0, gsum_mat)) {
    model$gamma <- 0
  }
  if (maxf < oldf) {
    model$beta1 <- old_beta1
    model$gamma <- old_gamma
  }
  model
}

.saint_spc_icm_z <- function(model) {
  pre_calc <- .saint_spc_precalculate_rep_logpdf(model)

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

.saint_spc_icms <- function(model, with_gamma, optimizer) {
  newllik <- .saint_spc_loglikelihood(model)
  model$gamma <- 0
  for (iter in seq_len(15)) {
    oldllik <- newllik
    model <- .saint_spc_icm_z(model)
    model <- if (with_gamma) {
      .saint_spc_optimize_mrf(model, optimizer)
    } else {
      .saint_spc_optimize_gamma0(model, optimizer)
    }
    newllik <- .saint_spc_loglikelihood(model)
    if (newllik >= oldllik && exp(newllik - oldllik) - 1 < 1e-3) {
      break
    }
  }
  model
}

.saint_spc_calculate_score <- function(model, R = 100L) {
  average_score <- matrix(0, nrow = model$nprey, ncol = model$nbait)
  max_score <- matrix(0, nrow = model$nprey, ncol = model$nbait)
  min_log_odds_score <- matrix(0, nrow = model$nprey, ncol = model$nbait)

  for (prey_row in seq_len(model$nprey)) {
    for (bait_col in seq_len(model$nbait)) {
      cell_score <- .saint_spc_calculate_cell_score(
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

.saint_spc_calculate_cell_score <- function(model, prey_row, bait_col, R) {
  y <- model$test_mat1[[prey_row, bait_col]]
  m <- model$n_rep_vec[[bait_col]]
  log_mrf <- .saint_int_cell_log_mrf(model, prey_row, bait_col)
  rep_scores <- .saint_spc_rep_scores(model, prey_row, y, log_mrf)
  tmp_scores <- sort(rep_scores$score, decreasing = TRUE)
  tmp_odds <- sort(rep_scores$odds, decreasing = TRUE)
  max_rep <- R

  average <- if (m > max_rep) {
    mean(tmp_scores[seq_len(max_rep)])
  } else {
    mean(tmp_scores)
  }
  maximum <- tmp_scores[[1]]
  if (max(y) == 1) {
    average <- 0
    maximum <- 0
  }

  list(
    average = average,
    maximum = maximum,
    odds = tmp_odds[[length(tmp_odds)]]
  )
}

.saint_spc_rep_scores <- function(model, prey_row, y, log_mrf) {
  score <- numeric(length(y))
  odds <- numeric(length(y))

  for (rep in seq_along(y)) {
    tmp_mean <- if (y[[rep]] >= 5) {
      min(
        model$eta[[prey_row]] + model$d[[prey_row]],
        5 * model$eta[[prey_row]]
      )
    } else {
      model$eta[[prey_row]] + model$d[[prey_row]]
    }
    true_log <- log_mrf[["true"]]
    false_log <- log_mrf[["false"]] +
      .saint_spc_gp_log_pmf_mean(
        max(y[[rep]], model$eta[[prey_row]]),
        model$eta[[prey_row]],
        model$lambda2_false[[prey_row]]
      ) -
      .saint_spc_gp_log_pmf_mean(
        min(y[[rep]], tmp_mean),
        tmp_mean,
        model$lambda2_true[[prey_row]]
      )
    score[[rep]] <- if (
      y[[rep]] <= 1 || y[[rep]] <= model$ctrl_mean[[prey_row]]
    ) {
      0
    } else {
      .saint_int_logit_probability(true_log, false_log)
    }
    odds[[rep]] <- true_log - false_log
  }

  list(score = score, odds = odds)
}

.saint_spc_computation <- function(model, with_gamma, optimizer) {
  fitted <- .saint_spc_icms(
    model,
    with_gamma = with_gamma,
    optimizer = optimizer
  )
  list(
    model = fitted,
    scores = .saint_spc_calculate_score(
      fitted,
      R = attr(model, "R", exact = TRUE) %||% 100L
    )
  )
}

.saint_spc_list_output <- function(model, parsed, scores, topo_scores) {
  out <- lapply(parsed$ui, function(ui) {
    .saint_spc_output_row(ui, model, parsed, scores, topo_scores)
  })
  if (!length(out)) {
    return(data.frame())
  }
  res <- do.call(rbind, out)
  numeric_cols <- c(
    "SpecSum",
    "AvgSpec",
    "AvgP",
    "MaxP",
    "TopoAvgP",
    "TopoMaxP",
    "SaintScore",
    "logOddsScore",
    "FoldChange",
    "BFDR"
  )
  for (col in numeric_cols) {
    res[[col]] <- round(res[[col]], 2)
  }
  res
}

.saint_spc_output_row <- function(ui, model, parsed, scores, topo_scores) {
  prey_row <- ui$row
  bait_col <- ui$baitcol
  counts <- model$test_mat[[prey_row, bait_col]]
  avg_spec <- mean(counts)
  ctrl_counts <- parsed$ctrl_mat[prey_row, ]
  avg_p <- scores$average[prey_row, bait_col]
  topo_avg_p <- topo_scores$average[prey_row, bait_col]

  data.frame(
    Bait = ui$baitId,
    Prey = ui$preyId,
    PreyGene = ui$preyGeneId,
    Spec = paste(counts, collapse = "|"),
    SpecSum = sum(counts),
    AvgSpec = avg_spec,
    NumReplicates = length(counts),
    ctrlCounts = paste(ctrl_counts, collapse = "|"),
    AvgP = avg_p,
    MaxP = scores$maximum[prey_row, bait_col],
    TopoAvgP = topo_avg_p,
    TopoMaxP = topo_scores$maximum[prey_row, bait_col],
    SaintScore = max(topo_avg_p, avg_p),
    logOddsScore = topo_scores$odds[prey_row, bait_col],
    FoldChange = avg_spec / model$eta[[prey_row]],
    BFDR = .saint_int_bfdr(avg_p, scores$average),
    boosted_by = .saint_int_boosted_by(model, parsed, prey_row, bait_col),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

saint_spc_r <- function(
  si,
  optimizer = c("base", "nloptr"),
  R = 100L,
  L = 100L,
  frequency = 0.5,
  p2p_mapping = NULL
) {
  optimizer <- match.arg(optimizer)
  parsed <- .saint_spc_read_tables(si)
  model <- .saint_spc_build_model(
    parsed,
    p2p_mapping = p2p_mapping,
    L = L,
    frequency = frequency
  )
  attr(model, "R") <- R

  gamma0 <- .saint_spc_computation(
    model,
    with_gamma = FALSE,
    optimizer = optimizer
  )
  if (!is.null(p2p_mapping)) {
    topo <- .saint_spc_computation(
      model,
      with_gamma = TRUE,
      optimizer = optimizer
    )
  } else {
    topo <- gamma0
  }
  .saint_spc_list_output(topo$model, parsed, gamma0$scores, topo$scores)
}
