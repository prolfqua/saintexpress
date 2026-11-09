# Intensity model (SAINTexpress-int): standardised log intensities follow
# normal distributions; see saint_engine.R for the meaning of the ops entries.

.saint_int_ops <- function() {
    list(
        transform = function(quant) log(as.numeric(quant)),
        missing = NA_real_,
        prepare = .saint_int_normalize,
        parameters = .saint_int_parameters,
        empty = function(y) all(is.na(y)),
        log_density = .saint_int_log_density,
        weight = function(log_mrf, log_density) exp(log_mrf) * exp(log_density),
        rep_scores = .saint_int_rep_scores,
        values = .saint_int_values,
        names = c(
            quant = "Intensity",
            sum = "IntensitySum",
            avg = "AvgIntensity",
            ctrl = "ctrlIntensity",
            odds = "OddsScore"
        ),
        digits = 3
    )
}

# Standardises all log intensities with their pooled mean and SD.
.saint_int_normalize <- function(parsed) {
    vals <- c(as.vector(parsed$test_mat), as.vector(parsed$ctrl_mat))
    vals <- vals[is.finite(vals) & vals > -1e5]
    if (length(vals) < 2) {
        stop(
            "not enough finite intensity values for SAINTexpress-int normalization"
        )
    }
    center <- mean(vals)
    scale <- sqrt(var(vals) * (length(vals) - 1) / length(vals))
    parsed$test_mat <- (parsed$test_mat - center) / scale
    parsed$ctrl_mat <- (parsed$ctrl_mat - center) / scale
    parsed
}

# Background from the controls: mean eta and SD sd_false per prey, with the
# value t standing in for missing intensities. True interactions have mean
# eta + d and SD sd_true.
.saint_int_parameters <- function(parsed, test_mat) {
    ctrl_mat <- parsed$ctrl_mat
    t_value <- .saint_int_missing_threshold(ctrl_mat)
    ctrl_lse <- .saint_int_control_lse(ctrl_mat, parsed$prey$preyId)
    priors <- .saint_int_control_priors(ctrl_mat, t_value, ctrl_lse$sd_NA, ctrl_lse$lse)
    Z <- .saint_list_matrix(nrow(test_mat), ncol(test_mat), function(prey_row, bait_col) {
        y <- test_mat[[prey_row, bait_col]]
        ifelse(is.na(y), FALSE, y > (priors$eta[[prey_row]] + log(5)))
    })
    list(
        t = t_value,
        eta = priors$eta,
        d = .saint_int_signal_shift(parsed$test_mat, priors$eta),
        sd_false = priors$sd_ctrl,
        sd_true = .saint_int_true_sd(test_mat, Z, priors$sd_ctrl, ctrl_lse$sd_NA),
        Z = Z
    )
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

.saint_int_var_mle <- function(x) {
    if (length(x) <= 1) {
        return(0)
    }
    mean((x - mean(x))^2)
}

# Regression of log control SD on control mean over complete control profiles,
# used to predict the SD of preys seen in fewer than 3 controls.
.saint_int_control_lse <- function(ctrl_mat, prey_ids) {
    complete <- complete.cases(ctrl_mat)
    complete_ctrl <- ctrl_mat[complete, , drop = FALSE]
    log_sd_ctrls <- apply(complete_ctrl, 1, function(ctrl) {
        log(sqrt(.saint_int_var_mle(ctrl)))
    })
    mean_ctrls <- rowMeans(complete_ctrl)
    if (!length(log_sd_ctrls)) {
        return(list(lse = c(0, 0), sd_NA = 1))
    }

    finite_sd <- is.finite(log_sd_ctrls)
    zero_variance_ids <- prey_ids[complete][!finite_sd]
    if (!any(finite_sd)) {
        stop(
            "SAINTexpress-int cannot estimate control variability: all ",
            length(log_sd_ctrls),
            " complete control profiles have zero variance (",
            paste(zero_variance_ids, collapse = ", "),
            ").",
            call. = FALSE
        )
    }

    sorted_log_sd_ctrls <- sort(log_sd_ctrls[finite_sd])
    median_sd <- exp(sorted_log_sd_ctrls[[length(sorted_log_sd_ctrls) %/% 2 + 1L]])
    if (length(zero_variance_ids)) {
        warning(
            "SAINTexpress-int: ",
            length(zero_variance_ids),
            " complete control ",
            if (length(zero_variance_ids) == 1L) "profile has" else "profiles have",
            " zero variance (",
            paste(zero_variance_ids, collapse = ", "),
            "). Using the median control SD, matching native SAINTexpress.",
            call. = FALSE
        )
        return(list(lse = c(log(median_sd), 0), sd_NA = median_sd))
    }
    if (length(log_sd_ctrls) < 2 || var(mean_ctrls) == 0) {
        return(list(
            lse = c(log(median(exp(log_sd_ctrls), na.rm = TRUE)), 0),
            sd_NA = median_sd
        ))
    }
    fit <- lm(log_sd_ctrls ~ mean_ctrls)
    list(lse = unname(coef(fit)), sd_NA = median_sd)
}

.saint_int_control_priors <- function(ctrl_mat, t_value, sd_NA, lse) {
    eta <- numeric(nrow(ctrl_mat))
    sd_ctrl <- numeric(nrow(ctrl_mat))
    for (prey_row in seq_len(nrow(ctrl_mat))) {
        detected <- ctrl_mat[prey_row, ]
        detected <- detected[!is.na(detected)]
        n_detected <- length(detected)
        if (n_detected == 0) {
            eta[[prey_row]] <- t_value
            sd <- sd_NA
        } else {
            eta[[prey_row]] <- mean(detected)
            sd <- if (n_detected >= 3) {
                sqrt(.saint_int_var_mle(detected))
            } else {
                exp(eta[[prey_row]] * lse[[2]] + lse[[1]])
            }
        }
        sd_ctrl[[prey_row]] <- max(sd_NA, sd)
    }
    list(eta = eta, sd_ctrl = sd_ctrl)
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

.saint_int_true_sd <- function(test_mat, Z, sd_false, sd_NA) {
    sd_true <- numeric(nrow(test_mat))
    for (prey_row in seq_len(nrow(test_mat))) {
        z1 <- unlist(lapply(seq_len(ncol(test_mat)), function(bait_col) {
            test_mat[[prey_row, bait_col]][Z[[prey_row, bait_col]]]
        }))
        sd <- if (length(z1) >= 2) sqrt(var(z1)) else sd_false[[prey_row]]
        sd_true[[prey_row]] <- max(sd_NA, sd)
    }
    sd_true
}

# Normal log-density, with x pulled back to within 5 SD of the mean.
.saint_int_log_pdf <- function(x, mu, sigma) {
    if (abs(x - mu) / sigma > 5) {
        x <- (if (x < 0) -1 else 1) * 5 * sigma + mu
    }
    dnorm(x, mean = mu, sd = sigma, log = TRUE)
}

.saint_int_log_density <- function(model, prey_row, y) {
    eta <- model$eta[[prey_row]]
    mu <- eta + model$d[[prey_row]]
    y[is.na(y)] <- model$t
    vals <- matrix(0, nrow = length(y), ncol = 2)
    for (rep in seq_along(y)) {
        vals[rep, 1] <- .saint_int_log_pdf(min(y[[rep]], mu), mu, model$sd_true[[prey_row]])
        vals[rep, 2] <- .saint_int_log_pdf(max(y[[rep]], eta), eta, model$sd_false[[prey_row]])
    }
    vals
}

.saint_int_rep_scores <- function(model, prey_row, y, log_mrf) {
    log_density <- .saint_int_log_density(model, prey_row, y)
    true_log <- log_mrf[["true"]] + log_density[, 1]
    false_log <- log_mrf[["false"]] + log_density[, 2]
    score <- mapply(.saint_logit_probability, true_log, false_log)
    score[is.na(y)] <- 0
    list(score = score, odds = true_log - false_log)
}

.saint_int_values <- function(model, parsed, prey_row, bait_col) {
    quant_orig <- exp(model$test_mat[[prey_row, bait_col]])
    quant <- quant_orig
    quant[is.na(quant)] <- exp(model$t)
    format_num <- function(x) ifelse(is.na(x), ".", sprintf("%.3f", x))
    list(
        quant = paste(format_num(quant_orig), collapse = "|"),
        sum = sum(quant),
        avg = mean(quant),
        n = length(quant),
        ctrl = paste(format_num(exp(parsed$ctrl_mat[prey_row, ])), collapse = "|"),
        fold_change = mean(quant) / exp(model$eta[[prey_row]])
    )
}
