# Spectral-count model (SAINTexpress-spc): counts follow a generalized Poisson
# distribution; see saint_engine.R for the meaning of the ops entries.

.saint_spc_ops <- function(L) {
    list(
        transform = function(quant) as.integer(as.numeric(quant)),
        missing = 0L,
        prepare = identity,
        parameters = function(parsed, test_mat) {
            .saint_spc_parameters(parsed$ctrl_mat, test_mat, L)
        },
        empty = function(y) all(y == 0),
        log_density = .saint_spc_log_density,
        weight = function(log_mrf, log_density) exp(log_mrf + log_density),
        rep_scores = .saint_spc_rep_scores,
        values = .saint_spc_values,
        names = c(
            quant = "Spec",
            sum = "SpecSum",
            avg = "AvgSpec",
            ctrl = "ctrlCounts",
            odds = "logOddsScore"
        ),
        digits = 2
    )
}

# Background (false interaction) from the L highest control counts of each
# prey: mean eta and dispersion lambda2_false. True interactions have mean
# eta + d with d = 4 eta.
.saint_spc_parameters <- function(ctrl_mat, test_mat, L) {
    L <- min(ncol(ctrl_mat), L)
    nprey <- nrow(ctrl_mat)
    eta <- numeric(nprey)
    lambda2_false <- rep(0.1, nprey)
    for (prey_row in seq_len(nprey)) {
        ctrl <- sort(ctrl_mat[prey_row, ], decreasing = TRUE)[seq_len(L)]
        eta[[prey_row]] <- max(mean(ctrl), 0.1)
        variance <- if (L > 1) var(ctrl) else 0
        if (variance > eta[[prey_row]]) {
            lambda2_false[[prey_row]] <- 1 - sqrt(eta[[prey_row]] / variance)
        }
    }
    list(
        eta = eta,
        d = 4 * eta,
        ctrl_mean = as.integer(rowMeans(ctrl_mat)),
        lambda2_true = rep(0.1, nprey),
        lambda2_false = lambda2_false,
        Z = .saint_list_matrix(nprey, ncol(test_mat), function(prey_row, bait_col) {
            test_mat[[prey_row, bait_col]] > (2 * eta[[prey_row]])
        })
    )
}

# Generalized Poisson log-probability of count k with the given mean and
# dispersion lambda2.
.saint_spc_log_pmf <- function(k, mean, lambda2) {
    k <- as.integer(k)
    lambda1 <- mean * (1 - lambda2)
    if (k == 0) {
        return(-lambda1)
    }
    tmp <- lambda1 + k * lambda2
    log(lambda1) + (k - 1) * log(tmp) - tmp - lgamma(k + 1)
}

.saint_spc_log_density <- function(model, prey_row, y) {
    eta <- model$eta[[prey_row]]
    mu <- eta + model$d[[prey_row]]
    vals <- matrix(0, nrow = length(y), ncol = 2)
    for (rep in seq_along(y)) {
        vals[rep, 1] <- .saint_spc_log_pmf(min(y[[rep]], mu), mu, model$lambda2_true[[prey_row]])
        vals[rep, 2] <- .saint_spc_log_pmf(max(y[[rep]], eta), eta, model$lambda2_false[[prey_row]])
    }
    vals
}

.saint_spc_rep_scores <- function(model, prey_row, y, log_mrf) {
    eta <- model$eta[[prey_row]]
    mu <- eta + model$d[[prey_row]]
    score <- numeric(length(y))
    odds <- numeric(length(y))
    for (rep in seq_along(y)) {
        tmp_mean <- if (y[[rep]] >= 5) min(mu, 5 * eta) else mu
        true_log <- log_mrf[["true"]]
        false_log <- log_mrf[["false"]] +
            .saint_spc_log_pmf(max(y[[rep]], eta), eta, model$lambda2_false[[prey_row]]) -
            .saint_spc_log_pmf(min(y[[rep]], tmp_mean), tmp_mean, model$lambda2_true[[prey_row]])
        score[[rep]] <- if (y[[rep]] <= 1 || y[[rep]] <= model$ctrl_mean[[prey_row]]) {
            0
        } else {
            .saint_logit_probability(true_log, false_log)
        }
        odds[[rep]] <- true_log - false_log
    }
    list(score = score, odds = odds)
}

.saint_spc_values <- function(model, parsed, prey_row, bait_col) {
    counts <- model$test_mat[[prey_row, bait_col]]
    list(
        quant = paste(counts, collapse = "|"),
        sum = sum(counts),
        avg = mean(counts),
        n = length(counts),
        ctrl = paste(parsed$ctrl_mat[prey_row, ], collapse = "|"),
        fold_change = mean(counts) / model$eta[[prey_row]]
    )
}
