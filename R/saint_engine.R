# Engine shared by the spectral-count and intensity models. The mode-specific
# parts come in `ops`, a list built by .saint_spc_ops() or .saint_int_ops():
#   transform(quant)       numeric scale of inter$quant used by the model
#   missing                value of an unobserved prey in an IP
#   prepare(parsed)        transformation of the observation matrices
#   parameters(parsed, test_mat)  background and signal parameters, initial Z
#   empty(y)               TRUE if a prey was not observed with a bait
#   log_density(model, prey_row, y)  replicates x (true, false) log densities
#   weight(log_mrf, log_density)     likelihood term of one replicate
#   rep_scores(model, prey_row, y, log_mrf)  replicate probabilities and odds
#   values(model, parsed, prey_row, bait_col)  mode-specific output values
#   names, digits          output column names and rounding

.saint_run <- function(si, ops, optimizer, R, p2p_mapping) {
    parsed <- ops$prepare(.saint_read_tables(si, ops))
    model <- .saint_build_model(parsed, p2p_mapping, ops)
    gamma0 <- .saint_computation(model, FALSE, optimizer, R, ops)
    topo <- if (is.null(p2p_mapping)) {
        gamma0
    } else {
        .saint_computation(model, TRUE, optimizer, R, ops)
    }
    .saint_list_output(topo$model, parsed, gamma0$scores, topo$scores, ops)
}

.saint_read_tables <- function(si, ops) {
    tables <- .saint_resolve_tables(si)
    prey <- tables$prey
    prey$row <- seq_len(nrow(prey))
    bait <- .saint_index_bait_replicates(tables$bait)

    inter <- tables$inter
    inter$quant <- ops$transform(inter$quant)
    inter$row <- prey$row[match(inter$preyId, prey$preyId)]
    bait_match <- match(inter$ipId, bait$ipId)
    inter$col <- bait$col[bait_match]
    inter$is_ctrl <- bait$is_ctrl[bait_match]

    ubait <- sort(unique(bait$baitId[!bait$is_ctrl]))
    ip_idx_to_bait_no <- match(bait$baitId[!bait$is_ctrl], ubait)
    ui <- .saint_unique_test_interactions(inter, prey, ip_idx_to_bait_no)

    test_mat <- matrix(ops$missing, nrow = nrow(prey), ncol = sum(!bait$is_ctrl))
    ctrl_mat <- matrix(ops$missing, nrow = nrow(prey), ncol = sum(bait$is_ctrl))
    for (i in seq_len(nrow(inter))) {
        if (inter$is_ctrl[[i]]) {
            ctrl_mat[inter$row[[i]], inter$col[[i]]] <- inter$quant[[i]]
        } else {
            test_mat[inter$row[[i]], inter$col[[i]]] <- inter$quant[[i]]
        }
    }
    for (entry in ui) {
        if (identical(entry$preyGeneId, entry$baitId)) {
            test_mat[entry$row, entry$colIds] <- ops$missing
        }
    }

    list(
        prey = prey,
        ubait = ubait,
        ip_idx_to_bait_no = ip_idx_to_bait_no,
        test_mat = test_mat,
        ctrl_mat = ctrl_mat,
        ui = ui
    )
}

# Numbers test and control IPs separately, in table order.
.saint_index_bait_replicates <- function(bait) {
    bait$is_ctrl <- bait$CorT != "T"
    bait$col <- integer(nrow(bait))
    bait$col[!bait$is_ctrl] <- seq_len(sum(!bait$is_ctrl))
    bait$col[bait$is_ctrl] <- seq_len(sum(bait$is_ctrl))
    bait
}

# One entry per observed (prey, bait) pair, in order of first appearance.
.saint_unique_test_interactions <- function(inter, prey, ip_idx_to_bait_no) {
    ui_key <- new.env(parent = emptyenv())
    ui <- list()
    test_inter <- inter[!inter$is_ctrl, , drop = FALSE]
    for (i in seq_len(nrow(test_inter))) {
        prey_row <- test_inter$row[[i]]
        bait_col <- ip_idx_to_bait_no[[test_inter$col[[i]]]]
        key <- paste(prey_row, bait_col, sep = "\r")
        idx <- ui_key[[key]]
        if (is.null(idx)) {
            idx <- length(ui) + 1L
            ui_key[[key]] <- idx
            ui[[idx]] <- list(
                baitId = test_inter$baitId[[i]],
                preyId = test_inter$preyId[[i]],
                preyGeneId = prey$preyGeneId[[prey_row]],
                row = prey_row,
                baitcol = bait_col,
                colIds = integer()
            )
        }
        ui[[idx]]$colIds <- c(ui[[idx]]$colIds, test_inter$col[[i]])
    }
    ui
}

.saint_build_model <- function(parsed, p2p_mapping, ops) {
    nprey <- nrow(parsed$prey)
    nbait <- length(parsed$ubait)
    replicates <- lapply(seq_len(nbait), function(bait_col) {
        which(parsed$ip_idx_to_bait_no == bait_col)
    })
    test_mat <- .saint_list_matrix(nprey, nbait, function(prey_row, bait_col) {
        parsed$test_mat[prey_row, replicates[[bait_col]]]
    })
    if (is.null(p2p_mapping)) {
        p2p_mapping <- rep(list(integer()), nprey)
    }
    c(
        list(
            nprey = nprey,
            nbait = nbait,
            test_mat = test_mat,
            n_rep_vec = lengths(replicates),
            p2p_mapping = p2p_mapping,
            beta0 = 0,
            beta1 = 0,
            gamma = 0
        ),
        ops$parameters(parsed, test_mat)
    )
}

.saint_list_matrix <- function(nrow, ncol, fn) {
    out <- vector("list", nrow * ncol)
    dim(out) <- c(nrow, ncol)
    for (row_idx in seq_len(nrow)) {
        for (col_idx in seq_len(ncol)) {
            out[[row_idx, col_idx]] <- fn(row_idx, col_idx)
        }
    }
    out
}

.saint_computation <- function(model, with_gamma, optimizer, R, ops) {
    fitted <- .saint_icms(model, with_gamma, optimizer, ops)
    list(model = fitted, scores = .saint_calculate_score(fitted, R, ops))
}

# Iterated conditional modes: alternates updating the indicators Z and the
# MRF parameters until the likelihood stops improving.
.saint_icms <- function(model, with_gamma, optimizer, ops) {
    pre_calc <- .saint_list_matrix(model$nprey, model$nbait, function(prey_row, bait_col) {
        ops$log_density(model, prey_row, model$test_mat[[prey_row, bait_col]])
    })
    newllik <- .saint_loglikelihood(model, pre_calc, ops)
    model$gamma <- 0
    for (iter in seq_len(15)) {
        oldllik <- newllik
        model <- .saint_icm_z(model, pre_calc)
        model <- .saint_optimize_mrf(model, with_gamma, optimizer, ops$empty)
        newllik <- .saint_loglikelihood(model, pre_calc, ops)
        if (newllik >= oldllik && exp(newllik - oldllik) - 1 < 1e-3) {
            break
        }
    }
    model
}

# Sum over the preys linked to `prey_row` of their mean Z with the bait.
.saint_gsum <- function(model, prey_row, bait_col) {
    links <- model$p2p_mapping[[prey_row]]
    if (!length(links)) {
        return(0)
    }
    sum(vapply(
        links,
        function(linked_prey) mean(model$Z[[linked_prey, bait_col]]),
        numeric(1)
    ))
}

.saint_cell_log_mrf <- function(model, prey_row, bait_col) {
    gsum <- if (model$gamma == 0) 0 else .saint_gsum(model, prey_row, bait_col)
    c(true = model$beta1 + model$gamma * gsum, false = model$beta0)
}

.saint_loglikelihood <- function(model, pre_calc, ops) {
    loglik <- 0
    for (prey_row in seq_len(model$nprey)) {
        for (bait_col in seq_len(model$nbait)) {
            if (ops$empty(model$test_mat[[prey_row, bait_col]])) {
                next
            }
            k <- model$Z[[prey_row, bait_col]]
            log_mrf <- .saint_cell_log_mrf(model, prey_row, bait_col)
            log_density <- pre_calc[[prey_row, bait_col]]
            norm <- exp(log_mrf[["true"]]) + exp(log_mrf[["false"]])
            prod_val <- 1
            for (rep in seq_along(k)) {
                term <- if (k[[rep]]) {
                    ops$weight(log_mrf[["true"]], log_density[rep, 1])
                } else {
                    ops$weight(log_mrf[["false"]], log_density[rep, 2])
                }
                prod_val <- prod_val * term / norm
            }
            loglik <- loglik + log(prod_val) / length(k)
        }
    }
    loglik
}

.saint_loglikelihood_z <- function(model, prey_row, bait_col, rep, pre_calc) {
    log_mrf <- .saint_cell_log_mrf(model, prey_row, bait_col)
    log_density <- pre_calc[[prey_row, bait_col]][rep, ]
    if (model$Z[[prey_row, bait_col]][[rep]]) {
        log_mrf[["true"]] + log_density[[1]]
    } else {
        log_mrf[["false"]] + log_density[[2]]
    }
}

.saint_icm_z <- function(model, pre_calc) {
    for (prey_row in seq_len(model$nprey)) {
        for (bait_col in seq_len(model$nbait)) {
            for (rep in seq_len(model$n_rep_vec[[bait_col]])) {
                first <- .saint_loglikelihood_z(model, prey_row, bait_col, rep, pre_calc)
                z <- model$Z[[prey_row, bait_col]]
                z[[rep]] <- !z[[rep]]
                model$Z[[prey_row, bait_col]] <- z
                second <- .saint_loglikelihood_z(model, prey_row, bait_col, rep, pre_calc)
                if (first > second) {
                    z[[rep]] <- !z[[rep]]
                    model$Z[[prey_row, bait_col]] <- z
                }
            }
        }
    }
    model
}

.saint_gsum_mat <- function(model) {
    gsum <- matrix(0, nrow = model$nprey, ncol = model$nbait)
    for (prey_row in seq_len(model$nprey)) {
        for (bait_col in seq_len(model$nbait)) {
            gsum[prey_row, bait_col] <- .saint_gsum(model, prey_row, bait_col)
        }
    }
    gsum
}

# MRF part of the likelihood as a function of beta1 and gamma.
.saint_llik_mrf <- function(model, beta1, gamma, gsum_mat, empty) {
    loglik <- 0
    mrf_false <- exp(0)
    for (prey_row in seq_len(model$nprey)) {
        for (bait_col in seq_len(model$nbait)) {
            if (empty(model$test_mat[[prey_row, bait_col]])) {
                next
            }
            k <- model$Z[[prey_row, bait_col]]
            mrf_true <- exp(beta1 + gamma * gsum_mat[prey_row, bait_col])
            prod <- prod(ifelse(k, mrf_true, mrf_false) / (mrf_true + mrf_false))
            loglik <- loglik + log(prod) / model$n_rep_vec[[bait_col]]
        }
    }
    loglik
}

# Maximises the MRF likelihood over beta1 (gamma fixed at 0) or over beta1
# and gamma, and keeps the old values when they were better.
.saint_optimize_mrf <- function(model, with_gamma, optimizer, empty) {
    gsum_mat <- if (with_gamma) {
        .saint_gsum_mat(model)
    } else {
        matrix(0, nrow = model$nprey, ncol = model$nbait)
    }
    llik <- function(beta1, gamma) .saint_llik_mrf(model, beta1, gamma, gsum_mat, empty)
    old <- if (with_gamma) c(model$beta1, model$gamma) else model$beta1
    oldf <- llik(model$beta1, model$gamma)
    fn <- if (with_gamma) {
        function(x) -llik(x[[1]], x[[2]])
    } else {
        function(x) -llik(x[[1]], 0)
    }
    res <- .saint_optimize(
        fn,
        x0 = old,
        lb = if (with_gamma) c(-15, 0) else -15,
        ub = if (with_gamma) c(15, 10) else 15,
        optimizer = optimizer
    )
    maxf <- -res$value
    model$beta1 <- res$par[[1]]
    if (with_gamma) {
        model$gamma <- res$par[[2]]
        if (maxf <= llik(res$par[[1]], 0)) {
            model$gamma <- 0
        }
    }
    if (maxf < oldf) {
        model$beta1 <- old[[1]]
        model$gamma <- if (with_gamma) old[[2]] else model$gamma
    }
    model
}

# Minimises fn within the bounds lb, ub: NLopt COBYLA (as native
# SAINTexpress) or stats::optim L-BFGS-B.
.saint_optimize <- function(fn, x0, lb, ub, optimizer) {
    if (optimizer == "nloptr") {
        res <- nloptr::nloptr(
            x0 = x0,
            eval_f = fn,
            lb = lb,
            ub = ub,
            opts = list(algorithm = "NLOPT_LN_COBYLA", ftol_abs = 1e-4, maxeval = 100)
        )
        return(list(par = res$solution, value = res$objective))
    }
    res <- optim(
        par = x0,
        fn = fn,
        method = "L-BFGS-B",
        lower = lb,
        upper = ub,
        control = list(maxit = 100)
    )
    list(par = res$par, value = res$value)
}

.saint_calculate_score <- function(model, R, ops) {
    average <- matrix(0, nrow = model$nprey, ncol = model$nbait)
    maximum <- average
    odds <- average
    for (prey_row in seq_len(model$nprey)) {
        for (bait_col in seq_len(model$nbait)) {
            y <- model$test_mat[[prey_row, bait_col]]
            log_mrf <- .saint_cell_log_mrf(model, prey_row, bait_col)
            reps <- ops$rep_scores(model, prey_row, y, log_mrf)
            scores <- sort(reps$score, decreasing = TRUE)
            average[prey_row, bait_col] <- mean(scores[seq_len(min(R, length(scores)))])
            maximum[prey_row, bait_col] <- scores[[1]]
            odds[prey_row, bait_col] <- sort(reps$odds)[[1]]
        }
    }
    list(average = average, maximum = maximum, odds = odds)
}

.saint_logit_probability <- function(true_log, false_log) {
    log_scale <- max(true_log, false_log)
    denom <- log_scale +
        log(exp(true_log - log_scale) + exp(false_log - log_scale))
    exp(true_log - denom)
}

.saint_list_output <- function(model, parsed, scores, topo_scores, ops) {
    rows <- lapply(parsed$ui, function(ui) {
        prey_row <- ui$row
        bait_col <- ui$baitcol
        values <- ops$values(model, parsed, prey_row, bait_col)
        avg_p <- scores$average[prey_row, bait_col]
        topo_avg_p <- topo_scores$average[prey_row, bait_col]
        data.frame(
            Bait = ui$baitId,
            Prey = ui$preyId,
            PreyGene = ui$preyGeneId,
            quant = values$quant,
            sum = values$sum,
            avg = values$avg,
            NumReplicates = values$n,
            ctrl = values$ctrl,
            AvgP = avg_p,
            MaxP = scores$maximum[prey_row, bait_col],
            TopoAvgP = topo_avg_p,
            TopoMaxP = topo_scores$maximum[prey_row, bait_col],
            SaintScore = max(topo_avg_p, avg_p),
            odds = topo_scores$odds[prey_row, bait_col],
            FoldChange = values$fold_change,
            BFDR = .saint_bfdr(avg_p, scores$average),
            boosted_by = .saint_boosted_by(model, parsed, prey_row, bait_col),
            stringsAsFactors = FALSE
        )
    })
    if (!length(rows)) {
        return(data.frame())
    }
    res <- do.call(rbind, rows)
    numeric_cols <- c("sum", "avg", "AvgP", "MaxP", "TopoAvgP", "TopoMaxP", "SaintScore", "odds", "FoldChange", "BFDR")
    for (col in numeric_cols) {
        res[[col]] <- round(res[[col]], ops$digits)
    }
    generic <- match(names(res), names(ops$names))
    names(res)[!is.na(generic)] <- ops$names[generic[!is.na(generic)]]
    res
}

# Bayesian FDR: 1 - mean AvgP of all pairs scoring higher.
.saint_bfdr <- function(avg_p, average_score) {
    all_average_values <- as.vector(average_score)
    greater <- all_average_values > avg_p
    denom <- sum(greater)
    if (denom == 0) {
        return(0)
    }
    1 - sum(all_average_values[greater]) / denom
}

.saint_boosted_by <- function(model, parsed, prey_row, bait_col) {
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
