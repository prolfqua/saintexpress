#' Validate a SAINT input list
#'
#' Checks that `si` is valid input for [run_saint()]: if `validate_saint_input()`
#' passes, `run_saint()` accepts the same `si` and `mode`. The checks are the
#' ones listed in the section "Input tables" below, applied after the columns
#' have been matched by name or position.
#'
#' @param si A named list with the data frames `inter`, `prey` and `bait`, in
#'   any order. See the section "Input tables".
#' @param mode `"spc"` for spectral counts, `"int"` for intensities. Sets the
#'   allowed range of `quant`.
#' @return `TRUE` invisibly on success; otherwise stops with a message naming
#'   the table and column at fault.
#' @inheritSection run_saint Input tables
#' @examples
#' si <- list(
#'   inter = data.frame(
#'     ipId = c("IP1", "IP2", "IP3"),
#'     baitId = c("BaitA", "BaitA", "Ctrl"),
#'     preyId = c("Prey1", "Prey1", "Prey1"),
#'     quant = c(12, 9, 1)
#'   ),
#'   prey = data.frame(preyId = "Prey1", preyLength = 500, preyGeneId = "GENE1"),
#'   bait = data.frame(
#'     ipId = c("IP1", "IP2", "IP3"),
#'     baitId = c("BaitA", "BaitA", "Ctrl"),
#'     CorT = c("T", "T", "C")
#'   )
#' )
#' validate_saint_input(si, mode = "spc")
#' @export
validate_saint_input <- function(si, mode = c("spc", "int")) {
    mode <- match.arg(mode)
    .saint_check_list(si)
    tables <- .saint_resolve_tables(si)
    .saint_check_ids(tables$inter, c("ipId", "baitId", "preyId"), "inter")
    .saint_check_ids(tables$prey, c("preyId", "preyGeneId"), "prey")
    .saint_check_ids(tables$bait, c("ipId", "baitId", "CorT"), "bait")
    .saint_check_quant(tables$inter$quant, mode)
    .saint_check_unique(tables$prey$preyId, "prey", "preyId")
    .saint_check_unique(tables$bait$ipId, "bait", "ipId")
    pair <- anyDuplicated(paste(tables$inter$ipId, tables$inter$preyId, sep = "\r"))
    if (pair) {
        stop(
            "inter: prey ",
            tables$inter$preyId[[pair]],
            " appears more than once in IP ",
            tables$inter$ipId[[pair]],
            "."
        )
    }
    .saint_check_cort(tables$bait$CorT)
    .saint_check_links(tables)
    invisible(TRUE)
}

.saint_check_list <- function(si) {
    elements <- c("inter", "prey", "bait")
    named_list <- is.list(si) && !is.data.frame(si) && length(si) == 3L
    if (!named_list || !setequal(names(si), elements)) {
        stop("si must be a named list with elements 'inter', 'prey', 'bait'.")
    }
    for (nm in elements) {
        if (!is.data.frame(si[[nm]])) {
            stop("si$", nm, " must be a data.frame.")
        }
        if (nrow(si[[nm]]) == 0L) {
            stop("si$", nm, " has no rows.")
        }
    }
}

.saint_check_quant <- function(quant, mode) {
    if (!is.numeric(quant) || !all(is.finite(quant))) {
        stop("inter: 'quant' must be numeric and finite, without missing values.")
    }
    if (mode == "spc" && any(quant < 0 | quant != round(quant))) {
        stop("inter: spectral counts in 'quant' must be whole numbers >= 0.")
    }
    if (mode == "int" && any(quant <= 0)) {
        stop("inter: intensities in 'quant' must be > 0.")
    }
}

.saint_check_cort <- function(cort) {
    if (!all(cort %in% c("T", "C"))) {
        stop("bait: 'CorT' must be 'T' (test) or 'C' (control).")
    }
    if (!all(c("T", "C") %in% cort)) {
        stop("bait: needs at least one test ('T') and one control ('C') IP.")
    }
}

.saint_check_links <- function(tables) {
    inter <- tables$inter
    bait <- tables$bait
    .saint_check_known(inter$ipId, bait$ipId, "inter", "ipId", "bait")
    .saint_check_known(inter$preyId, tables$prey$preyId, "inter", "preyId", "prey")
    bait_of_ip <- bait$baitId[match(inter$ipId, bait$ipId)]
    mismatch <- as.character(inter$baitId) != as.character(bait_of_ip)
    if (any(mismatch)) {
        stop(
            "inter: baitId differs from the bait table for ipId ",
            inter$ipId[which(mismatch)[[1]]],
            "."
        )
    }
}

# Matches the required columns of a table by name when all are present,
# otherwise by position. A column already named like a required column must
# then sit at that column's position, so named columns are never reassigned.
.saint_resolve_columns <- function(df, required, table) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (all(required %in% names(df))) {
        return(df)
    }
    if (ncol(df) < length(required)) {
        stop(
            table,
            ": needs columns ",
            paste(required, collapse = ", "),
            ", by name or in this order."
        )
    }
    position <- match(required, names(df))
    misplaced <- !is.na(position) & position != seq_along(required)
    if (any(misplaced)) {
        stop(
            table,
            ": column '",
            required[misplaced][[1]],
            "' is not at position ",
            which(misplaced)[[1]],
            ". Name all of ",
            paste(required, collapse = ", "),
            ", or give them in this order."
        )
    }
    names(df)[seq_along(required)] <- required
    df
}

# Resolved tables, with all id columns as character.
.saint_resolve_tables <- function(si) {
    as_character <- function(df, columns) {
        df[columns] <- lapply(df[columns], as.character)
        df
    }
    inter <- .saint_resolve_columns(si$inter, c("ipId", "baitId", "preyId", "quant"), "inter")
    prey <- .saint_normalize_prey_table(as.data.frame(si$prey, stringsAsFactors = FALSE))
    bait <- .saint_resolve_columns(si$bait, c("ipId", "baitId", "CorT"), "bait")
    list(
        inter = as_character(inter, c("ipId", "baitId", "preyId")),
        prey = as_character(prey, c("preyId", "preyGeneId")),
        bait = as_character(bait, c("ipId", "baitId", "CorT"))
    )
}

# Prey tables named preyId, preyGeneId are kept; otherwise column 1 is
# preyId and columns 2-3 are preyLength and/or preyGeneId.
.saint_normalize_prey_table <- function(prey) {
    if (all(c("preyId", "preyGeneId") %in% names(prey))) {
        return(prey)
    }
    slot <- c(preyId = 1L, preyLength = 2L, preyGeneId = if (ncol(prey) == 2L) 2L else 3L)
    position <- match(names(slot), names(prey))
    misplaced <- !is.na(position) & position != slot
    if (any(misplaced)) {
        stop(
            "prey: column '",
            names(slot)[misplaced][[1]],
            "' is not at position ",
            slot[misplaced][[1]],
            ". Name both preyId and preyGeneId, or give ",
            "preyId, preyLength, preyGeneId in this order."
        )
    }
    out <- data.frame(
        preyId = as.character(prey[[1]]),
        preyLength = 0,
        preyGeneId = as.character(prey[[1]]),
        stringsAsFactors = FALSE
    )
    if (ncol(prey) == 2) {
        prey_length <- .saint_parse_length_column(prey[[2]])
        if (prey_length$is_length) {
            out$preyLength <- prey_length$values
        } else {
            out$preyGeneId <- as.character(prey[[2]])
        }
    } else if (ncol(prey) > 2) {
        out$preyLength <- as.numeric(prey[[2]])
        out$preyGeneId <- as.character(prey[[3]])
    }
    out
}

.saint_parse_length_column <- function(x) {
    if (is.numeric(x)) {
        return(list(is_length = !anyNA(x), values = x))
    }
    x <- trimws(as.character(x))
    numeric_pattern <- "^[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?$"
    is_length <- all(nzchar(x)) && all(grepl(numeric_pattern, x))
    list(is_length = is_length, values = if (is_length) as.numeric(x) else numeric())
}

.saint_check_ids <- function(df, columns, table) {
    for (column in columns) {
        x <- df[[column]]
        if (!is.atomic(x) || anyNA(x)) {
            stop(table, ": '", column, "' must be an id column without missing values.")
        }
    }
}

.saint_check_unique <- function(x, table, column) {
    if (anyDuplicated(x)) {
        stop(table, ": duplicated ", column, " ", x[anyDuplicated(x)], ".")
    }
}

.saint_check_known <- function(x, known, table, column, other) {
    unknown <- setdiff(x, known)
    if (length(unknown)) {
        stop(table, ": ", column, " ", unknown[[1]], " is not in the ", other, " table.")
    }
}

.saint_check_count <- function(x, name) {
    scalar <- is.numeric(x) && length(x) == 1L
    if (!scalar || !isTRUE(is.finite(x) && x >= 1 && x == round(x))) {
        stop(name, " must be a single positive integer.")
    }
}

.saint_check_p2p_mapping <- function(p2p_mapping, nprey) {
    if (is.null(p2p_mapping)) {
        return(invisible(TRUE))
    }
    valid_links <- function(links) {
        is.numeric(links) &&
            !anyNA(links) &&
            all(links == round(links)) &&
            all(links >= 1 & links <= nprey)
    }
    if (
        !is.list(p2p_mapping) ||
            length(p2p_mapping) != nprey ||
            !all(vapply(p2p_mapping, valid_links, logical(1)))
    ) {
        stop(
            "p2p_mapping must be NULL or a list with one integer vector per ",
            "row of the prey table, holding the row numbers (1 to ",
            nprey,
            ") of the linked preys."
        )
    }
    invisible(TRUE)
}
