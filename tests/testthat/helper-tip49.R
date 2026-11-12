# TIP49 subset (inst/extdata/tip49) and cached run_saint() results, shared
# between test files.

tip49_dir <- function() {
    system.file("extdata", "tip49", package = "SAINTexpress")
}

read_tip49 <- function() {
    read_dat <- function(file) {
        utils::read.delim(
            file.path(tip49_dir(), file),
            header = FALSE,
            stringsAsFactors = FALSE
        )
    }
    list(
        inter = read_dat("inter.dat"),
        prey = read_dat("prey.dat"),
        bait = read_dat("bait.dat")
    )
}

read_tip49_native <- function(mode) {
    utils::read.delim(
        file.path(tip49_dir(), paste0("list_", mode, ".txt")),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

# Known prey-prey links for TIP49: members of the INO80 and SRCAP complexes.
tip49_p2p_mapping <- function(si) {
    complexes <- list(
        c(
            "RUVBL1",
            "RUVBL2",
            "ACTR5",
            "ACTR8",
            "INO80B",
            "INO80C",
            "INO80D",
            "INO80E",
            "NFRKB",
            "TFPT",
            "MCRS1",
            "UCHL5",
            "YY1",
            "ACTL6A"
        ),
        c("SRCAP", "ZNHIT1", "VPS72", "H2AFZ", "DMAP1", "ACTL6A", "RUVBL1", "RUVBL2")
    )
    prey_ids <- si$prey[[1]]
    lapply(prey_ids, function(prey_id) {
        in_complex <- vapply(complexes, function(members) prey_id %in% members, logical(1))
        linked <- setdiff(unique(unlist(complexes[in_complex])), prey_id)
        sort(match(linked, prey_ids))
    })
}

tip49_cache <- new.env(parent = emptyenv())

# run_saint() on TIP49, computed once per argument set.
run_tip49 <- function(mode, optimizer = "base", p2p = FALSE, ...) {
    key <- paste(mode, optimizer, p2p, paste(names(list(...)), list(...), collapse = ","))
    if (is.null(tip49_cache[[key]])) {
        si <- read_tip49()
        p2p_mapping <- if (p2p) tip49_p2p_mapping(si) else NULL
        tip49_cache[[key]] <- run_saint(
            si,
            mode = mode,
            optimizer = optimizer,
            p2p_mapping = p2p_mapping,
            ...
        )
    }
    tip49_cache[[key]]
}
