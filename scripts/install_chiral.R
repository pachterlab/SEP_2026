#!/usr/bin/env Rscript
# Install CHIRAL into this R environment, never the user's system R library.
revision <- "1f78a9a60b6d18facfedf5923d43bae9fa79f0b0"
lib <- .Library
description <- file.path(lib, "CHIRAL", "DESCRIPTION")
installed <- if (file.exists(description)) read.dcf(description) else NULL
if (is.null(installed) || !"RemoteSha" %in% colnames(installed) ||
    installed[1, "RemoteSha"] != revision) {
    remotes::install_github(
        paste0("naef-lab/CHIRAL@", revision), subdir = "Pkg/CHIRAL",
        lib = lib, dependencies = FALSE, upgrade = "never", build_vignettes = FALSE
    )
}
stopifnot(requireNamespace("CHIRAL", lib.loc = lib, quietly = TRUE))
cat("CHIRAL", as.character(packageVersion("CHIRAL", lib.loc = lib)),
    "installed in", lib, "\n")
