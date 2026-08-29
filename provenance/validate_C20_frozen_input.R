#!/usr/bin/env Rscript

# Optional provenance audit for the one auxiliary field needed by Appendix C20.
# The publication replication starts from the frozen merged inputs and does not
# run this script.  Run it only to independently reconstruct and validate
# data/h3_additional_moderators.csv from the recorded V-Dem v15 package payload.

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

required_packages <- c("vdemdata", "countrycode")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) {
  stop(sprintf("Optional C20 provenance audit requires: %s",
               paste(missing_packages, collapse = ", ")), call. = FALSE)
}
if (!identical(as.character(utils::packageVersion("vdemdata")), "15.0")) {
  stop("C20 provenance is frozen to vdemdata 15.0", call. = FALSE)
}

sha256_file <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  command <- if (nzchar(Sys.which("shasum"))) "shasum" else "sha256sum"
  command_args <- if (identical(command, "shasum")) c("-a", "256", shQuote(path)) else shQuote(path)
  output <- system2(command, command_args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop(paste(output, collapse = " "), call. = FALSE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
max_diff <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (!any(keep)) return(NA_real_)
  max(abs(x[keep] - y[keep]))
}

df_path <- file.path(root, "data", "df_final.csv")
frozen_path <- file.path(root, "data", "h3_additional_moderators.csv")
df <- read.csv(df_path, stringsAsFactors = FALSE, check.names = FALSE)
frozen <- read.csv(frozen_path, stringsAsFactors = FALSE, check.names = FALSE)
key <- c("CAMPAIGN", "LOCATION")
if (anyDuplicated(df[key]) || anyDuplicated(frozen[key])) {
  stop("C20 campaign-location keys must be unique", call. = FALSE)
}

vdem_columns <- c("COWcode", "year", "v2x_civlib", "v2x_clpol", "v2x_frassoc_thick")
vdem <- as.data.frame(vdemdata::vdem[, vdem_columns])
vdem$iso3 <- suppressWarnings(countrycode::countrycode(
  vdem$COWcode, origin = "cown", destination = "iso3c"
))
vdem$iso3[vdem$COWcode == 345] <- "SRB"
vdem$iso3[vdem$COWcode == 347] <- "XKX"

reconstructed <- lapply(seq_len(nrow(df)), function(index) {
  start_year <- max(2009L, as.integer(df$BYEAR[[index]]))
  end_year <- min(2019L, as.integer(df$EYEAR[[index]]))
  if (!is.finite(start_year) || !is.finite(end_year) || start_year > end_year) {
    stop(sprintf("Invalid campaign window at input row %d", index), call. = FALSE)
  }
  rows <- vdem$iso3 == df$iso3[[index]] & vdem$year >= start_year & vdem$year <= end_year
  data.frame(
    CAMPAIGN = df$CAMPAIGN[[index]], LOCATION = df$LOCATION[[index]],
    iso3 = df$iso3[[index]], BYEAR = df$BYEAR[[index]], EYEAR = df$EYEAR[[index]],
    v2x_civlib_check = mean_na(vdem$v2x_civlib[rows]),
    v2x_clpol_check = mean_na(vdem$v2x_clpol[rows]),
    v2x_frassoc_thick = mean_na(vdem$v2x_frassoc_thick[rows]),
    stringsAsFactors = FALSE
  )
})
reconstructed <- do.call(rbind, reconstructed)

order_df <- do.call(order, df[key])
order_reconstructed <- do.call(order, reconstructed[key])
order_frozen <- do.call(order, frozen[key])
df <- df[order_df, , drop = FALSE]
reconstructed <- reconstructed[order_reconstructed, , drop = FALSE]
frozen <- frozen[order_frozen, , drop = FALSE]
if (!identical(df[key], reconstructed[key]) || !identical(df[key], frozen[key])) {
  stop("C20 reconstruction did not preserve campaign-location membership", call. = FALSE)
}

civlib_diff <- max_diff(reconstructed$v2x_civlib_check, df$v2x_civlib)
clpol_diff <- max_diff(reconstructed$v2x_clpol_check, df$v2x_clpol)
frassoc_diff <- max_diff(reconstructed$v2x_frassoc_thick, frozen$v2x_frassoc_thick)
if (any(!is.finite(c(civlib_diff, clpol_diff, frassoc_diff))) ||
    any(c(civlib_diff, clpol_diff, frassoc_diff) > 1e-8)) {
  stop(sprintf("C20 provenance mismatch: civlib=%.12g, clpol=%.12g, frassoc=%.12g",
               civlib_diff, clpol_diff, frassoc_diff), call. = FALSE)
}

package_root <- find.package("vdemdata")
audit <- data.frame(
  field = c(
    "audit_status", "vdemdata_version", "vdemdata_payload_sha256",
    "vdemdata_description_sha256", "country_mapping", "country_overrides",
    "campaign_window", "aggregation", "campaign_rows", "frozen_input_sha256",
    "max_abs_diff_v2x_civlib", "max_abs_diff_v2x_clpol",
    "max_abs_diff_v2x_frassoc_thick"
  ),
  value = c(
    "PASS", as.character(utils::packageVersion("vdemdata")),
    sha256_file(file.path(package_root, "data", "Rdata.rdb")),
    sha256_file(file.path(package_root, "DESCRIPTION")),
    "countrycode COW numeric to ISO3", "COW 345=SRB; COW 347=XKX",
    "max(2009,BYEAR) through min(2019,EYEAR), inclusive",
    "arithmetic mean over nonmissing campaign-years", nrow(reconstructed),
    sha256_file(frozen_path), format(civlib_diff, scientific = TRUE, digits = 17),
    format(clpol_diff, scientific = TRUE, digits = 17),
    format(frassoc_diff, scientific = TRUE, digits = 17)
  ),
  stringsAsFactors = FALSE
)
write.csv(audit, file.path(root, "provenance", "C20_validation.csv"),
          row.names = FALSE, na = "")
cat("C20 frozen-input provenance: PASS\n")
