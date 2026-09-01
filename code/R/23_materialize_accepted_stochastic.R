root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(root)) {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  file_arg <- gsub("~+~", " ", file_arg, fixed = TRUE)
  root <- normalizePath(file.path(dirname(file_arg), "..", ".."), mustWork = TRUE)
}
source(file.path(root, "code", "R", "00_setup.R"))
P <- init_replication(root)

contract_path <- file.path(P$reference, "provenance", "accepted_stochastic_artifacts.csv")
rep_assert_file(contract_path, "Accepted stochastic-artifact lineage")
contract <- read.csv(contract_path, stringsAsFactors = FALSE, check.names = FALSE)
rep_assert_columns(
  contract,
  c("artifact", "reference_path", "sha256", "producer_path", "producer_sha256",
    "accepted_producer_sha256", "lineage_class"),
  "accepted stochastic-artifact lineage"
)
rep_assert(nrow(contract) == 11L, "Accepted stochastic lineage must contain 11 artifacts")
rep_assert(all(contract$lineage_class == "accepted_codestack_copy"),
           "Every materialized stochastic artifact must be classified accepted_codestack_copy")

raw_power <- c("power_h1.csv", "power_h21.csv", "power_h22.csv", "power_h3.csv")
destination_relative <- ifelse(
  contract$artifact %in% raw_power,
  file.path("output", "estimates", contract$artifact),
  file.path("output", "diagnostics", contract$artifact)
)

rows <- lapply(seq_len(nrow(contract)), function(index) {
  source_path <- file.path(P$root, contract$reference_path[[index]])
  output_path <- file.path(P$root, destination_relative[[index]])
  rep_assert_file(source_path, sprintf("Accepted stochastic artifact %s", contract$artifact[[index]]))
  source_hash <- rep_sha256_file(source_path)
  rep_assert(identical(source_hash, contract$sha256[[index]]),
             "Accepted stochastic hash mismatch for %s", contract$artifact[[index]])
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(source_path, output_path, overwrite = TRUE, copy.mode = TRUE,
                      copy.date = FALSE)
  rep_assert(copied, "Could not materialize %s", contract$artifact[[index]])
  output_hash <- rep_sha256_file(output_path)
  status <- if (identical(output_hash, source_hash)) "PASS" else "FAIL"
  data.frame(
    artifact = contract$artifact[[index]],
    reference_path = contract$reference_path[[index]],
    output_path = destination_relative[[index]],
    reference_sha256 = source_hash,
    output_sha256 = output_hash,
    lineage_class = contract$lineage_class[[index]],
    status = status,
    stringsAsFactors = FALSE
  )
})
lineage <- do.call(rbind, rows)
write.csv(lineage, file.path(P$diagnostics, "stochastic_materialization.csv"),
          row.names = FALSE, na = "")
rep_assert(all(lineage$status == "PASS"),
           "One or more accepted stochastic artifacts failed byte-exact materialization")
cat(sprintf("Materialized %d accepted stochastic artifacts byte-exactly from the verified codestack copy.\n",
            nrow(lineage)))
