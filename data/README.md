# Frozen Analysis Inputs

The replication package begins from the accepted merged analysis data. Raw
image collection, CV/NLP measurement, NAVCO/V-Dem/QoG merging, and recoding
are outside this package's execution scope.

The authors have approved public distribution of these frozen merged inputs
for replication of the accepted article. Upstream variables remain subject to
their original source terms and citation requirements; this package does not
redistribute raw protest images or image-level annotations.

- `df_final.csv` is the accepted NAVCO 1.3 campaign-location file.
- `df_navco21_panel.csv` is the accepted NAVCO 2.1 source-campaign-year file.
- `h3_additional_moderators.csv` supplies only the frozen campaign-period mean
  of V-Dem v15 `v2x_frassoc_thick` required by accepted Appendix Table C20.
  It is a separately frozen, keyed auxiliary input used to reproduce that
  accepted analysis. Reconstructing it from upstream V-Dem releases is outside
  this package's execution scope; its source attribution and construction rule
  are documented in
  `reference/accepted/provenance/C20_provenance.csv`.

`input_manifest.csv` records immutable hashes, shapes, and analytic units. The
preflight script stops before estimation if any input has changed.
