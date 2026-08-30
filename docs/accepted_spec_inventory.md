# Accepted Statistical Specification Inventory

## Scope

This package begins from the two frozen merged analysis files. It does not
rebuild the Google Images, CV/NLP, NAVCO, V-Dem, QoG, WDI, AidData, or Polity
merges. Its purpose is to reproduce every statistical number displayed in the
accepted main article and Appendix C.

Exact Word-table layout is outside scope. Table-body CSVs are convenience
views; the blocking replication surface is the underlying model/check output
(model ID, term, coefficient, standard error, p-value, adjusted p-value,
confidence interval, sample/event count, diagnostics, and plotted values).

## Frozen inputs

| File | Shape | Analytic unit |
|---|---:|---|
| `data/df_final.csv` | 101 x 66 | campaign-location |
| `data/df_navco21_panel.csv` | 65 x 79 | NAVCO 2.1 source-campaign-year |

Appendix C20 additionally requires the campaign-period mean of V-Dem
`v2x_frassoc_thick`, which is absent from the accepted `df_final.csv`. The
package therefore freezes one keyed, pre-merged auxiliary input containing
that moderator. No live V-Dem download or merge occurs during replication.

## Main article

| Exhibit | Accepted models | Producer |
|---|---|---|
| Table 1 / H1 | `M1`, `M2`, `M4`, `M5` | `10_fit_h1.R`, then main-table builder |
| Table 2 / H2.1 | `M1`, `M2`, `M3`, `M4`, `M7`, `M8`, `M9`, `M10` | `11_fit_h21.R` |
| Table 3 / H2.2 | `M1/M5`, `M2/M6`, `M3/M7`, `M4/M8` | `12_fit_h22.R` |
| Table 4 / H3 | `M1/M5`, `M2/M6`, `M3/M7`, `M4/M8` | `13_fit_h3.R` |
| Figure 7 | H1 `M1` Weibull AFT and `M4` Cox PH | H1 figure builder |
| Figure 8 | H2.2 trade `M2` Weibull and `M6` Firth | H2.2 figure builder |
| Figure 9 | H3 `M1/M5` CSO repression and `M3/M7` political liberty | H3 figure builder |

H1 multiplicity uses a fixed 10-model family. H2.1 uses Holm `k=8` for the
main panel and BH `k=12` for the core appendix. H2.2 and H3 each use Holm
`k=8` for the main panel and BH `k=32` for the core appendix. Labeled
robustness rows are outside those families.

## Appendix C

| Table | Accepted contents |
|---|---|
| C1 | H1 strict-success and ordered models |
| C2 | H1 distributional and clustering sensitivity |
| C3 | H2.1 curvilinear concession extension |
| C4 | H2.1 strict-success twins |
| C5 | H2.1 limited-success models |
| C6 | H2.1 ordered, recurrent-Cox, and GEE rows |
| C7 | H2.2 campaign-level broad-success models |
| C8 | H2.2 campaign-level strict-success models |
| C9 | H2.2 campaign-level squared-English interactions |
| C10 | H2.2 campaign-year broad-success models |
| C11 | H2.2 campaign-year strict-success models |
| C12 | H2.2 campaign-year squared-English interactions |
| C13 | H2.2 ordered models |
| C14 | H2.2 recurrent-Cox and GEE robustness |
| C15 | H3 campaign-level strict-success models |
| C16 | H3 ordered models |
| C17 | H3 campaign-year broad-success models |
| C18 | H3 campaign-year strict-success models |
| C19 | H3 recurrent-Cox and GEE robustness |
| C20 | H3 political-liberty parity and freedom-of-association extension |
| C21 | Strict-success robustness across all families |
| C22 | Nonviolent-only robustness, including exactly one `NVH1` row |
| C23 | Ten `sensemakr` robustness-value rows |
| C24 | Eight ongoing-campaign sensitivity models |
| C25 | Twenty-eight reader-facing power rows; power column precedes false-positive rate; no `Source` column |
| C26 | 34 Cox PH rows, 3 Cox-power quarantine rows, and 11 Weibull-calibration rows |

Appendix Figures C1 and C2 are the H1 curvature and H3 interaction forest
plots. A duplicate identification section from an earlier draft is not part of
the accepted appendix and is intentionally omitted.

## Expected output sizes

- Model results: H1 14; H2.1 24; H2.2 48; H3 48; H2.1-curvilinear 20;
  H2.1-limited 8; H3-additional 4; ongoing 96; nonviolent 289;
  `sensemakr` 10.
- Main table bodies: 28, 28, 34, and 32 rows for Tables 1-4.
- C1-C20 bodies: 29, 27, 35, 29, 29, 29, 34, 34, 44, 34, 34, 44,
  34, 34, 37, 39, 35, 35, 35, and 13 rows.
- C21-C24: 44, 31, 10, and 8 rows. C22's first row is the accepted
  `NVH1` curvilinear nonviolent specification.
- C25: 28 rows. C26 panels: 34, 3, and 11 rows.
