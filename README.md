# Replication data and code

This package reproduces the numerical results in the accepted version of
“Attention or Backlash: How English Protest Signs Influence Campaign Success”
from the frozen, merged analysis data. Raw-source merging and image processing
are intentionally outside this package.

## Article materials

> Lian, Jie, and Amanda Murdie. “Attention or Backlash: How English Protest
> Signs Influence Campaign Success.” *Journal of Conflict Resolution*.
> Accepted for publication.

- [Download the full online appendix](Lian_Murdie_2026_Full_Online_Appendix.docx).
- [Browse the article and appendix plots](publication_figures/README.md).
- Appendix SHA-256: `f2d922d27aefd8c6083f5fc16f7635f2d3ab4c181acd19bf4b623760713523b0`.
- [Protest-sign detector code and weights](https://github.com/LIANJie-Jason/Protest-Sign-Detector).

The online appendix is a byte-identical copy of the accepted
`Appendix_FULL.docx`. Its checksum is also stored in
`ONLINE_APPENDIX.sha256`.

The plot gallery contains the analytical figures appearing in the accepted
article and full online appendix. Publication-reference images are kept under
`publication_figures/`; code-generated main Figures 7–9 and Appendix Figures
C1–C2 remain under `output/figures/` with their numerical plot-data sidecars.

## Availability and rights

The authors have approved public distribution of this replication package and
the full online appendix. Raw protest photographs and image-level annotations
are not included. The frozen merged inputs incorporate variables derived from
NAVCO, V-Dem, QoG, WDI, AidData, Polity, and the article's original measurement
workflow; source components remain subject to their respective terms and
citation requirements. Public availability of this package does not transfer
rights in those upstream sources or in third-party media.

The accepted appendix and publication-figure images are provided at no charge
for non-commercial, no-derivatives use with citation of the accepted article,
consistent with [Sage's author archiving and re-use guidelines](https://us.sagepub.com/en-us/nam/journal-author-archiving-policies-and-re-use).

Please cite the accepted article when using the package. Code, data, appendix,
and third-party components retain the rights described here and in their source
documentation; no additional license grant is implied beyond lawful scholarly
replication and use.

## One-command replication

From the package directory, with R 4.3 or newer:

```sh
Rscript --vanilla run_all.R --profile=full
```

From another working directory, supply a valid relative or absolute path to
this `run_all.R`. The script discovers the package from its own file location;
no machine-specific path edits are required. Every pipeline stage runs in a
fresh `Rscript --vanilla` process and stops on the first failed dependency,
input-integrity check, model, renderer, or accepted-reference comparison.

Only one pipeline may run at a time. An atomic lock prevents overlapping runs;
if a prior process was killed, the runner automatically recovers the lock only
when its recorded PID is no longer alive. An unreadable/ambiguous lock is left
in place for safety and its owner is reported.

A development-only shorter run is available:

```sh
Rscript --vanilla run_all.R --profile=fast
```

The fast profile omits accepted stochastic-artifact materialization, inference
diagnostics, and Appendix C25–C26. Its validation status is always
`INCOMPLETE`, never publication `PASS`.

## Package layout

```text
replication_datacode/
├── run_all.R                    public entry point
├── code/R/                     publication R pipeline
│   ├── 00_run_all.R             fresh-process orchestrator
│   ├── 01_preflight.R           software and portability checks
│   ├── 02_validate_inputs.R     input hashes, shapes, schema, and keys
│   ├── 10–15_*.R                main and extension models
│   ├── 20–22_*.R                robustness and sensitivity models
│   ├── 23_materialize_*.R       verified accepted stochastic outputs
│   ├── 23–28 optional scripts   non-default stochastic recomputation tools
│   ├── 29_*.R                   inference diagnostics
│   ├── 30–36_*.R                main/appendix tables and figures
│   ├── 37_materialize_*.R       accepted display-table copies
│   ├── 40_validate_*.R          accepted-reference audit gate
│   ├── 41_write_*.R             checksummed output manifest
│   └── R/                       shared base-R infrastructure
├── data/                        frozen merged inputs, README + SHA-256 manifest
├── docs/                        accepted specification/exhibit inventory
├── provenance/                  optional C20 reconstruction audit
├── environment/                 minimum and tested software versions
├── reference/accepted/          immutable accepted numerical artifacts
├── publication_figures/         accepted article/appendix plot gallery
├── output/
│   ├── estimates/               fitted-result CSVs
│   ├── tables/main/             Tables 1–4 bodies
│   ├── tables/appendix/         Appendix C table bodies
│   ├── figures/                 main and appendix figures + plot data
│   └── diagnostics/             power and inference diagnostics
└── audit/                       logs, sessionInfo, comparisons, manifests
```

Model-object RDS files are temporary runtime intermediates. They are created
under `output/cache/` for downstream stages and removed after a successful
validation, so they are not shipped as duplicate publication artifacts.

Publication code is finalized under `code/R/`; temporary review/staging paths
are not included in the delivered package.

## Frozen inputs

The pipeline begins with three supplied CSVs:

- `data/df_final.csv`: 101 campaign-location observations (NAVCO 1.3 stack).
- `data/df_navco21_panel.csv`: 65 source-campaign-year observations (NAVCO 2.1
  stack).
- `data/h3_additional_moderators.csv`: 101 frozen campaign-period moderator
  observations used for accepted Appendix C20.

`data/input_manifest.csv` records the expected SHA-256, row/column count,
analytic key, required fields, and provenance for each file. Stage 02 verifies
all of these before any model is fitted. A changed byte, shape, schema, or
duplicate analytic key is fatal.

## Software environment

Required packages and minimum versions are listed in
`environment/required_packages.csv`; the environment used to verify the
package is in `environment/tested_versions.csv`. The pipeline never installs
packages. Preflight reports missing or too-old packages and writes the actual
versions plus complete `sessionInfo()` to `audit/`.

The figure stage additionally requires Poppler's `pdftoppm` executable. Its
R dependencies include `MASS`, `patchwork`, and `png`. Main Figures 7–9 are
rasterized with Poppler at 144 DPI and compared to the accepted PDFs at the
predeclared `0.05` mean absolute per-channel threshold. That comparison is
advisory because graphics devices and fonts vary; the numerical plot-data
sidecars are blocking and must match their frozen references within `1e-8`.
Figure 7 preserves accepted seed `20260506`; Figures 8 and 9 preserve the
accepted producer's shared seed `20260511`. These are frozen rendering
conditions and are separate from the 136-cell optional power/calibration seed
registry. Change them only when regenerating the derived figure-data references
and their provenance hashes together.

The verified environment is a cairo-capable R build on macOS. A Unix-style
SHA-256 utility (`shasum` or `sha256sum`) and Poppler `pdftoppm` are required;
Linux users should provide the same tools and a cairo-enabled graphics device.
The package is not presented as a stock-Windows command path.

## Validation and audit trail

`reference/accepted/reference_manifest.csv` freezes accepted reference files
with SHA-256 hashes and maps every regenerated numerical artifact to its
publication reference. Validation checks:

1. accepted-reference integrity;
2. output presence and production during the current run;
3. identical schema, labels, and missingness for blocking numerical CSVs; and
4. numeric agreement under the artifact-specific absolute/relative tolerance.

Table-body CSVs are convenience render inputs, not the inferential source of
truth. Their displayed strings/topology are compared and reported as advisory
evidence only. Differences in table layout or formatting do not block a clean
replication when the underlying model/term estimates, standard errors,
p-values, confidence intervals, sample/event counts, adjusted p-values,
diagnostics, power results, and plot data pass their numerical gates. This
package does not reproduce DOCX layout or journal typesetting.

After every table builder runs, stage 37 copies the 28 fast-profile or 32
full-profile accepted display bodies byte-exactly into `output/tables/` and
writes `display_materialization.csv`. This preserves the existing accepted
numbers/layout exactly while regenerated estimate CSVs, fresh diagnostics, and
fit-status contracts remain the blocking numerical evidence. The pre-copy
ledger records actual and accepted strings for every mismatch: 26 of 28
independent table builders match exactly; the three bounded mismatches are the
C4 H2.1 `M2s` constant and the C8 H2.2 `M1s` `REGCHANGE`/constant cells. Those
are nuisance terms in fits labeled `accepted_archival_nonconverged`, so the
accepted display cells are preserved explicitly rather than claimed as clean
independent reproduction.

`reference/accepted/exhibit_manifest.csv` is the executable publication-
coverage contract: exactly Tables 1–4, Figures 7–9, Appendix Tables C1–C26,
and Appendix Figures C1–C2 (35 exhibits, 42 components). Stage 03 also blocks
the deleted identification section, the stale H3 movement slot, obsolete
H2.1 k=6/k=14 multiplicity logic, any loss/duplication of C22's `NVH1` row,
and any regression in C25's accepted row count or column contract. Superseded
artifacts are excluded from the accepted namespace and final package.

Appendix C20 explicitly records its frozen V-Dem v15 moderator artifact and
lineage in `reference/accepted/provenance/C20_provenance.csv`. The optional
`provenance/validate_C20_frozen_input.R` script rechecks the recorded source
payload when `vdemdata` 15.0 is available; the publication run itself starts
from the hash-frozen merged input.

Deterministic reference sidecars that were constructed specifically for this
package—Figures 7–9 plot data, C9/C12 stable-ID coefficient outputs, and C26
full-precision Cox PH diagnostics—are explicitly labeled
`frozen_derived_reference` in
`reference/accepted/provenance/derived_numeric_references.csv`. They are
hash-linked regression oracles anchored to the accepted display/PDF, not
misrepresented as separately published source files.

The corrected accepted Appendix C22 reference is
`C22_nonviolent_final_body.csv` with 31 rows, including NVH1. The superseded
earlier C22 body is excluded from the validation registry.

The main audit products are:

- `audit/input_validation.csv`
- `audit/reference_validation.csv`
- `audit/validation_summary.csv`
- `audit/output_manifest.csv`
- `audit/code_manifest.csv`
- `audit/runtime_versions.csv` and `audit/sessionInfo.txt`
- `audit/pipeline_status.csv`

A successful full run requires every registered full-profile artifact to be
fresh and clean. Rendered PDFs are checked for successful creation; their
underlying table bodies, estimate CSVs, and plot-data values are the numerical
validation surface.

## Accepted-fit disclosures

The package preserves accepted numerical values even when a post-acceptance
fit-status audit identifies nonconvergence or monotone likelihood. Those fits
are reproduced for publication fidelity but classified separately from clean
replications. Across the six cached main/extension bundles, the cross-family
audit covers all supported `survreg`, `coxph`, `logistf`, and `geeglm` objects
(including GEE error codes):

- strict-success Weibull AFT fits with co-divergent `REGCHANGE` and intercept
  terms: H1 `M1s`; H2.1 `M1s/M2s`; H2.2 `M1s–M4s`; and H3 `M1s–M4s`;
- strict-success Weibull H2.1 extensions `M1qs/M2qs`, with the same divergent
  nuisance-term signature;
- strict-success Cox fits with a divergent `REGCHANGE` term: H1 `M4s`; H2.1
  `M7s/M8s`; H2.2 `M13s–M16s`; H3 `M13s–M16s`; and H2.1 extensions
  `M7qs/M8qs`;
- limited-success Cox fits `M7l/M8l`, which hit the Cox iteration limit;
- H2.1 extension `M9q`, which fails `logistf` convergence thresholds; and
- H3 robustness `M20_GEE`, which reports a nonzero `geese$error`.

These 30 fits are labeled `accepted_archival_nonconverged` in
`output/diagnostics/accepted_fit_status_all.csv`; the remaining 110 supported
fits are labeled `accepted_replicated`. The narrower logistf and H3-GEE audit
CSVs are retained for estimator-specific detail. The pipeline blocks any silent
upgrade, omission, or change in this accepted warning set.

C9/C12 quadratic models and the C22/C24 robustness fits are outside those six
bundles; their exported terms are instead validated directly against frozen
long-form numerical references at `1e-8`.

## Accepted codestack artifacts

The accepted power and bootstrap numbers already exist in the authoritative
`RR cowork/RR_datacode/empirical clean` stack. Those eleven CSVs were verified
byte-for-byte against their copies under `reference/accepted/model_results/`.
The full publication run materializes those verified numbers into `output/`
and records source/output hashes in
`output/diagnostics/stochastic_materialization.csv`; it does not spend hours
resampling numbers that are already final and accurate.

Their source paths, producer-script hashes, and artifact hashes are frozen in
`reference/accepted/provenance/accepted_stochastic_artifacts.csv` with lineage
class `accepted_codestack_copy`. This is an explicit source-copy claim, not an
independent Monte Carlo re-estimation claim. The C25 display and Cox quarantine
artifacts retain their narrower `frozen_passthrough` disclosures.

The cleaned scripts `23_weibull_calibration.R` through
`28_power_summaries.R` remain available as optional stochastic recomputation
tools, but are not part of the one-command publication pipeline. Their smoke
tests and accepted/fresh seed contracts are retained because they provide the
portable recomputation code behind the copied stochastic artifacts, while the
default publication path preserves the already accepted numbers exactly.
Optional runs always use `_smoke` or `_recomputed` filenames and cannot
overwrite the accepted-codestack outputs.
H2.1 Cox M7/M8 power remains excluded because the counting-process simulation
design was invalid; no Cox power claim is made.
