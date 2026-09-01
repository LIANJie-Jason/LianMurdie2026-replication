# Frozen Analysis Inputs

The replication package begins from the accepted merged analysis data. Raw
image collection, CV/NLP measurement, NAVCO/V-Dem/QoG merging, and recoding
are outside this package's execution scope.

The authors intend these frozen merged inputs for scholarly replication of the
accepted article. Upstream variables remain subject to their original source
terms and citation requirements; this package does not redistribute raw protest
images or image-level annotations.

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

## Source provenance and terms

The merged files contain campaign-level or campaign-year derivatives, not raw
provider datasets. Cite the article and the applicable sources below.

| Source | Release used | Terms relevant to this package |
|---|---|---|
| NAVCO | [NAVCO 1.3](https://doi.org/10.7910/DVN/ON9XND) and [NAVCO 2.1](https://doi.org/10.7910/DVN/MHOXDV) | Cite the dataset versions; NAVCO 2.1 is deposited under CC0 1.0. |
| V-Dem | [Country-year v14 and v15](https://www.v-dem.net/data/dataset-archive/) | CC BY-SA 4.0; attribution and share-alike apply to V-Dem-derived fields. |
| QoG | [Standard Dataset, January 2023](https://doi.org/10.18157/qogstdjan23) | The economic and institutional fields were accessed through this compilation. QoG permits academic, non-commercial use but states that redistribution is not allowed. |
| World Bank | [World Development Indicators](https://databank.worldbank.org/source/world-development-indicators) | WDI fields were delivered through QoG; World Bank open data are generally CC BY 4.0. |
| AidData | Core Research Release 3.1; [data-use terms](https://www.aiddata.org/pages/data-user-guide) | Cite AidData and the original providers; use is limited to private or personal, non-commercial purposes. |
| Polity | [Polity5, 1800–2018](https://www.systemicpeace.org/inscrdata.html) | The provider requires prior written permission for reproduction or redistribution. |
| Fraser Institute | Economic Freedom of the World field included in QoG January 2023; the article cites [Gwartney and Murphy (2024)](https://www.fraserinstitute.org/studies/economic-freedom-of-the-world-2024-annual-report) | The trade-freedom field was delivered through QoG; cite the source and QoG. |
| Authoritarian Regime Dataset | Hadenius and Teorell (2007), as [distributed through QoG](https://datafinder.qog.gu.se/dataset/ht) | The colonial-history indicator was delivered through QoG; cite the original source and QoG. |

The authors' distribution approval covers materials they control and does not
override upstream terms. Public redistribution of the merged CSVs therefore
requires source permission or repository access conditions that honor the QoG,
AidData, and Polity restrictions.
