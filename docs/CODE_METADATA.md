# SoftwareX code metadata: amr-clonalshare

The row labels and ordering follow the SoftwareX original-software-publication
template version 6 (March 2026). The product version is `1.0.0` and there is no
candidate revision: a build from the same commit produces the same version.
Three public fields depend on actions the author has to take and fail closed
rather than pointing at an untagged or mutable location.

| Nr | Code metadata description | Metadata |
|---|---|---|
| C1 | Current code version | 1.0.0 |
| C2 | Permanent link to code/repository used for this code version | **Author action before submission.** The private GitHub repository and the `v1.0.0` tag are prepared by `scripts/prepare_release.py`; the permanent link is the Zenodo DOI minted from that tag, with deposit metadata in `.zenodo.json`. |
| C3 | Legal code license | MIT License (`LICENSE` and the SoftwareX-required `Licence.txt` are byte-identical). |
| C4 | Code versioning system used | git |
| C5 | Software code languages, tools and services used | Python; YAML; pytest; GitHub Actions. |
| C6 | Compilation requirements, operating environments and dependencies | No compilation step. Python 3.10-3.12; NumPy >= 1.21, pandas >= 1.5, SciPy >= 1.7, scikit-learn >= 1.0, PyYAML >= 5.4, pandera >= 0.18. Optional DendroPy >= 4.5 for phylogenetic tests and Matplotlib >= 3.5 for figures. |
| C7 | If available, link to developer documentation/manual | **Author action before submission.** `README.md` and `docs/methodology.md` ship with the package and are the source for the published pages. |
| C8 | Support email for questions | **Author action before submission.** |
