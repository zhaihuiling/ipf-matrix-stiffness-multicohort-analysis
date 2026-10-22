# IPF matrix-stiffness multicohort analysis

Reproducible analysis materials for a cross-cohort assessment of an in-vitro matrix-stiffness-associated fibroblast transcriptional program in idiopathic pulmonary fibrosis (IPF).

## Scope

This repository contains the code and audit materials used to evaluate whether a matrix-stiffness-associated fibroblast transcriptional score shows reproducible associations across public transcriptomic cohorts. The score is a transcriptional similarity measure; it is not a direct measurement of tissue stiffness and does not establish causality.

## Contents

- IPF_matrix_stiffness_combat_sensitivity.py — exploratory ComBat batch-correction sensitivity analysis.
- ComBat_sensitivity_audit.json — audit summary and direction-preservation checks.
- ComBat_sensitivity_corrected_common_signature_matrix.csv — corrected common-signature matrix used in the sensitivity analysis.
- ComBat_sensitivity_score_comparison.csv — score and contrast comparison before and after correction.
- ComBat_sensitivity_sample_manifest.csv — cohort, group, and batch labels used in the analysis.
- GSE134692_library_to_donor_manifest.csv — library-to-donor aggregation manifest.
- GSE124685_sample_manifest.csv — sample-level manifest for the public GEO cohort.
- IPF_matrix_stiffness_supporting_data_审计更新版.xlsx — supporting tables and audit notes.
- IPF_证据链审计说明.md — evidence-chain audit notes.
- IPF_matrix_stiffness_code_archive_v1.0.0.zip — versioned archive of the code package.

## Data provenance

The analysis uses publicly available GEO-derived cohort summaries and manifests. No raw participant-level clinical database is included in this repository. GSE47460 was not included in the ComBat sensitivity analysis because an independent expression input suitable for reanalysis was not available.

## Reproducibility

The script requires Python and the packages listed in requirements.txt, including pycombat==0.20. Input/output paths and cohort labels should be checked before rerunning the analysis. The reported ComBat analysis is a sensitivity analysis, not an independent validation cohort.

## Citation and versioning

Please cite the associated manuscript when available and refer to the repository URL and release/commit used for any reproduction. GitHub provides the persistent code URL; a DOI is not assigned by GitHub itself.

## License

No license has been assigned yet. Please contact the corresponding author before redistribution or reuse beyond evaluation and reproduction of the reported analysis.
