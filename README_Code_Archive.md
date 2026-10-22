# Cross-cohort matrix-stiffness signature analysis

This archive contains the frozen analysis code and supporting files for the manuscript:

**Cross-cohort assessment of an in-vitro matrix-stiffness-associated fibroblast transcriptional program in idiopathic pulmonary fibrosis**

## Contents

- `IPF_matrix_stiffness_combat_sensitivity.py`: reproducible exploratory ComBat sensitivity analysis.
- `IPF_matrix_stiffness_supporting_data_审计更新版.xlsx`: supporting numerical results, sample manifests, and audit notes.
- `ComBat_sensitivity_audit.json`: software versions, input hashes, analysis scope, and results.
- `ComBat_sensitivity_score_comparison.csv`: pre/post-ComBat contrast comparison.
- `ComBat_sensitivity_corrected_common_signature_matrix.csv`: corrected 31-gene common matrix.
- `ComBat_sensitivity_sample_manifest.csv`: analysis-level sample manifest.
- `GSE134692_library_to_donor_manifest.csv`: library-to-donor mapping.
- `GSE124685_sample_manifest.csv`: GEO sample manifest used for the regional analysis.
- `requirements.txt`: Python package versions used for the sensitivity analysis.

## Data provenance

All expression inputs were obtained from the public NCBI GEO accessions GSE150910, GSE134692, and GSE124685. The archive does not redistribute the original GEO expression files. SHA-256 hashes of the input files are recorded in `ComBat_sensitivity_audit.json`.

## Interpretation

The ComBat analysis is exploratory sensitivity analysis rather than independent validation. The direction of four predefined case–control contrasts was preserved, while effect magnitudes changed. GSE47460 was not independently reanalyzed and contributes no present-study result.
