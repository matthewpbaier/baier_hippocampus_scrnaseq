
# Baier et al. — Hippocampal scRNA-seq Analysis Pipeline

## Required input data

All dataset files should be placed in `<ROOT>/input` (see bottom of readme for directory structure).

### 1. baier_mouse

**Study source:** Baier et al. (present study), hippocampal scRNA-seq, 13 male C57BL/6N mouse libraries.

**GEO accession:** To be provided upon manuscript acceptance.

**Mouse matrices** (`<ROOT>/input/baier_mouse/<sequence_id>/`):

- For each `<sequence_id>`, both Cell Ranger v8.0.1 exports:
  - `filtered_feature_bc_matrix.tar`
  - `raw_feature_bc_matrix.tar`

**Sample metadata** (`<ROOT>/input/baier_mouse/metadata/`):

- `behavior.csv`: The 44-mouse parent cohort; `mouse_id`, `age_group`, `subgroup_re80`, `InitialEntriesToEighty`, `ReversalEntriesToEighty`, `CF51to61`, `PreWeight`, and `PostWeight`.
- `sample_metadata.csv`: The 13-library sample sheet with `mouse_id`, `sequence_id`, `sequence_batch`, `age_group`, `primary_group`, `cf51to61`, `pre_weight`, and `post_weight`. `sequence_id` identifies the Cell Ranger directory (`H1`, etc.).

**Allen MapMyCells reference** (`<ROOT>/input/baier_mouse/reference/MapMyCells/`):

Obtained from the [Allen WMB-10X 20240831 reference](https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/mapmycells/WMB-10X/20240831/).

- `mouse_markers_230821.json`
- `precomputed_stats_ABC_revision_230821.h5`

**Other reference files** (`<ROOT>/input/baier_mouse/reference/`):

- `Mouse.MitoCarta3.0.xls` — [MitoCarta3.0 mouse reference](https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls)
- `mouse_mm39_gene_lengths.csv` — Two columns: `gene`, `length`; one row per unique mouse gene symbol, positive mm39 gene span in bases. The pipeline can derive this table locally when mm39 `TxDb`, `GenomicFeatures`, `AnnotationDbi`, and `org.Mm.eg.db` are available. The derived file is saved in the Script 08 output directory.

### 2. zemke_human

**Study source:** Zemke et al., *Science* (2026), DOI: 10.1126/science.adt8307. Forty-donor normative human hippocampal aging.

**GEO accession:** [GSE278576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278576)

**QC-filtered cell metadata:**

- `GSE278576_hippocampus_RNA_seurat_object_filtered_cells_metadata.tsv.gz`
  - Includes the original column name `bacrode` [sic].

**Original donor-specific 10x feature-barcode matrices:**

- `GSE278576_<orig.ident>_raw_feature_bc_matrix.h5`
  - One file corresponding to each of the 40 donor-specific `<orig.ident>` values.

### 3. SEAAD_human

**Study source:** Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD), Allen Institute for Brain Science. Human hippocampal single-nucleus RNA-seq with donor-level Alzheimer's disease neuropathology and longitudinal cognitive measures.

**Source repository:** [SEA-AD public data portal](https://brain-map.org/consortia/sea-ad/our-data)

**Hippocampal expression matrix** (`<ROOT>/input/SEAAD_human/`):

- `HIP-10X-raw.h5ad` — Hippocampal 10x single-nucleus RNA-seq raw count matrix.

**Cell, library, and gene metadata:**

- `cell_metadata.csv` — Cell-level metadata and identifiers.
- `library.csv` — Sequencing-library metadata.
- `gene.csv` — Gene-level annotations corresponding to the expression matrix.

**Cell-type reference annotations:**

- `cell_to_cluster_membership.csv` — Cell-to-cluster assignments.
- `cluster_to_cluster_annotation_membership.csv` — Mappings between transcriptomic clusters and their annotated cell identities.

**Donor-level clinical and neuropathological metadata:**

- `SEAAD_donor_metadata.xlsx` — Donor characteristics and clinical covariates.
- `SEAAD_global_local_CPS.xlsx` — Global and local continuous neuropathological scores used for pathology-associated analyses.
- `SEAAD_harmonized_cognition.xlsx` — Harmonized donor-level longitudinal cognitive measures used for pathology-adjusted cognitive analyses.

### 4. ROSMAP_human

**Study source:** Mathys et al., *Nature* (2024), DOI: 10.1038/s41586-024-07606-7. Multiregion Alzheimer's disease astrocyte transcriptomic data.

**Source repository:** [Mathys et al. multiregion DEG data](https://personal.broadinstitute.org/cboix/multiregion/DEGs/AD_DEGs/aggregated_fullset.Ast_Ast.tsv.gz)

**Astrocyte dataset** (`<ROOT>/input/ROSMAP_human/`):

- `aggregated_fullset.Ast_Ast.tsv.gz`

**Other reference files** (`<ROOT>/input/ROSMAP_human/reference/`):

- `Human.MitoCarta3.0.xls` — [Human MitoCarta3.0](https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls)
- `human.mitopathways3.0.html` — [Human mitochondrial pathway annotations](https://personal.broadinstitute.org/scalvo/MitoCarta_Download/human.mitopathways3.0.html)

## Root and directory structure

`HCA_OUTPUT_ROOT` defaults to `E:/baier_scrnaseq`, and `HCA_INPUT_ROOT` defaults to `<HCA_OUTPUT_ROOT>/input`.

The expected directory structure is:

```text
<ROOT>/
├── scripts/
│   ├── 00_paths.R
│   ├── 01_...
│   └── 14_supplemental_tables.R
└── input/
    ├── INPUT_REQUIREMENTS.csv
    ├── baier_mouse/
    │   ├── metadata/
    │   │   ├── behavior.csv
    │   │   └── sample_metadata.csv
    │   ├── H1/
    │   │   ├── filtered_feature_bc_matrix
    │   │   └── raw_feature_bc_matrix
    │   ├── <other 12 sequence_id directories>/
    │   └── reference/
    │       ├── Mouse.MitoCarta3.0.xls
    │       ├── mouse_mm39_gene_lengths.csv
    │       └── MapMyCells/
    │           ├── mouse_markers_230821.json
    │           └── precomputed_stats_ABC_revision_230821.h5
    ├── zemke_human/
    │   ├── GSE278576_hippocampus_RNA_seurat_object_filtered_cells_metadata.tsv.gz
    │   └── GSE278576_<donor_id>_raw_feature_bc_matrix.h5
    ├── SEAAD_human/
    │   ├── HIP-10X-raw.h5ad
    │   ├── cell_metadata.csv
    │   ├── library.csv
    │   ├── gene.csv
    │   ├── cell_to_cluster_membership.csv
    │   ├── cluster_to_cluster_annotation_membership.csv
    │   ├── SEAAD_donor_metadata.xlsx
    │   ├── SEAAD_global_local_CPS.xlsx
    │   └── SEAAD_harmonized_cognition.xlsx
    └── ROSMAP_human/
        ├── aggregated_fullset.Ast_Ast.tsv.gz
        └── reference/
            ├── Human.MitoCarta3.0.xls
            └── human.mitopathways3.0.html
```
