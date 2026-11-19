# baier_hippocampus_scrnaseq
scRNA-seq analysis pipeline for Baier et al., 2026

#Baier et al.
#Hippocampal scRNA-seq pipeline, v1.0


## Required input data
All dataset files should be placed in "<ROOT>/input".


#1. "baier_mouse"
*Study source: Baier et al. (present study), hippocampal scRNA-seq, 13 male C57BL/6N mouse libraries. 
*GEO accession: <to be provided upon manuscript acceptance>

*Mouse matrices (<root>/baier_mouse/<sequence_id>)
  -For each <sequence_id>, both Cell Ranger v8.0.1 directory exports "filtered_feature_bc_matrix.tar" and "raw_feature_bc_matrix.tar"

*Sample metadata (<root>/baier_mouse/metadata): 
  -"/behavior.csv": the 44-mouse parent cohort cohort; `mouse_id`, `age_group`, `subgroup_re80`, `InitialEntriesToEighty`, `ReversalEntriesToEighty`, `CF51to61`, `PreWeight`, and `PostWeight`.
  -"/sample_metadata.csv": 13-library sample sheet with `mouse_id`, `sequence_id`, `sequence_batch`, `age_group`, `primary_group`, `cf51to61`, `pre_weight`, and `post_weight`. `sequence_id` identifies the Cell Ranger directory (`H1`, etc.).

*Allen MapMyCells reference (<root>/baier_mouse/reference/mapmycells):
  -Obtained from Allen WMB-10X 20240831 reference directory (https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/mapmycells/WMB-10X/20240831/)
  
  -"/mouse_markers_230821.json"
  -"/precomputed_stats_ABC_revision_230821.h5"

*Other reference files (<root>/baier_mouse/reference):
  -"/Mouse.MitoCarta3.0.xls" (https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls)
  -"/mouse_mm39_gene_lengths.csv" Two columns: `gene`, `length`; one row per unique mouse gene symbol, positive mm39 gene span in bases. The pipeline will derive this table locally when mm39 `TxDb` plus `GenomicFeatures`, `AnnotationDbi`, and `org.Mm.eg.db` are available and saved in Script 08 output folder.


#2. "zemke_human"
*Study source: Zemke et al., Science 2026, DOI: 10.1126/science.adt8307, 40-donor normative human hippocampal aging.
*GEO accession: GSE278576 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278576)

*QC-filtered cell metadata (includes the original column name `bacrode` [sic]).
  -"/GSE278576_hippocampus_RNA_seurat_object_filtered_cells_metadata.tsv.gz"

*Original donor-specific 10x feature-barcode matrix.
  -"/GSE278576_<orig.ident>_raw_feature_bc_matrix.h5" (corresponding to the 40 <orig.ident> donors)


#3. "SEAAD_human"
*Study source: Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD), Allen Institute for Brain Science. Human hippocampal single-nucleus RNA-seq with donor-level Alzheimer's disease neuropathology and longitudinal cognitive measures.
*Source repository: SEA-AD public data portal (https://brain-map.org/consortia/sea-ad/our-data)

*Hippocampal expression matrix (<root>/SEAAD_human):
  -"/HIP-10X-raw.h5ad": hippocampal 10x single-nucleus RNA-seq raw count matrix.

*Cell, library, and gene metadata (<root>/SEAAD_human):
  -"/cell_metadata.csv": cell-level metadata and identifiers.
  -"/library.csv": sequencing-library metadata.
  -"/gene.csv": gene-level annotations corresponding to the expression matrix.

*Cell-type reference annotations (<root>/SEAAD_human):
  -"/cell_to_cluster_membership.csv": cell-to-cluster assignments.
  -"/cluster_to_cluster_annotation_membership.csv": mappings between transcriptomic clusters and their annotated cell identities.

*Donor-level clinical and neuropathological metadata (<root>/SEAAD_human):
  -"/SEAAD_donor_metadata.xlsx": donor characteristics and clinical covariates.
  -"/SEAAD_global_local_CPS.xlsx": global and local continuous neuropathological scores used for pathology-associated analyses.
  -"/SEAAD_harmonized_cognition.xlsx": harmonized donor-level longitudinal cognitive measures used for pathology-adjusted cognitive analyses.


#4. "ROSMAP_human"
*Study source: Mathys et al., Nature 2024, DOI: 10.1038/s41586-024-07606-7, multiregion Alzheimer's disease (astrocytes).
*Source repository: https://personal.broadinstitute.org/cboix/multiregion/DEGs/AD_DEGs/aggregated_fullset.Ast_Ast.tsv.gz

*Astrocyte dataset:
  -"/aggregated_fullset.Ast_Ast.tsv.gz"

*Other reference files (<root>/ROSMAP_human/reference):
  -"/Human.MitoCarta3.0.xls" (https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls)
  -"/human.mitopathways3.0.html" (https://personal.broadinstitute.org/scalvo/MitoCarta_Download/human.mitopathways3.0.html)


##Root and directory structure
`HCA_OUTPUT_ROOT` defaults to `E:/baier_scrnaseq`, and `HCA_INPUT_ROOT` defaults to `<HCA_OUTPUT_ROOT>/input`.

The directory structure for this pipeline is as follows:
<ROOT>
  scripts/                     
  input/
    README.md
    INPUT_REQUIREMENTS.csv
    baier_mouse/
      metadata/
        behavior.csv
        sample_metadata.csv
      H1/
        filtered_feature_bc_matrix/
          matrix.mtx.gz
          barcodes.tsv.gz
          features.tsv.gz
        raw_feature_bc_matrix/
          matrix.mtx.gz
          barcodes.tsv.gz
          features.tsv.gz
      <other 12 sequence_id folders following same structure>/
      reference/
        Mouse.MitoCarta3.0.xls
        mouse_mm39_gene_lengths.csv     
        MapMyCells/
          mouse_markers_230821.json
          precomputed_stats_ABC_revision_230821.h5
    zemke_human/
      GSE278576_hippocampus_RNA_seurat_object_filtered_cells_metadata.tsv.gz
      GSE278576_<donor_id>_raw_feature_bc_matrix.h5   #<40 donor files>
    SEAAD_human/
      README.md
      HIP-10X-raw.h5ad
      cell_metadata.csv
      library.csv
      gene.csv
      cell_to_cluster_membership.csv
      cluster_to_cluster_annotation_membership.csv
      SEAAD_donor_metadata.xlsx
      SEAAD_global_local_CPS.xlsx
      SEAAD_harmonized_cognition.xlsx
    ROSMAP_human/
      aggregated_fullset.Ast_Ast.tsv.gz
      reference/
        Human.MitoCarta3.0.xls
        human.mitopathways3.0.html
