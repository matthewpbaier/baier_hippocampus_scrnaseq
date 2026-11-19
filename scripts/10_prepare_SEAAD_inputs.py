#!/usr/bin/env python3
"""Reconstruct SEA-AD HIP donor pseudobulk from locally acquired Allen releases
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


LINEAGES = ("Astrocyte", "Microglia", "OPC", "Oligodendrocyte", "DG neuron")
STEMS = {"DG neuron": "DG_neuron"}
CLINICAL = ["Donor ID", "Primary Study Name", "Age at Death", "Sex",
            "Years of education", "PMI", "RIN", "APOE Genotype", "Cognitive Status"]
SLOPES = ["slope_zmem0", "slope_zexf0", "slope_zlan0", "slope_zvsp0",
          "slope_zcasi_irt0"]
CPS_COLS = ["CPS_Global", "CPS_Global_ABeta", "CPS_Global_pTau",
            "CPS_Local", "CPS_Local_ABeta", "CPS_Local_pTau"]


def existing_table(folder, stem):
    for suffix in (".csv", ".csv.gz", ".parquet"):
        path = folder / (stem + suffix)
        if path.is_file():
            return path
    raise FileNotFoundError(f"Missing Allen metadata table {stem}.csv or .parquet in {folder}")


def read_table(folder, stem, columns=None):
    path = existing_table(folder, stem)
    if path.suffix == ".parquet":
        return pd.read_parquet(path, columns=columns)
    return pd.read_csv(path, usecols=columns, low_memory=False)


def stream_table(folder, stem, columns, batch_size=300000):
    path = existing_table(folder, stem)
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq
        for batch in pq.ParquetFile(path).iter_batches(batch_size=batch_size, columns=columns):
            yield batch.to_pandas()
    else:
        yield from pd.read_csv(path, usecols=columns, chunksize=batch_size, dtype=str)


def require_unique(frame, columns, description):
    if frame.duplicated(columns).any():
        raise ValueError(f"Nonunique {description}: {columns}")


def prepare_cohort(source, out):
    library = read_table(source, "library")
    hip_lib = library.loc[library["Brain Region"].eq("HIP"),
                          ["library_label", "donor_label"]].copy()
    hip_lib = hip_lib.astype(str)
    require_unique(hip_lib, ["library_label"], "HIP library")
    hip_libraries = set(hip_lib.library_label)

    cps = pd.read_excel(source / "SEAAD_global_local_CPS.xlsx", sheet_name="Travaglini 2026")
    cps = cps.loc[cps["Brain Region"].eq("HIP"), ["Donor ID"] + CPS_COLS]
    cognitive = pd.read_excel(source / "SEAAD_harmonized_cognition.xlsx", sheet_name="in")
    cognitive = cognitive[["Donor ID"] + SLOPES].drop_duplicates()
    clinical = pd.read_excel(source / "SEAAD_donor_metadata.xlsx",
                             sheet_name="SEA-AD_Cohort_Metadata")[CLINICAL]
    for frame, description in ((cps, "HIP CPS"), (cognitive, "cognition"),
                                (clinical, "clinical")):
        frame["Donor ID"] = frame["Donor ID"].astype(str)
        require_unique(frame, ["Donor ID"], description + " donor")
    hip_donors = pd.DataFrame({"Donor ID": sorted(hip_lib.donor_label.unique())})
    cohort = (hip_donors.merge(cps, on="Donor ID", how="left", validate="one_to_one")
              .merge(cognitive, on="Donor ID", how="left", validate="one_to_one")
              .merge(clinical, on="Donor ID", how="left", validate="one_to_one"))
    cohort["primary_slope_available"] = cohort["slope_zcasi_irt0"].notna()
    cognitive_donors = set(cohort.loc[cohort.primary_slope_available, "Donor ID"])

    hip_chunks = []
    for chunk in stream_table(source, "cell_metadata", ["cell_label", "library_label", "donor_label"]):
        selected = chunk[chunk.library_label.isin(hip_libraries) &
                         chunk.donor_label.isin(cognitive_donors)]
        if not selected.empty:
            hip_chunks.append(selected.copy())
    hip_cells = pd.concat(hip_chunks, ignore_index=True)
    require_unique(hip_cells, ["cell_label"], "HIP cell")
    print(f"HIP libraries: {len(hip_libraries)}; cognition-complete donors: {len(cognitive_donors)}")

    hip_ids = set(hip_cells.cell_label)
    cluster_parts = []
    for chunk in stream_table(source, "cell_to_cluster_membership", ["cell_label", "cluster_alias"]):
        selected = chunk[chunk.cell_label.isin(hip_ids)]
        if not selected.empty:
            cluster_parts.append(selected.copy())
    membership = pd.concat(cluster_parts, ignore_index=True)
    require_unique(membership, ["cell_label"], "cell-to-cluster assignment")
    annotations = read_table(source, "cluster_to_cluster_annotation_membership")
    annotations = annotations[["cluster_alias", "cluster_annotation_term_set_name",
                               "cluster_annotation_term_name"]].drop_duplicates()
    require_unique(annotations, ["cluster_alias", "cluster_annotation_term_set_name"],
                   "cluster annotation level")
    wide = annotations.pivot(index="cluster_alias", columns="cluster_annotation_term_set_name",
                              values="cluster_annotation_term_name").reset_index()
    levels = list(wide.columns)
    class_name = next(x for x in levels if str(x).lower() == "class")
    subclass_name = next(x for x in levels if str(x).lower() == "subclass")
    fine = [x for x in levels if x not in ("cluster_alias", class_name, subclass_name)]
    if len(fine) != 1:
        raise ValueError(f"Expected one fine taxonomy level, found {fine}")
    wide = wide.rename(columns={class_name: "class_name", subclass_name: "subclass_name",
                                fine[0]: "fine_type_name"})
    wide["cluster_alias"] = wide.cluster_alias.astype(str).str.replace(r"\.0$", "", regex=True)
    membership["cluster_alias"] = membership.cluster_alias.astype(str).str.replace(r"\.0$", "", regex=True)
    hip = (hip_cells.merge(membership, on="cell_label", how="left", validate="one_to_one")
           .merge(wide, on="cluster_alias", how="left", validate="many_to_one"))

    sub = hip.subclass_name.fillna("")
    fine_type = hip.fine_type_name.fillna("")
    hip["crossstudy_celltype"] = None
    for source_name, target in (("Astrocyte", "Astrocyte"), ("OPC", "OPC"),
                                ("Oligodendrocyte", "Oligodendrocyte"),
                                ("DG", "DG neuron")):
        hip.loc[sub.eq(source_name), "crossstudy_celltype"] = target
    hip.loc[sub.eq("Immune") & fine_type.str.startswith("Micro-PVM"),
            "crossstudy_celltype"] = "Microglia"

    primary = hip[hip.crossstudy_celltype.isin(LINEAGES)].copy()
    eligible = sorted(primary.donor_label.unique())
    analysis_cohort = cohort[cohort["Donor ID"].isin(eligible)].sort_values("Donor ID")
    excluded = sorted(cognitive_donors.difference(eligible))
    print(f"Expression-eligible donors: {len(eligible)}; excluded: {excluded}")
    if len(analysis_cohort) != 31:
        raise ValueError(f"Expected 31 SEA-AD primary donors, found {len(analysis_cohort)}")
    if analysis_cohort[["Sex", "Age at Death", "CPS_Local", "slope_zcasi_irt0",
                        "slope_zmem0", "PMI", "RIN"]].isna().any().any():
        raise ValueError("Missing covariates in the reconstructed 31-donor cohort")
    analysis_cohort.to_csv(out / "10_seaad_alldonors_cognition_eligible_cohort.csv", index=False)
    (cohort[cohort["Donor ID"].isin(excluded)]
     .to_csv(out / "10_seaad_alldonors_cognition_excluded_donors.csv", index=False))
    primary = primary[primary.donor_label.isin(eligible)].copy()
    coverage = (primary.groupby(["donor_label", "crossstudy_celltype"]).size()
                .rename("n_cells").reset_index())
    full = (pd.MultiIndex.from_product([eligible, LINEAGES],
               names=["donor_label", "crossstudy_celltype"]).to_frame(index=False)
            .merge(coverage, how="left", on=["donor_label", "crossstudy_celltype"]))
    full["n_cells"] = full.n_cells.fillna(0).astype(int)
    full.to_csv(out / "10_seaad_alllineages_cognition_donor_coverage.csv", index=False)
    if (full.n_cells == 0).any():
        raise ValueError("Selected donor has a zero-cell primary lineage")
    primary.to_csv(out / "10_seaad_alllineages_cognition_selected_cell_annotations.csv.gz", index=False,
                   compression="gzip")
    print(f"Selected annotated nuclei: {len(primary):,}")
    return primary, eligible


def aggregate_raw(source, out, selected, donors):
    import anndata as ad
    import scipy.sparse as sp

    path = source / "HIP-10X-raw.h5ad"
    if not path.is_file():
        raise FileNotFoundError(f"Download the Allen HIP-10X/raw H5AD and save it as {path}")
    adata = ad.read_h5ad(path, backed="r")
    try:
        obs_id = pd.Index(adata.obs_names.astype(str))
        selected_ids = pd.Index(selected.cell_label.astype(str))
        if not selected_ids.isin(obs_id).all():
            raise ValueError(f"{(~selected_ids.isin(obs_id)).sum()} selected cells missing from H5AD index")
        if not obs_id.is_unique:
            raise ValueError("Duplicate H5AD cell IDs")
        position = pd.Series(np.arange(adata.n_obs, dtype=np.int64), index=obs_id)
        selected["expression_row"] = position.loc[selected_ids].to_numpy()
        gene_table = read_table(source, "gene")
        var_names = pd.Index(adata.var_names.astype(str))
        gene_ids = var_names.to_numpy()
        symbol_column = next((x for x in ("gene_symbol", "gene_name", "gene", "symbol")
                              if x in adata.var.columns), None)
        if symbol_column is not None:
            genes = adata.var[symbol_column].astype(str).to_numpy()
        else:
            gene_lookup = gene_table.drop_duplicates("gene_identifier").set_index("gene_identifier")
            genes = gene_lookup.reindex(var_names).gene_symbol
            if genes.isna().any():
                raise ValueError("H5AD genes not fully represented in Allen gene metadata")
            genes = genes.astype(str).to_numpy()
        pd.DataFrame({"gene_identifier": gene_ids, "gene_symbol": genes}).to_csv(
            out / "10_seaad_alllineages_cognition_raw_gene_map.csv", index=False)

        keys = pd.MultiIndex.from_product([LINEAGES, donors],
                  names=["cell_type", "donor_label"]).to_frame(index=False)
        keys["group_code"] = np.arange(len(keys))
        code_map = {(r.donor_label, r.cell_type): r.group_code
                    for r in keys.itertuples(index=False)}
        selected["group_code"] = [code_map[(d, t)] for d, t in
                                   zip(selected.donor_label, selected.crossstudy_celltype)]
        row_group = np.full(adata.n_obs, -1, dtype=np.int32)
        row_group[selected.expression_row.to_numpy(dtype=np.int64)] = selected.group_code.to_numpy()
        counts = np.zeros((len(keys), adata.n_vars), dtype=np.int64)
        for start in range(0, adata.n_obs, 4096):
            stop = min(start + 4096, adata.n_obs)
            valid = np.flatnonzero(row_group[start:stop] >= 0)
            if not len(valid):
                continue
            matrix = adata.X[start:stop, :][valid, :]
            for code in np.unique(row_group[start:stop][valid]):
                part = matrix[row_group[start:stop][valid] == code, :]
                summed = np.asarray(part.sum(axis=0)).ravel()
                if not np.allclose(summed, np.rint(summed), atol=1e-6):
                    raise ValueError("Noninteger count values in HIP raw expression")
                counts[int(code), :] += np.rint(summed).astype(np.int64)
            if start % (4096 * 20) == 0:
                print(f"  raw expression rows: {stop:,} / {adata.n_obs:,}")
        for ct in LINEAGES:
            idx = keys.index[keys.cell_type == ct].to_numpy()
            donor_names = keys.loc[idx, "donor_label"].tolist()
            frame = pd.DataFrame(counts[idx, :].T, columns=donor_names)
            frame.insert(0, "gene_symbol", genes)
            frame.insert(0, "gene_identifier", gene_ids)
            out_name = f"10_seaad_{STEMS.get(ct, ct).lower()}_cognition_raw_pseudobulk_counts.csv.gz"
            frame.to_csv(out / out_name, index=False, compression="gzip")
            print(f"Saved {out_name}: {len(genes):,} genes x {len(donor_names)} donors")
        keys["n_cells"] = selected.groupby(["crossstudy_celltype", "donor_label"]).size().reindex(
            pd.MultiIndex.from_frame(keys[["cell_type", "donor_label"]]), fill_value=0).to_numpy()
        keys["library_size_UMI"] = counts.sum(axis=1)
        keys.to_csv(out / "10_seaad_alllineages_cognition_pseudobulk_qc.csv", index=False)
    finally:
        adata.file.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path,
                        help="Manually acquired original SEA-AD data directory")
    parser.add_argument("--output", required=True, type=Path,
                        help="New derived files directory, outside input/")
    parser.add_argument("--metadata-only", action="store_true",
                        help="Reconstruct and audit cohort without opening the large H5AD")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source = args.input
    for name in ("SEAAD_global_local_CPS.xlsx", "SEAAD_harmonized_cognition.xlsx",
                 "SEAAD_donor_metadata.xlsx"):
        if not (source / name).is_file():
            raise FileNotFoundError(source / name)
    selected, donors = prepare_cohort(source, args.output)
    if not args.metadata_only:
        aggregate_raw(source, args.output, selected, donors)


if __name__ == "__main__":
    main()
