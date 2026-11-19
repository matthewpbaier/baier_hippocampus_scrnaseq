from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import anndata as ad
import h5py
import numpy as np
from cell_type_mapper.utils.anndata_utils import pivot_sparse_h5ad


# Allen Institute WMB-10X 20240831 reference-file SHA-256 values.
REFERENCE_SHA256 = {
    "mouse_markers_230821.json":
        "3e2a6a7beb9ab265ff545e9f8dc993b2c0ac6ae66b5e8825c42f8166a2f51035",
    "precomputed_stats_ABC_revision_230821.h5":
        "b21ca985652fb25f9608f99005139a40757133a76fbe845ae5b175c5c26a447b",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_reference(reference_dir: Path, name: str) -> Path:
    target = reference_dir / name
    if not target.is_file():
        raise FileNotFoundError(
            f"Missing Allen MapMyCells reference: {target}\n"
            "Download the pinned reference manually; see input/README.md."
        )
    digest = sha256(target)
    if digest != REFERENCE_SHA256[name]:
        raise ValueError(f"Reference checksum mismatch for {target}: {digest}")
    return target


def read_x_encoding(path: Path) -> str:
    with h5py.File(path, "r") as h5:
        value = h5["X"].attrs["encoding-type"]
        return value.decode() if isinstance(value, bytes) else str(value)


def validate_csr_conversion(source: Path, target: Path) -> None:
    original = ad.read_h5ad(source, backed="r")
    converted = ad.read_h5ad(target, backed="r")
    try:
        assert original.shape == converted.shape
        assert original.n_obs == 198122
        assert original.obs_names.is_unique and converted.obs_names.is_unique
        assert original.var_names.is_unique and converted.var_names.is_unique
        assert np.array_equal(original.obs_names, converted.obs_names)
        assert np.array_equal(original.var_names, converted.var_names)
        cells = np.sort(np.random.default_rng(1234).choice(original.n_obs, 100, replace=False))
        old = original.X[cells, :]
        new = converted.X[cells, :]
        old = old.toarray() if hasattr(old, "toarray") else np.asarray(old)
        new = new.toarray() if hasattr(new, "toarray") else np.asarray(new)
        assert np.array_equal(old, new), "CSC/CSR cell counts differ"
    finally:
        original.file.close()
        converted.file.close()


def main(root: Path) -> None:
    if not root.is_dir():
        raise FileNotFoundError(f"Data root does not exist: {root}")
    mapping_dir = root / "04_annotation" / "MapMyCells"
    input_root = Path(os.environ.get("HCA_INPUT_ROOT", str(root / "input")))
    reference_dir = input_root / "baier_mouse" / "reference" / "MapMyCells"
    result_dir = mapping_dir / "pipeline_results"
    tmp_dir = mapping_dir / "local" / "tmp_csr"
    for folder in (result_dir, tmp_dir):
        folder.mkdir(parents=True, exist_ok=True)

    source = mapping_dir / "04_baiermouse_allcells_mapmycells_input_counts.h5ad"
    csr = mapping_dir / "04_baiermouse_allcells_mapmycells_csr_counts.h5ad"
    csr_temp = mapping_dir / "04_baiermouse_allcells_mapmycells_csr_partial.h5ad"
    markers = load_reference(reference_dir, "mouse_markers_230821.json")
    stats = load_reference(reference_dir, "precomputed_stats_ABC_revision_230821.h5")

    print("Converting MapMyCells input to CSR", flush=True)
    encoding = read_x_encoding(source)
    csr_temp.unlink(missing_ok=True)
    if encoding == "csc_matrix":
        # Exact conversion call used in the original analysis.
        pivot_sparse_h5ad(
            src_path=str(source), dst_path=str(csr_temp), tmp_dir=str(tmp_dir),
            n_processors=4, max_gb=24, compression=True, layer="X",
        )
    elif encoding == "csr_matrix":
        shutil.copyfile(source, csr_temp)
    else:
        raise ValueError(f"Expected CSC or CSR sparse count matrix, found {encoding}")
    validate_csr_conversion(source, csr_temp)
    os.replace(csr_temp, csr)

    # All substantive settings match the successful original 1.7.4 run.
    csv_path = result_dir / "04_baiermouse_allcells_mapmycells_assignments.csv"
    json_path = result_dir / "04_baiermouse_allcells_mapmycells_assignments.json"
    config = {
        "query_path": str(csr),
        "map_to_ensembl": False,
        "query_markers": {
            "serialized_lookup": str(markers),
            "collapse_markers": False,
            "log_level": "ERROR",
        },
        "verbose_stdout": True,
        "hdf5_result_path": None,
        "drop_level": "CCN20230722_SUPT",
        "verbose_csv": False,
        "obsm_key": None,
        "extended_result_dir": None,
        "log_path": None,
        "gene_mapping": None,
        "csv_result_path": str(csv_path),
        "obsm_clobber": False,
        "cloud_safe": False,
        "summary_metadata_path": None,
        "log_level": "ERROR",
        "nodes_to_drop": None,
        "flatten": False,
        "precomputed_stats": {"path": str(stats), "log_level": "ERROR"},
        "tmp_dir": None,
        "type_assignment": {
            "bootstrap_iteration": 100,
            "n_processors": 4,
            "bootstrap_factor_lookup": None,
            "bootstrap_factor": 0.5,
            "chunk_size": 10000,
            "normalization": "raw",
            "n_runners_up": 5,
            "log_level": "ERROR",
            "algorithm": "hierarchical",
            "rng_seed": 11235813,
            "min_markers": 10,
        },
        "extended_result_path": str(json_path),
        "max_gb": 100.0,
        "query_gene_id_col": None,
    }
    config_path = result_dir / "04_baiermouse_allcells_mapmycells_configuration.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

    print("Running MapMyCells CCN20230722 (see pipeline_results log)", flush=True)
    log_path = result_dir / "04_baiermouse_allcells_mapmycells_run.log"
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(
            [sys.executable, "-m", "cell_type_mapper.cli.from_specified_markers",
             "--input_json", str(config_path)],
            stdout=log, stderr=subprocess.STDOUT, check=True,
        )
    print(f"MapMyCells output: {csv_path}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python 04_mapmycells.py <HCA_DATA_ROOT>")
    main(Path(sys.argv[1]).expanduser().resolve())
