# Epitope prediction benchmark

Reproducible analysis scripts for benchmarking HLA class I and HLA class II
epitope-prediction workflows. The repository covers recovered/unrecovered pair
generation, threshold comparisons, recovery benchmarking, FRANK analysis,
McNemar tests, similarity analyses, and publication figures.

## Repository layout

```text
epitope-prediction-benchmark/
├── data/                 # Versioned study inputs, intermediate outputs, and figures
│   ├── p24/
│   ├── pp65/
│   └── PtxS1/
├── results/              # Optional location for collected final figures
├── scripts/              # Analysis and visualization scripts
├── requirements.txt
└── README.md
```

The scripts use paths relative to the repository and therefore do not depend on
the original author's computer. By default, `scripts/project_paths.py` resolves
the input root as `data/`. To keep large data elsewhere, set `RM_ROOT`:

```bash
export RM_ROOT=/absolute/path/to/benchmark-data
```

The directory supplied through `RM_ROOT` must contain the `p24`, `pp65`, and
`PtxS1` folders described in [data/README.md](data/README.md).

## Installation

Python 3.10 or newer is recommended.

```bash
git clone https://github.com/Asya050595/epitope-prediction-benchmark.git
cd epitope-prediction-benchmark
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The UMAP/PCA workflow additionally downloads the ProtT5 model from Hugging Face
on first use. This step requires substantial disk space and memory; GPU support
is optional.

## Input data

The repository contains the experimental reference files, prediction outputs,
RM/UP workbooks, selected intermediate results, and publication figures under
`data/<antigen>/`. The original directory hierarchy is retained within the
portable `p24`, `pp65`, and `PtxS1` directories described in
[data/README.md](data/README.md). Generated workbooks are written alongside the
corresponding antigen data, matching the validated analysis workflow.

`data/SHA256SUMS.txt` records a checksum for every versioned data file. The
repository does not distribute prediction-tool executables, model weights, or
local software caches. Benchmarking scripts stop with an error when a required
input is missing; they never report missing inputs as zero-valued results.

## Workflow overview

Run commands from the repository root. The scripts directory is added to the
module search path automatically when a script is executed directly.

### 1. Generate RM and UP files

Run the relevant `*_search_for_matches.py` scripts after placing the raw
reference datasets and prediction-tool outputs under `data/`.

```bash
python scripts/IEDB_I_Consensus_search_for_matches.py
python scripts/NetMHCpan_4_1_search_for_matches.py
python scripts/IEDB_II_Consensus_search_for_matches.py
python scripts/NetMHCIIpan_4_1_search_for_matches.py
```

Equivalent scripts are provided for the remaining IEDB modes, NetMHC, NetMHCII,
and NetCTL. Each script writes recovered matches (`RM`) and unrecovered pairs
(`UP`) to the corresponding `Matches MHC I` or `Matches MHC II` directory.

### 2. Compare prediction thresholds

```bash
python scripts/iedb_i_match_comparison.py
python scripts/iedb_ii_match_comparison.py
python scripts/netmhcpan_i_match_comparison.py
python scripts/netmhcpan_ii_match_comparison.py
```

### 3. Benchmark recovery

```bash
python scripts/benchmarking_mhc_i.py
python scripts/benchmarking_mhc_ii.py
python scripts/benchmarking_mhc_i_per_allele.py
python scripts/benchmarking_mhc_ii_per_allele.py
python scripts/micro_macro_recovery.py
```

`allele_sets.py` defines the evaluable and full reference allele sets used by
the benchmarking and statistical workflows.

### 4. Statistical and rank-based analyses

```bash
python scripts/McNemar_mhc_i_sensitivity.py
python scripts/McNemar_mhc_ii_sensitivity.py
python scripts/frank_mhc_i.py
python scripts/frank_mhc_ii.py
```

### 5. Similarity analyses

```bash
python scripts/jaccard_mhc_i.py
python scripts/jaccard_mhc_ii.py
python scripts/overlap_mhc_i.py
python scripts/overlap_mhc_ii.py
python scripts/plot_similarity_delta_combined.py
```

### 6. Publication figures

```bash
python scripts/pie_chart_RM.py
python scripts/donut_chart_RM.py
python scripts/ring_chart_RM.py
python scripts/upset_plot_RM.py
python scripts/iedb_netmhcpan_combined_visualisation.py
python scripts/frank_mcnemar_combined.py
python scripts/umap_pca_mhc_i_ii_combined.py
```

The combined figures are saved at 600 dpi. Most final images are written to
`data/p24/Visualization/` to preserve compatibility with the analysis pipeline.

## Reproducibility notes

- Analyses operate on unique HLA allele–peptide pairs.
- RM denotes recovered matches; UP denotes unrecovered pairs.
- Benchmarking scripts distinguish the common evaluable allele set from the
  full reference set.
- IEDB prediction modes are combined by set union where the script defines IEDB
  as one logical workflow; duplicates are removed before calculation.
- Statistical multiplicity correction and plotting settings are defined inside
  the corresponding scripts.
- The UMAP workflow uses a fixed random seed where specified in the code.

## Verification

Before committing changes, run:

```bash
python -m compileall -q scripts
python scripts/check_repository.py
```

The repository check rejects absolute paths tied to the original workstation,
known path-conversion corruption, and Cyrillic comments or module docstrings.
Runtime labels inherited from the validated analysis outputs are not rewritten
mechanically because some are Excel sheet names or keys used downstream.

## Citation

See [CITATION.cff](CITATION.cff). Add the publication title, complete author
list, and ORCID identifiers before the final release, then archive a tagged
release in Zenodo to obtain a citable DOI.

## License

This repository is distributed under the MIT License. Prediction-tool software,
model weights, and input datasets may have separate licenses or terms of use.
