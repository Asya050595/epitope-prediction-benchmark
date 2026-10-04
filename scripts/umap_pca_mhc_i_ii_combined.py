#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combined ProtT5 + UMAP/PCA visualisation for HLA classes I and II.

The script preserves the data-processing logic of the original
``umap_mhc_i.py`` and ``umap_mhc_ii.py`` scripts:

* RM and UP allele-peptide pairs are read from the same input hierarchy;
* IEDB modes are merged into one logical tool within each HLA class;
* every unique peptide is represented by a mean ProtT5 residue embedding;
* PCA pre-reduction is applied before UMAP (50 components by default);
* UMAP and PCA coordinates depend only on peptide sequence;
* overlapping RM observations with identical coordinates are separated using
  a deterministic custom beeswarm-style packing algorithm, while UP
  coordinates remain unchanged;
* all RM observations are shown by tool, while pooled UP observations form a
  neutral log-scaled hexbin density background.

One 600-dpi PNG is produced with twelve panels arranged in four rows:

  A  HLA I UMAP: p24, pp65, PtxS1
  B  HLA II UMAP: p24, pp65, PtxS1
  C  HLA I PCA:  p24, pp65, PtxS1
  D  HLA II PCA: p24, pp65, PtxS1

The ProtT5 model is loaded once: embeddings are calculated for the union of
all peptides, while every HLA-class/antigen projection is still fitted
independently exactly as in the original per-antigen workflow.

Default output:
  PROJECT_ROOT/data/
  p24/Visualization/UMAP_PCA_HLA_I_II_combined.png

Examples:
  python3 umap_pca_mhc_i_ii_combined.py --dry-run
  python3 umap_pca_mhc_i_ii_combined.py
  python3 umap_pca_mhc_i_ii_combined.py --save-coords

Dependencies:
  pip install pandas openpyxl numpy matplotlib umap-learn torch transformers \
      sentencepiece protobuf scikit-learn
"""

from __future__ import annotations

from project_paths import DATA_ROOT

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Tuple

try:
    import numpy as np
    import pandas as pd
except ImportError as exc:
    sys.exit(
        "ERROR: install the basic dependencies:\n"
        "  pip install pandas openpyxl numpy\n"
        f"Details: {exc}"
    )

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BASE = Path(f"{DATA_ROOT}")
ANTIGENS = ("p24", "pp65", "PtxS1")
DEFAULT_OUTPUT_DIR = BASE / "p24" / "Visualization"
DEFAULT_OUTPUT_NAME = "UMAP_PCA_HLA_I_II_combined.png"

# A 320-mm canvas gives each of the three plots in a row enough horizontal
# room for axis labels, a density scale and the count inset at 10-point text.
FIGURE_WIDTH_MM = 320.0
FIGURE_HEIGHT_MM = 430.0

CLASS_CONFIG = {
    "I": {
        "tools": (
            {
                "label": "IEDB I",
                "sources": (
                    ("Matches IEDB_I_Consensus", "IEDB_I_Consensus"),
                    ("Matches IEDB_I_NetMHCpan_4.1_BA", "IEDB_I_NetMHCpan_4.1_BA"),
                    ("Matches IEDB_I_NetMHCpan_4.1_EL", "IEDB_I_NetMHCpan_4.1_EL"),
                ),
            },
            {"label": "NetCTL", "sources": (("Matches NetCTL_1.2", "NetCTL_1.2"),)},
            {"label": "NetMHC", "sources": (("Matches NetMHC_4.0", "NetMHC_4.0"),)},
            {"label": "NetMHCpan", "sources": (("Matches NetMHCpan_4.1", "NetMHCpan_4.1"),)},
        ),
        "colors": {
            "IEDB I": "#C47E8F",
            "NetCTL": "#D9A441",
            "NetMHC": "#76A78E",
            "NetMHCpan": "#7E8FC4",
        },
    },
    "II": {
        "tools": (
            {
                "label": "IEDB II",
                "sources": (
                    ("Matches IEDB_II_Consensus", "IEDB_II_Consensus"),
                    ("Matches IEDB_II_NetMHCIIpan_4.1_BA", "IEDB_II_NetMHCIIpan_4.1_BA"),
                    ("Matches IEDB_II_NetMHCIIpan_4.1_EL", "IEDB_II_NetMHCIIpan_4.1_EL"),
                ),
            },
            {"label": "NetMHCII", "sources": (("Matches NetMHC_II_2.3", "NetMHC_II_2.3"),)},
            {"label": "NetMHCIIpan", "sources": (("Matches NetMHCIIpan_4.1", "NetMHCIIpan_4.1"),)},
        ),
        "colors": {
            "IEDB II": "#C47E8F",
            "NetMHCII": "#76A78E",
            "NetMHCIIpan": "#7E8FC4",
        },
    },
}

UP_DENSITY_CMAP = LinearSegmentedColormap.from_list(
    "up_density_light_grey", ["#FFFFFF", "#9A9A9A"]
)

SHEET_PRIORITY = {
    "RM": ["Union_OR", "RM", "Recovered Matches", "Recovered_Matches", "Matches"],
    "UP": ["Union_OR", "UP", "Unrecovered Predictions", "Unrecovered_Predictions", "Predictions"],
}

AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXBZUO]+$", re.IGNORECASE)
KNOWN_HLA_GENES = [
    "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5",
    "A", "B", "C", "E", "G",
]


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

def clean_colname(value: object) -> str:
    return str(value).strip().lower().replace("_", " ").replace("-", " ")


def guess_allele_col(cols: Sequence[object]) -> Optional[object]:
    for col in cols:
        if clean_colname(col) in {"allele", "hla allele"}:
            return col
    for col in cols:
        name = clean_colname(col)
        if "allele" in name and "source" not in name:
            return col
    return None


def guess_peptide_col(cols: Sequence[object]) -> Optional[object]:
    for col in cols:
        if clean_colname(col) == "peptide":
            return col
    for col in cols:
        name = clean_colname(col)
        if "peptide" in name and "source" not in name:
            return col
    return None


def sheet_has_pair_cols(xls: pd.ExcelFile, sheet_name: str) -> bool:
    try:
        header = pd.read_excel(xls, sheet_name=sheet_name, nrows=0)
    except Exception:
        return False
    return (
        guess_allele_col(header.columns) is not None
        and guess_peptide_col(header.columns) is not None
    )


def find_pair_sheet(xls: pd.ExcelFile, kind: str) -> Optional[str]:
    for sheet in SHEET_PRIORITY[kind]:
        if sheet in xls.sheet_names and sheet_has_pair_cols(xls, sheet):
            return sheet
    for sheet in xls.sheet_names:
        if sheet_has_pair_cols(xls, sheet):
            return sheet
    return None


def normalize_allele(value: object) -> str:
    """Convert HLA-A01:01/A0101/DRB1_0101 to A*01:01/DRB1*01:01."""
    allele = str(value).strip().upper()
    if not allele or allele == "NAN":
        return ""
    allele = allele.replace("HLA-", "").replace("HLA_", "").replace("HLA", "")
    allele = allele.replace(" ", "")

    if "*" in allele and ":" in allele:
        gene, rest = allele.split("*", 1)
        return f"{gene}*{rest.replace('_', '')}"

    raw = allele.replace("*", "").replace("_", "").replace(":", "")
    for gene in sorted(KNOWN_HLA_GENES, key=len, reverse=True):
        if raw.startswith(gene):
            rest = raw[len(gene):]
            if len(rest) >= 4 and rest[:4].isdigit():
                return f"{gene}*{rest[:2]}:{rest[2:4]}"
    return allele


def normalize_peptide(value: object) -> str:
    peptide = str(value).strip().upper()
    if not peptide or peptide == "NAN":
        return ""
    peptide = re.sub(r"\s+", "", peptide)
    return re.sub(r"[^A-Z]", "", peptide)


def find_input_file(
    base: Path,
    antigen: str,
    mhc_class: str,
    kind: str,
    source_dir: str,
    suffix: str,
) -> Path:
    folder = (
        base
        / f"{antigen}"
        / f"Matches MHC {mhc_class}"
        / source_dir
    )
    exact = folder / f"{kind}_{antigen}_{suffix}.xlsx"
    if exact.exists():
        return exact
    candidates = sorted(
        path
        for path in folder.glob(f"{kind}_{antigen}_{suffix}*.xlsx")
        if not path.name.startswith("~$")
    )
    return candidates[0] if candidates else exact


def load_pair_file(
    path: Path,
    *,
    kind: str,
    antigen: str,
    mhc_class: str,
    tool_label: str,
    source_suffix: str,
    strict: bool = False,
) -> pd.DataFrame:
    columns = [
        "mhc_class", "antigen", "tool", "source", "kind", "allele",
        "peptide", "input_path", "sheet",
    ]
    if not path.exists():
        message = f"[missing] {kind} | HLA {mhc_class} | {antigen} | {tool_label} | {path}"
        if strict:
            raise FileNotFoundError(message)
        print("  [!]", message)
        return pd.DataFrame(columns=columns)

    try:
        xls = pd.ExcelFile(path)
        sheet = find_pair_sheet(xls, kind)
        if sheet is None:
            raise ValueError(
                f"no sheet with allele/peptide columns; sheets={xls.sheet_names}"
            )
        frame = pd.read_excel(xls, sheet_name=sheet, keep_default_na=False)
    except Exception as exc:
        message = f"Could not read {path}: {exc}"
        if strict:
            raise RuntimeError(message) from exc
        print("  [!]", message)
        return pd.DataFrame(columns=columns)

    allele_col = guess_allele_col(frame.columns)
    peptide_col = guess_peptide_col(frame.columns)
    if allele_col is None or peptide_col is None:
        message = (
            f"Allele/peptide columns not found: {path}, sheet={sheet}, "
            f"columns={list(frame.columns)}"
        )
        if strict:
            raise ValueError(message)
        print("  [!]", message)
        return pd.DataFrame(columns=columns)

    rows = []
    for _, row in frame.iterrows():
        allele = normalize_allele(row[allele_col])
        peptide = normalize_peptide(row[peptide_col])
        if not allele or not peptide:
            continue
        if not AA_RE.match(peptide):
            print(f"  [!] suspicious peptide skipped: {peptide!r} in {path.name}")
            continue
        rows.append(
            (
                mhc_class, antigen, tool_label, source_suffix, kind, allele,
                peptide, str(path), sheet,
            )
        )

    result = pd.DataFrame(rows, columns=columns)
    if not result.empty:
        result = result.drop_duplicates(
            subset=["mhc_class", "antigen", "tool", "kind", "allele", "peptide"]
        )
    print(
        f"  {kind:2s} | HLA {mhc_class:2s} | {antigen:5s} | "
        f"{tool_label:12s} | {path.name} | sheet={sheet} | pairs={len(result)}"
    )
    return result


def load_all_pairs(
    base: Path, strict: bool = False
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rm_frames = []
    up_frames = []
    for mhc_class, config in CLASS_CONFIG.items():
        for antigen in ANTIGENS:
            print(f"\n=== HLA {mhc_class} | antigen={antigen} ===")
            for tool in config["tools"]:
                label = tool["label"]
                for source_dir, suffix in tool["sources"]:
                    rm_path = find_input_file(
                        base, antigen, mhc_class, "RM", source_dir, suffix
                    )
                    up_path = find_input_file(
                        base, antigen, mhc_class, "UP", source_dir, suffix
                    )
                    rm_frames.append(
                        load_pair_file(
                            rm_path,
                            kind="RM",
                            antigen=antigen,
                            mhc_class=mhc_class,
                            tool_label=label,
                            source_suffix=suffix,
                            strict=strict,
                        )
                    )
                    up_frames.append(
                        load_pair_file(
                            up_path,
                            kind="UP",
                            antigen=antigen,
                            mhc_class=mhc_class,
                            tool_label=label,
                            source_suffix=suffix,
                            strict=strict,
                        )
                    )

    rm = pd.concat(rm_frames, ignore_index=True) if rm_frames else pd.DataFrame()
    up = pd.concat(up_frames, ignore_index=True) if up_frames else pd.DataFrame()
    subset = ["mhc_class", "antigen", "tool", "allele", "peptide"]
    if not rm.empty:
        rm = rm.drop_duplicates(subset=subset)
    if not up.empty:
        up = up.drop_duplicates(subset=subset)
    return rm, up


def build_display_dataframe(rm: pd.DataFrame, up: pd.DataFrame) -> pd.DataFrame:
    frames = []
    if not rm.empty:
        rm_part = rm.copy()
        rm_part["status"] = "RM"
        frames.append(rm_part)
    if not up.empty:
        up_part = up.copy()
        up_part["status"] = "UP"
        frames.append(up_part)
    if not frames:
        return pd.DataFrame(
            columns=["mhc_class", "antigen", "tool", "allele", "peptide", "status"]
        )
    result = pd.concat(frames, ignore_index=True, sort=False)
    return result.drop_duplicates(
        subset=["mhc_class", "antigen", "tool", "allele", "peptide", "status"]
    )


# ---------------------------------------------------------------------------
# ProtT5 embeddings and projections
# ---------------------------------------------------------------------------

def prot_t5_input(peptide: str) -> str:
    sequence = re.sub(r"[UZOB]", "X", peptide.upper())
    return " ".join(sequence)


def compute_prott5_embeddings(
    peptides: Sequence[str],
    *,
    model_name: str,
    batch_size: int,
    device_arg: str,
) -> Dict[str, np.ndarray]:
    try:
        import torch
        from transformers import T5EncoderModel, T5Tokenizer
    except ImportError as exc:
        sys.exit(
            "ERROR: ProtT5 requires torch/transformers/sentencepiece:\n"
            "  pip install torch transformers sentencepiece protobuf\n"
            f"Details: {exc}"
        )

    device = device_arg
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"

    unique_peptides = sorted(set(peptides))
    print(
        f"[ProtT5] model={model_name} | device={device} | "
        f"peptides={len(unique_peptides)}"
    )
    tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(model_name)
    model.eval()
    model.to(device)
    if device.startswith("cuda"):
        model.half()

    result: Dict[str, np.ndarray] = {}
    with torch.no_grad():
        for start in range(0, len(unique_peptides), batch_size):
            batch_peptides = unique_peptides[start:start + batch_size]
            encoded = tokenizer.batch_encode_plus(
                [prot_t5_input(peptide) for peptide in batch_peptides],
                add_special_tokens=True,
                padding=True,
                return_tensors="pt",
            )
            input_ids = encoded["input_ids"].to(device)
            attention_mask = encoded["attention_mask"].to(device)
            hidden = model(
                input_ids=input_ids, attention_mask=attention_mask
            ).last_hidden_state.detach().cpu()
            mask = attention_mask.detach().cpu()

            for index, peptide in enumerate(batch_peptides):
                seq_len = max(int(mask[index].sum().item()) - 1, 1)
                embedding = hidden[index, :seq_len, :].mean(dim=0).float().numpy()
                result[peptide] = embedding.astype(np.float32)

            done = min(start + batch_size, len(unique_peptides))
            print(f"[ProtT5] embedded {done}/{len(unique_peptides)}")
    return result


def run_projection(
    embeddings: Mapping[str, np.ndarray],
    *,
    n_neighbors: int,
    min_dist: float,
    spread: float,
    seed: int,
    pca_components: int,
) -> Tuple[pd.DataFrame, Dict[str, Optional[float]]]:
    try:
        import umap
    except ImportError as exc:
        sys.exit(
            "ERROR: install UMAP:\n  pip install umap-learn\n"
            f"Details: {exc}"
        )

    peptides = sorted(embeddings)
    if len(peptides) < 2:
        raise ValueError("At least two unique peptides are required for UMAP.")

    raw = np.vstack([embeddings[peptide] for peptide in peptides]).astype(np.float32)
    reduced = raw
    max_components = min(raw.shape[0], raw.shape[1])
    pc1 = pc2 = None
    pc1_var: Optional[float] = None
    pc2_var: Optional[float] = None

    if pca_components > 0 and max_components > 1:
        try:
            from sklearn.decomposition import PCA
        except ImportError as exc:
            sys.exit(
                "ERROR: install scikit-learn for PCA (or use --pca-components 0):\n"
                "  pip install scikit-learn\n"
                f"Details: {exc}"
            )
        n_components = min(pca_components, max_components - 1)
        pca = PCA(n_components=n_components, random_state=seed)
        reduced = pca.fit_transform(raw)
        explained = float(pca.explained_variance_ratio_.sum())
        print(
            f"[pca] {len(peptides)} peptides: 1024 -> {n_components} dims, "
            f"explained_variance={explained:.1%}"
        )
        pc1 = reduced[:, 0]
        pc1_var = float(pca.explained_variance_ratio_[0])
        if n_components >= 2:
            pc2 = reduced[:, 1]
            pc2_var = float(pca.explained_variance_ratio_[1])
    else:
        print(
            f"[pca] disabled (pca_components={pca_components}); "
            f"using raw {reduced.shape[1]}-dim embeddings"
        )

    if pc2 is None and max_components > 1:
        try:
            from sklearn.decomposition import PCA as PCA2
        except ImportError:
            PCA2 = None
        if PCA2 is not None:
            fallback_components = min(2, max_components - 1)
            if fallback_components >= 1:
                pca2 = PCA2(n_components=fallback_components, random_state=seed)
                pca_coords = pca2.fit_transform(raw)
                pc1 = pca_coords[:, 0]
                pc1_var = float(pca2.explained_variance_ratio_[0])
                if fallback_components >= 2:
                    pc2 = pca_coords[:, 1]
                    pc2_var = float(pca2.explained_variance_ratio_[1])

    neighbors = min(n_neighbors, max(2, len(peptides) - 1))
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=neighbors,
        min_dist=min_dist,
        spread=spread,
        metric="cosine",
        random_state=seed,
    )
    umap_coords = reducer.fit_transform(reduced)

    data = {
        "peptide": peptides,
        "UMAP1": umap_coords[:, 0],
        "UMAP2": umap_coords[:, 1],
    }
    if pc1 is not None:
        data["PC1"] = pc1
    if pc2 is not None:
        data["PC2"] = pc2
    return pd.DataFrame(data), {"pc1_var": pc1_var, "pc2_var": pc2_var}


def beeswarm_offsets(n: int, min_dist: float = 1.0) -> list[Tuple[float, float]]:
    if n <= 1:
        return [(0.0, 0.0)]
    golden_angle = math.radians(137.50776405003785)
    placed: list[Tuple[float, float]] = [(0.0, 0.0)]
    for _ in range(1, n):
        candidate_index = 0
        while True:
            candidate_index += 1
            radius = min_dist * math.sqrt(candidate_index)
            angle = candidate_index * golden_angle
            candidate = (radius * math.cos(angle), radius * math.sin(angle))
            far_enough = all(
                (candidate[0] - x) ** 2 + (candidate[1] - y) ** 2
                >= (min_dist * 0.98) ** 2
                for x, y in placed
            )
            if far_enough or candidate_index > 4000:
                placed.append(candidate)
                break
    return placed


def _round_robin_indices_by_tool(frame: pd.DataFrame, indices: Sequence[int]) -> list[int]:
    """Return a stable tool-interleaved order for one coordinate group."""
    buckets: Dict[str, list[int]] = {}
    ordered = frame.loc[list(indices)].sort_values(
        ["tool", "allele", "peptide"], kind="stable"
    )
    for row_index, row in ordered.iterrows():
        buckets.setdefault(str(row["tool"]), []).append(row_index)

    result: list[int] = []
    tool_names = sorted(buckets)
    while any(buckets.values()):
        for tool in tool_names:
            if buckets[tool]:
                result.append(buckets[tool].pop(0))
    return result


def add_beeswarm_packing(
    frame: pd.DataFrame,
    jitter: float,
    x_col: str,
    y_col: str,
) -> pd.DataFrame:
    """
    Separate RM observations with identical projection coordinates.

    UMAP/PCA coordinates are peptide-level coordinates, so the same peptide
    can produce several coincident RM observations from different tools and
    alleles.  Only these RM markers are moved.  UP observations stay at their
    original coordinates so that the hexbin density represents the unmodified
    UP distribution.

    A deterministic golden-angle beeswarm is used within every group of
    identical coordinates.  Rows are interleaved by tool before packing so
    colours from different tools are distributed around the common centre
    instead of being hidden by plotting order.
    """
    result = frame.copy()
    result["plot_x"] = result[x_col]
    result["plot_y"] = result[y_col]
    if result.empty or jitter <= 0:
        return result

    x_span = max(float(result[x_col].max() - result[x_col].min()), 1.0)
    y_span = max(float(result[y_col].max() - result[y_col].min()), 1.0)
    unit_x = jitter * x_span
    unit_y = jitter * y_span

    rm = result[result["status"] == "RM"]
    # Values originating from the same peptide are bit-identical after the
    # merge. Rounding only protects the grouping from inconsequential floating
    # representation noise without merging visibly distinct points.
    coordinate_groups = rm.groupby(
        [rm[x_col].round(10), rm[y_col].round(10)], sort=False
    ).groups

    for indices in coordinate_groups.values():
        indices = _round_robin_indices_by_tool(result, list(indices))
        if len(indices) <= 1:
            continue
        for row_index, (dx, dy) in zip(indices, beeswarm_offsets(len(indices))):
            result.at[row_index, "plot_x"] = result.at[row_index, x_col] + dx * unit_x
            result.at[row_index, "plot_y"] = result.at[row_index, y_col] + dy * unit_y
    return result


def prepare_panels(
    display_df: pd.DataFrame,
    embeddings: Mapping[str, np.ndarray],
    *,
    n_neighbors: int,
    min_dist: float,
    spread: float,
    seed: int,
    pca_components: int,
    jitter: float,
) -> Dict[Tuple[str, str, str], dict]:
    panels: Dict[Tuple[str, str, str], dict] = {}
    for mhc_class in ("I", "II"):
        for antigen in ANTIGENS:
            subset = display_df[
                (display_df["mhc_class"] == mhc_class)
                & (display_df["antigen"] == antigen)
            ].copy()
            peptides = sorted(subset["peptide"].dropna().unique())
            if len(peptides) < 2:
                print(
                    f"[skip] HLA {mhc_class} {antigen}: "
                    f"need at least 2 unique peptides, found {len(peptides)}"
                )
                panels[(mhc_class, antigen, "UMAP")] = {"df": pd.DataFrame()}
                panels[(mhc_class, antigen, "PCA")] = {"df": pd.DataFrame()}
                continue

            print(f"\n[projection] HLA {mhc_class} | {antigen} | peptides={len(peptides)}")
            coords, pca_meta = run_projection(
                {peptide: embeddings[peptide] for peptide in peptides},
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                spread=spread,
                seed=seed,
                pca_components=pca_components,
            )

            umap_df = subset.merge(coords, on="peptide", how="left")
            umap_df = add_beeswarm_packing(
                umap_df, jitter=jitter, x_col="UMAP1", y_col="UMAP2"
            )
            panels[(mhc_class, antigen, "UMAP")] = {
                "df": umap_df,
                "xlabel": "UMAP-1",
                "ylabel": "UMAP-2",
            }

            if "PC2" in coords.columns:
                pca_df = subset.merge(
                    coords[["peptide", "PC1", "PC2"]], on="peptide", how="left"
                )
                pca_df = add_beeswarm_packing(
                    pca_df, jitter=jitter, x_col="PC1", y_col="PC2"
                )
                pc1_var = pca_meta.get("pc1_var")
                pc2_var = pca_meta.get("pc2_var")
                panels[(mhc_class, antigen, "PCA")] = {
                    "df": pca_df,
                    "xlabel": f"PC1 ({pc1_var:.1%} var)" if pc1_var is not None else "PC1",
                    "ylabel": f"PC2 ({pc2_var:.1%} var)" if pc2_var is not None else "PC2",
                }
            else:
                panels[(mhc_class, antigen, "PCA")] = {"df": pd.DataFrame()}
    return panels


# ---------------------------------------------------------------------------
# Combined figure
# ---------------------------------------------------------------------------

def plot_panel(
    fig: plt.Figure,
    ax: plt.Axes,
    panel: dict,
    *,
    mhc_class: str,
    antigen: str,
    up_hexbin_gridsize: int,
):
    frame = panel.get("df", pd.DataFrame())
    ax.set_title(antigen, fontsize=10, fontweight="bold", pad=6)
    if frame.empty:
        ax.text(
            0.5, 0.5, "Insufficient data", transform=ax.transAxes,
            ha="center", va="center", fontsize=10,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        return None

    up_all = frame[frame["status"] == "UP"]
    density = None
    if not up_all.empty:
        density = ax.hexbin(
            up_all["plot_x"],
            up_all["plot_y"],
            gridsize=up_hexbin_gridsize,
            cmap=UP_DENSITY_CMAP,
            mincnt=1,
            bins="log",
            edgecolors="none",
            rasterized=True,
            zorder=1,
        )

    config = CLASS_CONFIG[mhc_class]
    tool_labels = [tool["label"] for tool in config["tools"]]
    for tool in tool_labels:
        subset = frame[(frame["tool"] == tool) & (frame["status"] == "RM")]
        if not subset.empty:
            ax.scatter(
                subset["plot_x"],
                subset["plot_y"],
                marker="o",
                s=42,
                c=config["colors"].get(tool, "#999999"),
                alpha=0.95,
                edgecolors="#222222",
                linewidths=0.45,
                rasterized=True,
                zorder=2,
            )

    ax.set_xlabel(panel.get("xlabel", ""), fontsize=10, labelpad=4)
    ax.set_ylabel(panel.get("ylabel", ""), fontsize=10, labelpad=4)
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(True, linewidth=0.45, alpha=0.30)
    ax.spines[["top", "right"]].set_visible(False)

    counts = (
        frame.groupby(["tool", "status"])
        .size()
        .unstack(fill_value=0)
        .reindex(tool_labels)
        .fillna(0)
        .astype(int)
    )
    count_lines = [
        f"{tool}: RM={int(counts.loc[tool].get('RM', 0))}, "
        f"UP={int(counts.loc[tool].get('UP', 0))}"
        for tool in counts.index
    ]
    ax.text(
        0.012,
        0.012,
        "\n".join(count_lines),
        transform=ax.transAxes,
        fontsize=7.5,
        va="bottom",
        ha="left",
        linespacing=1.12,
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": "white",
            "edgecolor": "#DDDDDD",
            "alpha": 0.84,
        },
        zorder=5,
    )
    return density


def make_combined_figure(
    panels: Mapping[Tuple[str, str, str], dict],
    out_png: Path,
    *,
    dpi: int,
    up_hexbin_gridsize: int,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "font.size": 10,
            "axes.titlesize": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 10,
        }
    )

    fig, axes = plt.subplots(
        4,
        3,
        figsize=(FIGURE_WIDTH_MM / 25.4, FIGURE_HEIGHT_MM / 25.4),
        squeeze=False,
    )
    fig.subplots_adjust(
        left=0.070,
        right=0.925,
        top=0.950,
        bottom=0.115,
        wspace=0.38,
        hspace=0.72,
    )

    row_specs = (
        ("A", "HLA I UMAP", "I", "UMAP"),
        ("B", "HLA II UMAP", "II", "UMAP"),
        ("C", "HLA I PCA", "I", "PCA"),
        ("D", "HLA II PCA", "II", "PCA"),
    )
    for row, (_, _, mhc_class, projection) in enumerate(row_specs):
        row_densities = []
        for col, antigen in enumerate(ANTIGENS):
            density = plot_panel(
                fig,
                axes[row, col],
                panels.get((mhc_class, antigen, projection), {"df": pd.DataFrame()}),
                mhc_class=mhc_class,
                antigen=antigen,
                up_hexbin_gridsize=up_hexbin_gridsize,
            )
            if density is not None:
                row_densities.append(density)

        # Use one log-scaled UP-density bar for the whole antigen row.  This
        # keeps the three data panels wide and prevents the scale label from
        # colliding with the next panel's y-axis label.
        if row_densities:
            row_max = max(
                float(np.nanmax(density.get_array()))
                for density in row_densities
                if len(density.get_array())
            )
            norm = LogNorm(vmin=1.0, vmax=max(row_max, 1.01))
            for density in row_densities:
                density.set_norm(norm)
            colorbar = fig.colorbar(
                row_densities[0],
                ax=list(axes[row, :]),
                fraction=0.014,
                pad=0.018,
                shrink=0.88,
            )
            colorbar.set_label("UP count per bin", fontsize=10, labelpad=5)
            colorbar.ax.tick_params(labelsize=9, length=2.5)

    # Add row titles only after colorbars have adjusted the axes positions.
    fig.canvas.draw()
    for row, (letter, title, _, _) in enumerate(row_specs):
        positions = [axes[row, col].get_position() for col in range(3)]
        left = min(position.x0 for position in positions)
        right = max(position.x1 for position in positions)
        top = max(position.y1 for position in positions)
        title_y = min(top + 0.024, 0.976)
        fig.text(
            0.022,
            title_y,
            letter,
            ha="left",
            va="center",
            fontsize=14,
            fontweight="bold",
        )
        fig.text(
            (left + right) / 2,
            title_y,
            title,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
        )

    legend_handles = []
    for mhc_class in ("I", "II"):
        config = CLASS_CONFIG[mhc_class]
        for tool in config["tools"]:
            label = tool["label"]
            legend_handles.append(
                Patch(facecolor=config["colors"][label], edgecolor="none", label=label)
            )
    legend_handles.extend(
        [
            Line2D(
                [0], [0], marker="o", linestyle="none", label="RM",
                markerfacecolor="#777777", markeredgecolor="#222222", markersize=7,
            ),
            Patch(facecolor="#BBBBBB", edgecolor="none", label="UP (density)"),
        ]
    )
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.022),
        ncol=5,
        frameon=False,
        fontsize=10,
        handlelength=1.3,
        columnspacing=1.5,
    )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, facecolor="white")
    plt.close(fig)
    print(f"[saved] {out_png}")


def save_coordinate_tables(
    panels: Mapping[Tuple[str, str, str], dict], out_dir: Path
) -> None:
    for (mhc_class, antigen, projection), panel in panels.items():
        frame = panel.get("df", pd.DataFrame())
        if frame.empty:
            continue
        out_csv = out_dir / f"{projection.lower()}_mhc_{mhc_class.lower()}_{antigen}_coords.csv"
        frame.to_csv(out_csv, index=False)
        print(f"[saved] {out_csv}")


def print_summary(
    rm: pd.DataFrame, up: pd.DataFrame, display_df: pd.DataFrame, out_png: Path
) -> None:
    print("\n=== SUMMARY ===")
    print(f"Raw RM pairs after per-tool union: {len(rm)}")
    print(f"Raw UP pairs after per-tool union: {len(up)}")
    print(f"Displayed observations: {len(display_df)}")
    print(f"Planned PNG: {out_png}")
    if not display_df.empty:
        summary = (
            display_df.groupby(["mhc_class", "antigen", "tool", "status"])
            .size()
            .rename("n")
            .reset_index()
            .sort_values(["mhc_class", "antigen", "tool", "status"])
        )
        print(summary.to_string(index=False))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combined ProtT5 UMAP/PCA figure for HLA I and HLA II"
    )
    parser.add_argument(
        "--base", default=str(BASE), help="Base directory containing '...'"
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory; default: BASE/p24/Visualization",
    )
    parser.add_argument("--output-name", default=DEFAULT_OUTPUT_NAME)
    parser.add_argument("--strict", action="store_true", help="Fail on missing/unreadable input files")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Read inputs and print counts without computing embeddings/projections",
    )
    parser.add_argument("--save-coords", action="store_true", help="Save coordinate CSV files")
    parser.add_argument("--model-name", default="Rostlab/prot_t5_xl_uniref50")
    parser.add_argument("--device", default="auto", help="auto, cpu, cuda, cuda:0, ...")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--pca-components", type=int, default=50)
    parser.add_argument("--min-dist", type=float, default=0.65)
    parser.add_argument("--spread", type=float, default=2.5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--jitter",
        type=float,
        default=0.032,
        help=(
            "Minimum beeswarm spacing for coincident RM markers, expressed as "
            "a fraction of the projection span; default: 0.032"
        ),
    )
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--up-hexbin-gridsize", type=int, default=60)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base = Path(args.base)
    out_dir = Path(args.output_dir) if args.output_dir else (
        base / "p24" / "Visualization"
    )
    out_png = out_dir / args.output_name

    rm, up = load_all_pairs(base, strict=args.strict)
    display_df = build_display_dataframe(rm, up)
    print_summary(rm, up, display_df, out_png)

    if args.dry_run:
        print("\n[dry-run] embeddings, projections and figure were not generated")
        return
    if display_df.empty:
        sys.exit("No data to plot.")

    peptides = sorted(display_df["peptide"].dropna().unique())
    embeddings = compute_prott5_embeddings(
        peptides,
        model_name=args.model_name,
        batch_size=args.batch_size,
        device_arg=args.device,
    )
    panels = prepare_panels(
        display_df,
        embeddings,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        spread=args.spread,
        seed=args.seed,
        pca_components=args.pca_components,
        jitter=args.jitter,
    )
    if args.save_coords:
        out_dir.mkdir(parents=True, exist_ok=True)
        save_coordinate_tables(panels, out_dir)
    make_combined_figure(
        panels,
        out_png,
        dpi=args.dpi,
        up_hexbin_gridsize=args.up_hexbin_gridsize,
    )


if __name__ == "__main__":
    main()
