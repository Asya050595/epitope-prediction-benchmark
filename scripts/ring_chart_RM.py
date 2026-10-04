#!/usr/bin/env python3
"""
ring_chart_RM.py — combined ring-chart figure for MHC I and MHC II RM matches.

The script creates one publication-ready 3 × 2 figure:
  • top row: HLA I p24, pp65 and PtxS1;
  • bottom row: HLA II p24, pp65 and PtxS1;
  • one ring per (tool × HLA locus) cell;
  • centre label: total RMs for that tool/locus combination;
  • ring split: shared (≥2 tools) vs tool-specific pairs.

Directory structure expected:
    PROJECT_ROOT/data/
        {antigen}/
            Matches MHC {I,II}/
                Matches {tool_dir}/
                    RM_{antigen}_{tool_suffix}.xlsx

Output:
    PROJECT_ROOT/data/p24/Visualization/
        ring_chart_all.png
"""

from project_paths import DATA_ROOT

import os
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Wedge
import numpy as np

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7,
})


# Implementation detail; see the repository documentation.
COL_SHARED  = "#A7ADB4"   # cool grey-blue; fixed Shared colour across figures
COL_UNIQUE  = "#A8333D"   
COL_EMPTY   = "#E8E0DC"   # light grey for "no data" rings
FONT_FAMILY = "sans-serif"

FIG_WIDTH_MM = 180
FIG_HEIGHT_MM = 110
OUTPUT_DPI = 600
PANEL_TITLE_SIZE = 10
PANEL_LETTER_SIZE = 10
HLA_I_LOCUS_SLOTS = 3
HLA_II_LOCUS_SLOTS = 1

# ── Base paths ────────────────────────────────────────────────────────────────
DEFAULT_ROOT = f"{DATA_ROOT}"
DEFAULT_OUT_DIR = os.path.join(
    DEFAULT_ROOT, "p24", "Visualization"
)
ROOT = Path(os.environ.get("RM_ROOT", DEFAULT_ROOT))
OUT_DIR = Path(os.environ.get("RM_OUT_DIR", DEFAULT_OUT_DIR))


# ── Tool definitions ─────────────────────────────────────────────────────────
# Each entry:  display_name, mhc_class ("I" or "II"), list of (dir_fragment, file_suffix)
# IEDB modes are merged into one logical tool named "IEDB I" / "IEDB II".

MHC_I_TOOLS = [
    {
        "label":   "IEDB I",
        "sources": [
            ("Matches IEDB_I_Consensus",            "IEDB_I_Consensus"),
            ("Matches IEDB_I_NetMHCpan_4.1_BA",     "IEDB_I_NetMHCpan_4.1_BA"),
            ("Matches IEDB_I_NetMHCpan_4.1_EL",     "IEDB_I_NetMHCpan_4.1_EL"),
        ],
    },
    {
        "label":   "NetCTL",
        "sources": [
            ("Matches NetCTL_1.2",                  "NetCTL_1.2"),
        ],
    },
    {
        "label":   "NetMHC",
        "sources": [
            ("Matches NetMHC_4.0",                  "NetMHC_4.0"),
        ],
    },
    {
        "label":   "NetMHCpan",
        "sources": [
            ("Matches NetMHCpan_4.1",               "NetMHCpan_4.1"),
        ],
    },
]

MHC_II_TOOLS = [
    {
        "label":   "IEDB II",
        "sources": [
            ("Matches IEDB_II_Consensus",               "IEDB_II_Consensus"),
            ("Matches IEDB_II_NetMHCIIpan_4.1_BA",      "IEDB_II_NetMHCIIpan_4.1_BA"),
            ("Matches IEDB_II_NetMHCIIpan_4.1_EL",      "IEDB_II_NetMHCIIpan_4.1_EL"),
        ],
    },
    {
        "label":   "NetMHCII",
        "sources": [
            ("Matches NetMHC_II_2.3",                   "NetMHC_II_2.3"),
        ],
    },
    {
        "label":   "NetMHCIIpan",
        "sources": [
            ("Matches NetMHCIIpan_4.1",                 "NetMHCIIpan_4.1"),
        ],
    },
]

# HLA class definitions: prefix → display row label
MHC_I_CLASSES  = {"A": "HLA-A",  "B": "HLA-B",  "C": "HLA-C"}
MHC_II_CLASSES = {
    "DRB": "HLA-DR",   # DRB1, DRB3, DRB4, DRB5
    "DQB": "HLA-DQ",
    "DPB": "HLA-DP",
}


# ── Helpers ──────────────────────────────────────────────────────────────────

def hla_class_mhcI(allele: str) -> str | None:
    """Return 'A', 'B', 'C', etc. for a MHC-I allele string."""
    a = str(allele).strip()
    # strip optional HLA- prefix
    if a.upper().startswith("HLA-"):
        a = a[4:]
    gene = a.split("*")[0].upper()
    return gene if gene in MHC_I_CLASSES else None


def hla_class_mhcII(allele: str) -> str | None:
    """Return 'DRB', 'DQB', or 'DPB' for a MHC-II allele string."""
    a = str(allele).strip()
    if a.upper().startswith("HLA-"):
        a = a[4:]
    gene = a.split("*")[0].upper()
    for prefix in MHC_II_CLASSES:
        if gene.startswith(prefix):
            return prefix
    return None


def find_rm_file(antigen: str, mhc_class: str, dir_fragment: str, file_suffix: str) -> Path:
    """Resolve an RM input using the same layout as pie_chart_RM.py."""
    path = (
        ROOT
        / f"Processing details{antigen}"
        / f"Matches MHC {mhc_class}"
        / dir_fragment
        / f"RM_{antigen}_{file_suffix}.xlsx"
    )

    if not path.is_file():
        raise FileNotFoundError(f"Required input or value was not found{path}")

    return path


def load_pairs_from_file(path: Path, class_fn) -> dict[str, set]:
    """
    Load unique (allele, peptide) pairs from an Excel RM file.

    If the file contains a Union_OR sheet, that sheet is used because it is the
    deduplicated union of all threshold-specific RM sheets. If Union_OR is not
    present, the first sheet is used, which is expected for single-threshold RM
    files with a single prediction mode, such as NetCTL.
    """
    try:
        xl = pd.ExcelFile(path)
        sheet = "Union_OR" if "Union_OR" in xl.sheet_names else xl.sheet_names[0]
        df = pd.read_excel(path, sheet_name=sheet)
    except Exception as e:
        raise RuntimeError(f"Processing details{path}: {e}") from e

    if df.empty or df.shape[1] == 0:
        return {}

    col_map = {str(c).strip().lower(): c for c in df.columns}
    allele_col = col_map.get("allele") or col_map.get("hla allele")
    peptide_col = col_map.get("peptide")

    if allele_col is None or peptide_col is None:
        raise ValueError(
            f"Required input or value was not found{path} "
            f"Processing details{sheet!r}: {list(df.columns)}"
        )

    gene_pairs: dict[str, set] = {}

    for _, row in df.iterrows():
        allele = str(row[allele_col]).strip()
        peptide = str(row[peptide_col]).strip()

        if not allele or not peptide:
            continue
        if allele.lower() == "nan" or peptide.lower() == "nan":
            continue

        gene = class_fn(allele)
        if gene is None:
            continue

        gene_pairs.setdefault(gene, set()).add((allele, peptide))

    return gene_pairs


def load_tool_pairs(antigen: str, mhc_class: str, sources: list, class_fn) -> dict[str, set]:
    """
    Load all (allele, peptide) pairs for one logical tool, possibly merging
    several elementary modes through union.

    Returns {hla_gene_prefix: set_of_(allele, peptide)_pairs}.

    Input problems are not swallowed: if a required file is missing, unreadable,
    or lacks allele/peptide columns, the script stops instead of producing a
    misleading figure from incomplete data.
    """
    gene_pairs: dict[str, set] = {}

    for dir_fragment, file_suffix in sources:
        path = find_rm_file(antigen, mhc_class, dir_fragment, file_suffix)
        file_gene_pairs = load_pairs_from_file(path, class_fn)

        for gene, pairs in file_gene_pairs.items():
            gene_pairs.setdefault(gene, set()).update(pairs)

    return gene_pairs


# ── Core computation ─────────────────────────────────────────────────────────

def compute_shared(tools_data: list[dict]) -> dict[str, dict[str, tuple[int, int]]]:
    """
    Given tools_data = [{label, gene_pairs}], compute for each (tool, gene):
        (total_count, shared_count)
    where shared = pair present in ≥2 tools' gene buckets.

    Returns dict[tool_label][gene_prefix] = (total, shared)
    """
    # collect all genes
    all_genes = set()
    for td in tools_data:
        all_genes.update(td["gene_pairs"].keys())

    result: dict[str, dict[str, tuple[int, int]]] = {}
    for gene in all_genes:
        # count how many tools have each pair for this gene
        pair_tool_count: dict[tuple, int] = {}
        for td in tools_data:
            for pair in td["gene_pairs"].get(gene, set()):
                pair_tool_count[pair] = pair_tool_count.get(pair, 0) + 1

        for td in tools_data:
            lbl = td["label"]
            pairs = td["gene_pairs"].get(gene, set())
            total  = len(pairs)
            shared = sum(1 for p in pairs if pair_tool_count[p] >= 2)
            result.setdefault(lbl, {})[gene] = (total, shared)

    return result


# ── Drawing ──────────────────────────────────────────────────────────────────

def draw_ring(ax, total: int, shared: int,
              col_shared=COL_SHARED, col_unique=COL_UNIQUE,
              col_empty=COL_EMPTY):
    """Draw a single donut ring on ax."""
    ax.set_aspect("equal")
    ax.axis("off")

    if total == 0:
        # empty grey ring
        ring = plt.Circle((0.5, 0.5), 0.42, color=col_empty, linewidth=0)
        hole = plt.Circle((0.5, 0.5), 0.26, color="white", linewidth=0)
        ax.add_patch(ring)
        ax.add_patch(hole)
        return

    unique = total - shared

    # angles
    fracs = []
    cols  = []
    if shared > 0:
        fracs.append(shared / total)
        cols.append(col_shared)
    if unique > 0:
        fracs.append(unique / total)
        cols.append(col_unique)

    start = 90.0   # start at top
    R_outer = 0.42
    R_inner = 0.26
    cx, cy  = 0.5, 0.5

    for frac, col in zip(fracs, cols):
        angle = frac * 360.0
        wedge = Wedge(
            center=(cx, cy),
            r=R_outer,
            theta1=start - angle,
            theta2=start,
            width=R_outer - R_inner,
            facecolor=col,
            edgecolor="white",
            linewidth=0.8,
        )
        ax.add_patch(wedge)
        start -= angle

    # centre label
    ax.text(cx, cy, str(total),
            ha="center", va="center",
            fontsize=6.5, fontweight="bold",
            fontfamily=FONT_FAMILY,
            color="#2B2B2B")


def collect_panel_data(antigen: str, mhc_class: str, tools: list,
                       class_map: dict, class_fn) -> dict:
    """Load and summarise all data needed for one panel."""
    tools_data = []
    for tool in tools:
        gene_pairs = load_tool_pairs(antigen, mhc_class, tool["sources"], class_fn)
        tools_data.append({"label": tool["label"], "gene_pairs": gene_pairs})

    present_genes = [
        gene for gene in class_map
        if any(gene in td["gene_pairs"] for td in tools_data)
    ]
    # Keep one placeholder row if a panel contains no pairs, so that the
    # combined 3 × 2 layout remains complete and visibly reports the absence.
    row_genes = present_genes or [next(iter(class_map))]

    return {
        "stats": compute_shared(tools_data),
        "row_genes": row_genes,
        "row_labels": [class_map[g] for g in row_genes],
        "tools": tools,
        "has_data": bool(present_genes),
        # HLA I keeps three fixed slots so HLA-A/HLA-B are aligned in all
        # three upper panels. HLA II needs only its single HLA-DR slot.
        "layout_rows": HLA_I_LOCUS_SLOTS if mhc_class == "I" else HLA_II_LOCUS_SLOTS,
    }


def draw_panel(fig, panel_spec, panel_data: dict, title: str, letter: str) -> None:
    """Draw one compact ring-chart panel inside an outer GridSpec cell."""
    row_genes = panel_data["row_genes"]
    row_labels = panel_data["row_labels"]
    tools = panel_data["tools"]
    stats = panel_data["stats"]
    n_rows = len(row_genes)
    n_cols = len(tools)

    # Upper HLA I panels use three identical locus slots so HLA-A/HLA-B stay
    # aligned across antigens. Lower HLA II panels use only their HLA-DR slot,
    # avoiding empty rows between the rings and the shared legend.
    layout_rows = max(panel_data["layout_rows"], n_rows)

    inner = panel_spec.subgridspec(
        layout_rows + 2,
        n_cols + 1,
        height_ratios=[0.30, 0.24] + [1.0] * layout_rows,
        width_ratios=[0.48] + [1.0] * n_cols,
        hspace=0.03,
        wspace=0.03,
    )

    ax_title = fig.add_subplot(inner[0, :])
    ax_title.axis("off")
    ax_title.text(
        0.02, 0.30, letter,
        ha="left", va="center",
        fontsize=PANEL_LETTER_SIZE, fontweight="bold",
    )
    ax_title.text(
        0.55, 0.74, title,
        ha="center", va="center",
        fontsize=PANEL_TITLE_SIZE, fontweight="bold",
    )

    for ci, tool in enumerate(tools):
        ax_header = fig.add_subplot(inner[1, ci + 1])
        ax_header.axis("off")
        ax_header.text(
            0.5, 0.50, tool["label"],
            ha="center", va="center",
            fontsize=5.5,
        )

    for ri, (gene, row_label) in enumerate(zip(row_genes, row_labels)):
        ax_label = fig.add_subplot(inner[ri + 2, 0])
        ax_label.axis("off")
        ax_label.text(
            0.98, 0.5, row_label,
            ha="right", va="center",
            fontsize=6.0,
        )

        for ci, tool in enumerate(tools):
            ax = fig.add_subplot(inner[ri + 2, ci + 1])
            total, shared = stats.get(tool["label"], {}).get(gene, (0, 0))
            draw_ring(ax, total, shared)

    if not panel_data["has_data"]:
        ax_title.text(
            0.55, 0.02, "No recovered pairs",
            ha="center", va="bottom", fontsize=5.5, color="#666666",
        )


def make_combined_figure(out_path: Path) -> None:
    """Create the final six-panel figure at exactly 180 mm width."""
    panel_specs = []
    for mhc_class, tools, class_map, class_fn in (
        ("I", MHC_I_TOOLS, MHC_I_CLASSES, hla_class_mhcI),
        ("II", MHC_II_TOOLS, MHC_II_CLASSES, hla_class_mhcII),
    ):
        for antigen in ("p24", "pp65", "PtxS1"):
            print(f"Loading HLA {mhc_class} {antigen}...")
            panel_specs.append((
                collect_panel_data(antigen, mhc_class, tools, class_map, class_fn),
                f"HLA {mhc_class} {antigen}",
            ))

    fig = plt.figure(
        figsize=(FIG_WIDTH_MM / 25.4, FIG_HEIGHT_MM / 25.4),
        facecolor="white",
    )
    outer = fig.add_gridspec(
        2, 3,
        left=0.018, right=0.992, top=0.99, bottom=0.115,
        # Match the relative heights of the internal grids: the upper row has
        # three locus slots, while the lower row has one. This preserves the
        # same title/header/ring spacing without reserving two empty HLA-II rows.
        height_ratios=[3.54, 1.54],
        wspace=0.13, hspace=0.10,
    )

    for idx, (panel_data, title) in enumerate(panel_specs):
        draw_panel(fig, outer[idx // 3, idx % 3], panel_data, title, chr(65 + idx))

    legend_handles = [
        mpatches.Patch(facecolor=COL_UNIQUE, edgecolor="none",
                       label="Tool-specific"),
        mpatches.Patch(facecolor=COL_SHARED, edgecolor="none",
                       label="Shared (≥2 tools)"),
        mpatches.Patch(facecolor=COL_EMPTY, edgecolor="none",
                       label="No recovered pairs"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center", bbox_to_anchor=(0.5, 0.022),
        ncol=3, frameon=False, fontsize=7,
        handlelength=1.2, columnspacing=1.6,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=OUTPUT_DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    out_dir = OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    make_combined_figure(out_dir / "ring_chart_all.png")


if __name__ == "__main__":
    main()
