#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Portable implementation of upset_plot_RM. See the repository README and data/README.md for inputs, outputs, and execution instructions."""

from project_paths import DATA_ROOT

import os
import re
import unicodedata

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7,
})


# ════════════════════════════════════════════════════════════════════
# Configuration.
# ════════════════════════════════════════════════════════════════════

DEFAULT_ROOT = rf"{DATA_ROOT}"
DEFAULT_OUT_DIR = os.path.join(
    DEFAULT_ROOT, "p24", "Visualization"
)
BASE_DIR = os.environ.get("RM_ROOT", DEFAULT_ROOT)
OUTPUT_DIR = os.environ.get("RM_OUT_DIR", DEFAULT_OUT_DIR)

# Implementation detail; see the repository documentation.
ANTIGENS = ["p24", "pp65", "PtxS1"]

# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
MHC_I_TOOLS = {
    "IEDB I": [
        "IEDB_I_Consensus",
        "IEDB_I_NetMHCpan_4.1_BA",
        "IEDB_I_NetMHCpan_4.1_EL",
    ],
    "NetCTL": ["NetCTL_1.2"],
    "NetMHC": ["NetMHC_4.0"],
    "NetMHCpan": ["NetMHCpan_4.1"],
}

MHC_II_TOOLS = {
    "IEDB II": [
        "IEDB_II_Consensus",
        "IEDB_II_NetMHCIIpan_4.1_BA",
        "IEDB_II_NetMHCIIpan_4.1_EL",
    ],
    "NetMHCII": ["NetMHC_II_2.3"],
    "NetMHCIIpan": ["NetMHCIIpan_4.1"],
}

# Formatting.
# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
FONT_FAMILY = "sans-serif"
UNIQUE_COLOR = "#A8333D"       # Implementation detail; see the repository documentation.
SHARED_COLOR = "#A7ADB4"       # fixed Shared colour for intersections of ≥2 tools
MEMBER_COLOR = "#000000"       # Implementation detail; see the repository documentation.
NON_MEMBER_COLOR = "#D9D9D9"   # Implementation detail; see the repository documentation.
SETSIZE_COLOR = "#9A9A9A"      # Implementation detail; see the repository documentation.
DOT_SIZE_MATRIX = 14

FIG_WIDTH_MM = 180
FIG_HEIGHT_MM = 145
OUTPUT_DPI = 600
PANEL_TITLE_SIZE = 10
PANEL_LETTER_SIZE = 10

# ════════════════════════════════════════════════════════════════════


def normalize_allele(value) -> str:
    """Portable implementation of upset_plot_RM. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    s = unicodedata.normalize("NFKC", str(value)).strip()
    s = re.sub(r"[\u2010-\u2015\u2212]", "-", s)
    if s.upper().startswith("HLA-"):
        s = s[4:]
    return s


def normalize_peptide(value) -> str:
    return str(value).strip().upper()


def find_column(columns, keyword: str) -> str:
    for c in columns:
        if keyword in str(c).lower():
            return c
    raise ValueError(f"Не найдена колонка, содержащая '{keyword}', среди {list(columns)}")


def load_pairs(path: str) -> set:
    """Portable implementation of upset_plot_RM. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Не найден input-файл: {path}")

    try:
        xl = pd.ExcelFile(path)
        sheet = "Union_OR" if "Union_OR" in xl.sheet_names else xl.sheet_names[0]
        df = pd.read_excel(path, sheet_name=sheet)
    except Exception as e:
        raise RuntimeError(f"Не удалось прочитать файл {path}: {e}") from e

    if df.shape[1] == 0:
        return set()

    allele_col = find_column(df.columns, "allele")
    peptide_col = find_column(df.columns, "peptide")

    pairs = set()
    for _, row in df.iterrows():
        allele = normalize_allele(row[allele_col])
        peptide = normalize_peptide(row[peptide_col])
        if allele and peptide and allele.lower() != "nan" and peptide.lower() != "nan":
            pairs.add((allele, peptide))

    return pairs


def build_tool_set(antigen: str, mhc_class: str, elementary_tools: list) -> set:
    """Portable implementation of upset_plot_RM. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    pairs = set()
    for tool_token in elementary_tools:
        path = os.path.join(
            BASE_DIR,
            f"{antigen}",
            f"Matches MHC {mhc_class}",
            f"Matches {tool_token}",
            f"RM_{antigen}_{tool_token}.xlsx",
        )
        pairs |= load_pairs(path)

    return pairs


def compute_exclusive_intersections(tool_sets: dict):
    """Portable implementation of upset_plot_RM. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    all_pairs = set().union(*tool_sets.values())
    counts = {}
    for pair in all_pairs:
        combo = tuple(name for name in tool_sets if pair in tool_sets[name])
        counts[combo] = counts.get(combo, 0) + 1
    return sorted(counts.items(), key=lambda kv: (-kv[1], len(kv[0]), kv[0]))


def degree_color(degree: int, n_tools: int) -> str:
    """Use one consistent colour for every shared intersection."""
    if degree <= 1:
        return UNIQUE_COLOR
    return SHARED_COLOR


def draw_upset_panel(fig, panel_spec, tool_sets: dict,
                     title: str, letter: str) -> None:
    """Draw one compact UpSet panel inside a cell of the outer 3 × 2 grid."""
    rows = compute_exclusive_intersections(tool_sets)
    if not rows:
        ax = fig.add_subplot(panel_spec)
        ax.axis("off")
        ax.text(0.02, 0.90, letter, transform=ax.transAxes,
                ha="left", va="top", fontsize=PANEL_LETTER_SIZE,
                fontweight="bold")
        ax.text(0.55, 0.98, title, transform=ax.transAxes,
                ha="center", va="top", fontsize=PANEL_TITLE_SIZE,
                fontweight="bold")
        ax.text(0.55, 0.48, "No recovered pairs", transform=ax.transAxes,
                ha="center", va="center", fontsize=7, color="#666666")
        return

    tool_order = sorted(tool_sets, key=lambda t: -len(tool_sets[t]))
    n_rows = n_tools = len(tool_order)
    n_combos = len(rows)
    max_count = max(count for _, count in rows)
    max_size = max(max(len(tool_sets[t]) for t in tool_order), 1)

    # A dedicated middle column for tool names prevents them from colliding
    # with either the Set size bars or the membership matrix.
    gs = panel_spec.subgridspec(
        3, 3,
        width_ratios=[0.68, 0.72, 1.45],
        height_ratios=[0.30, 1.42, 0.72 + 0.17 * n_rows],
        wspace=0.04, hspace=0.08,
    )

    ax_title = fig.add_subplot(gs[0, :])
    ax_bars = fig.add_subplot(gs[1, 2])
    ax_matrix = fig.add_subplot(gs[2, 2], sharex=ax_bars)
    ax_setsize = fig.add_subplot(gs[2, 0])
    ax_labels = fig.add_subplot(gs[2, 1], sharey=ax_matrix)

    ax_title.axis("off")
    ax_title.text(0.02, 0.30, letter, ha="left", va="center",
                  fontsize=PANEL_LETTER_SIZE, fontweight="bold")
    ax_title.text(0.55, 0.74, title, ha="center", va="center",
                  fontsize=PANEL_TITLE_SIZE, fontweight="bold")

    for idx, (combo, count) in enumerate(rows):
        color = degree_color(len(combo), n_tools)
        ax_bars.bar(idx, count, color=color, width=0.7, zorder=2,
                    edgecolor="black", linewidth=0.35)
        # Vertical labels remain separated even when all 15 combinations of
        # four tools are present in a narrow journal-width panel.
        ax_bars.text(
            idx, count + max_count * 0.025, str(count),
            ha="center", va="bottom", rotation=90,
            fontsize=5.0, clip_on=False,
        )

        member_rows = [j for j, t in enumerate(tool_order) if t in combo]
        if len(member_rows) > 1:
            ax_matrix.plot([idx, idx], [min(member_rows), max(member_rows)],
                           color=MEMBER_COLOR, linewidth=0.8, zorder=1)
        for j, t in enumerate(tool_order):
            is_member = t in combo
            ax_matrix.scatter(idx, j, s=DOT_SIZE_MATRIX,
                               color=MEMBER_COLOR if is_member else NON_MEMBER_COLOR,
                               zorder=2, edgecolors="none")

    ax_bars.set_ylabel("Intersection size", fontsize=6.0, labelpad=1.5)
    ax_bars.set_xlim(-0.6, n_combos - 0.4)
    ax_bars.set_ylim(0, max_count * 1.30)
    ax_bars.set_xticks([])
    ax_bars.tick_params(axis="y", labelsize=5.3, width=0.5, length=2)
    for side in ("top", "right"):
        ax_bars.spines[side].set_visible(False)
    ax_bars.spines["left"].set_linewidth(0.5)
    ax_bars.spines["bottom"].set_linewidth(0.5)

    ax_matrix.set_xlim(-0.6, n_combos - 0.4)
    ax_matrix.set_ylim(-0.7, n_rows - 0.3)
    ax_matrix.invert_yaxis()
    ax_matrix.set_yticks([])
    ax_matrix.set_xticks([])
    for spine in ax_matrix.spines.values():
        spine.set_visible(False)

    sizes = [len(tool_sets[t]) for t in tool_order]
    bars_setsize = ax_setsize.barh(range(n_rows), sizes, color=SETSIZE_COLOR, height=0.55,
                                    edgecolor="black", linewidth=0.35)
    ax_setsize.set_xlim(0, max_size * 1.38)
    ax_setsize.set_ylim(-0.7, n_rows - 0.3)
    ax_setsize.invert_yaxis()
    ax_setsize.invert_xaxis()
    ax_setsize.bar_label(bars_setsize, padding=1.5, fontsize=5.2)
    ax_setsize.set_yticks([])
    ax_setsize.set_xlabel("Set size", fontsize=6.0, labelpad=1.5)
    ax_setsize.tick_params(axis="x", labelsize=5.0, width=0.5, length=2)
    for side in ("top", "left"):
        ax_setsize.spines[side].set_visible(False)
    ax_setsize.spines["bottom"].set_linewidth(0.5)
    ax_setsize.spines["right"].set_linewidth(0.5)

    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(-0.7, n_rows - 0.3)
    ax_labels.invert_yaxis()
    ax_labels.axis("off")
    for j, tool_name in enumerate(tool_order):
        ax_labels.text(0.5, j, tool_name, ha="center", va="center", fontsize=5.7)


def load_panel_sets(antigen: str, mhc_class: str, tools_config: dict) -> dict:
    """Load all logical tool sets for one antigen/HLA-class panel."""
    print(f"[{antigen} / HLA {mhc_class}]")
    tool_sets = {}
    for tool_name, elementary in tools_config.items():
        pairs = build_tool_set(antigen, mhc_class, elementary)
        tool_sets[tool_name] = pairs
        print(f"    {tool_name}: {len(pairs)} уникальных allele-peptide пар")
    return tool_sets


def make_combined_figure(output_path: str) -> None:
    """Create the final six-panel UpSet figure at exactly 180 mm width."""
    panels = []
    for mhc_class, tools_config in (("I", MHC_I_TOOLS), ("II", MHC_II_TOOLS)):
        for antigen in ANTIGENS:
            panels.append((
                load_panel_sets(antigen, mhc_class, tools_config),
                f"HLA {mhc_class} {antigen}",
            ))

    fig = plt.figure(
        figsize=(FIG_WIDTH_MM / 25.4, FIG_HEIGHT_MM / 25.4),
        facecolor="white",
    )
    outer = fig.add_gridspec(
        2, 3,
        left=0.02, right=0.995, top=0.992, bottom=0.125,
        wspace=0.14, hspace=0.15,
    )

    for idx, (tool_sets, title) in enumerate(panels):
        draw_upset_panel(
            fig, outer[idx // 3, idx % 3], tool_sets, title, chr(65 + idx)
        )

    from matplotlib.patches import Patch
    fig.legend(
        handles=[
            Patch(facecolor=UNIQUE_COLOR, edgecolor="black", linewidth=0.35,
                  label="Tool-specific"),
            Patch(facecolor=SHARED_COLOR, edgecolor="black", linewidth=0.35,
                  label="Shared (≥2 tools)"),
        ],
        loc="lower center", bbox_to_anchor=(0.5, 0.022),
        ncol=2, frameon=False, fontsize=7,
        handlelength=1.2, columnspacing=1.8,
    )

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=OUTPUT_DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved: {output_path}")


def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    make_combined_figure(os.path.join(OUTPUT_DIR, "upset_plot_all.png"))


if __name__ == "__main__":
    main()
