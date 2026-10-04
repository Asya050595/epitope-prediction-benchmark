#!/usr/bin/env python3
"""English documentation for this module or helper is provided in the repository README and in the surrounding code."""

from project_paths import DATA_ROOT

import os
import math
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# ---------------------------------------------------------------------------
# Formatting.
# ---------------------------------------------------------------------------
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]

TITLE_FONTSIZE    = 10
TITLE_FONTWEIGHT  = "bold"
LETTER_FONTSIZE   = 10
LETTER_FONTWEIGHT = "bold"
LABEL_FONTSIZE    = 8
LEGEND_FONTSIZE   = 8

# ---------------------------------------------------------------------------
# Configuration.
# ---------------------------------------------------------------------------

DEFAULT_ROOT = f"{DATA_ROOT}"
DEFAULT_OUT_DIR = os.path.join(
    DEFAULT_ROOT, "p24", "Visualization"
)
ROOT = os.environ.get("RM_ROOT", DEFAULT_ROOT)
OUT_DIR = os.environ.get("RM_OUT_DIR", DEFAULT_OUT_DIR)

ANTIGENS = ["p24", "pp65", "PtxS1"]  # Implementation detail; see the repository documentation.

COLORS_I = {
    "IEDB":      "#FFB7B2",
    "NetCTL":    "#FFDAC1",
    "NetMHC":    "#B5EAD7",
    "NetMHCpan": "#C7CEEA",
}

COLORS_II = {
    "IEDB II":     "#FFB7B2",
    "NetMHCII":    "#B5EAD7",
    "NetMHCIIpan": "#C7CEEA",
}

WIDTH_MM  = 179   # Implementation detail; see the repository documentation.
HEIGHT_MM = 158   # Implementation detail; see the repository documentation.
DPI       = 600

# ---------------------------------------------------------------------------
# Input loading.
# ---------------------------------------------------------------------------

def guess_allele_col(cols):
    for c in cols:
        if c.strip().lower() in ("allele", "hla allele"):
            return c
    return None


def guess_peptide_col(cols):
    for c in cols:
        if c.strip().lower() in ("peptide",):
            return c
    return None


def load_unique_pairs(paths):
    pairs = set()
    for path in paths:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Required input or value was not found{path}")
        try:
            xl = pd.ExcelFile(path)
            # Threshold handling.
            # Threshold handling.
            sheet = "Union_OR" if "Union_OR" in xl.sheet_names else xl.sheet_names[0]
            df = pd.read_excel(path, sheet_name=sheet)
        except FileNotFoundError:
            raise
        except Exception as e:
            raise RuntimeError(f"Processing details{path}: {e}") from e
        cols = list(df.columns)
        a_col = guess_allele_col(cols)
        p_col = guess_peptide_col(cols)
        if a_col is None or p_col is None:
            raise ValueError(f"Required input or value was not found{path}: {cols}")
        for _, row in df.iterrows():
            a = str(row[a_col]).strip()
            p = str(row[p_col]).strip()
            if a and p and a != "nan" and p != "nan":
                pairs.add((a, p))
    return len(pairs)


def build_tool_map_I(base, antigen):
    m = os.path.join(base, "Matches MHC I")
    def rm(*parts): return os.path.join(m, *parts)
    ag = antigen
    tools = {
        "IEDB": [
            rm("Matches IEDB_I_Consensus",         f"RM_{ag}_IEDB_I_Consensus.xlsx"),
            rm("Matches IEDB_I_NetMHCpan_4.1_BA",  f"RM_{ag}_IEDB_I_NetMHCpan_4.1_BA.xlsx"),
            rm("Matches IEDB_I_NetMHCpan_4.1_EL",  f"RM_{ag}_IEDB_I_NetMHCpan_4.1_EL.xlsx"),
        ],
        "NetCTL":    [rm("Matches NetCTL_1.2",    f"RM_{ag}_NetCTL_1.2.xlsx")],
        "NetMHC":    [rm("Matches NetMHC_4.0",    f"RM_{ag}_NetMHC_4.0.xlsx")],
        "NetMHCpan": [rm("Matches NetMHCpan_4.1", f"RM_{ag}_NetMHCpan_4.1.xlsx")],
    }
    result = {}
    for label, paths in tools.items():
        cnt = load_unique_pairs(paths)
        print(f"  MHC I | {label}: {cnt}Processing details")
        result[label] = cnt
    return result


def build_tool_map_II(base, antigen):
    m = os.path.join(base, "Matches MHC II")
    def rm(*parts): return os.path.join(m, *parts)
    ag = antigen
    tools = {
        "IEDB II": [
            rm("Matches IEDB_II_Consensus",           f"RM_{ag}_IEDB_II_Consensus.xlsx"),
            rm("Matches IEDB_II_NetMHCIIpan_4.1_BA",  f"RM_{ag}_IEDB_II_NetMHCIIpan_4.1_BA.xlsx"),
            rm("Matches IEDB_II_NetMHCIIpan_4.1_EL",  f"RM_{ag}_IEDB_II_NetMHCIIpan_4.1_EL.xlsx"),
        ],
        "NetMHCII":    [rm("Matches NetMHC_II_2.3",   f"RM_{ag}_NetMHC_II_2.3.xlsx")],
        "NetMHCIIpan": [rm("Matches NetMHCIIpan_4.1", f"RM_{ag}_NetMHCIIpan_4.1.xlsx")],
    }
    result = {}
    for label, paths in tools.items():
        cnt = load_unique_pairs(paths)
        print(f"  MHC II | {label}: {cnt}Processing details")
        result[label] = cnt
    return result


def get_counts_for_antigen(antigen):
    base = os.path.join(ROOT, f"Processing details{antigen}")

    print(f"Processing details{antigen}...")
    try:
        counts_i = build_tool_map_I(base, antigen)
    except Exception as e:
        print(f"Error{antigen}: {e}")
        counts_i = {}

    print(f"Processing details{antigen}...")
    try:
        counts_ii = build_tool_map_II(base, antigen)
    except Exception as e:
        print(f"Error{antigen}: {e}")
        counts_ii = {}

    return counts_i, counts_ii


# ---------------------------------------------------------------------------
# Normalization.
# ---------------------------------------------------------------------------

def circular_diff(a, b):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    return (a - b + 180) % 360 - 180


def spread_clustered_angles(angles, min_gap):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    n = len(angles)
    if n <= 1:
        return list(angles)

    idx_sorted = sorted(range(n), key=lambda i: angles[i])
    sorted_angles = [angles[i] for i in idx_sorted]

    gaps = [(sorted_angles[(i + 1) % n] - sorted_angles[i]) % 360 for i in range(n)]
    seam = gaps.index(max(gaps))
    order = [idx_sorted[(seam + 1 + i) % n] for i in range(n)]
    lin_angles = [angles[i] for i in order]

    unwrapped = [lin_angles[0]]
    for a in lin_angles[1:]:
        while a < unwrapped[-1]:
            a += 360
        unwrapped.append(a)

    clusters = [[0]]
    for i in range(1, n):
        if unwrapped[i] - unwrapped[i - 1] < min_gap:
            clusters[-1].append(i)
        else:
            clusters.append([i])

    result = list(unwrapped)
    for cluster in clusters:
        if len(cluster) == 1:
            continue
        mean_angle = sum(unwrapped[i] for i in cluster) / len(cluster)
        k = len(cluster)
        start = mean_angle - min_gap * (k - 1) / 2
        for j, i in enumerate(cluster):
            result[i] = start + j * min_gap

    for i in range(1, n):
        if result[i] - result[i - 1] < min_gap - 1e-6:
            result[i] = result[i - 1] + min_gap

    final = [None] * n
    for pos, orig_i in enumerate(order):
        final[orig_i] = result[pos] % 360
    return final


def normalize_percentages(values, total, decimals=1):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if total <= 0 or not values:
        return [0.0 for _ in values]
    scale = 10 ** decimals
    raw = [v / total * 100 * scale for v in values]
    floored = [math.floor(r) for r in raw]
    remainder = round(100 * scale) - sum(floored)
    order = sorted(range(len(values)), key=lambda i: raw[i] - floored[i], reverse=True)
    for i in range(max(remainder, 0)):
        floored[order[i % len(order)]] += 1
    for i in range(max(-remainder, 0)):
        floored[order[-(i + 1)]] -= 1
    if sum(floored) != 100 * scale:
        raise ArithmeticError("Processing details")
    return [f / scale for f in floored]


# ---------------------------------------------------------------------------
# Implementation detail; see the repository documentation.
# ---------------------------------------------------------------------------

def draw_pie_ax(ax, counts, colors, title, letter, label_overrides=None):
    ax.set_aspect("equal", adjustable="box")

    labels = [k for k, v in counts.items() if v > 0]
    values = [v for v in counts.values() if v > 0]
    clrs   = [colors[k] for k, v in counts.items() if v > 0]

    if not values or sum(values) == 0:
        ax.text(0.5, 0.5, "Processing details", transform=ax.transAxes,
                ha='center', va='center', fontsize=LABEL_FONTSIZE)
        ax.axis('off')
        return

    if label_overrides is None:
        label_overrides = {}

    total = sum(values)

    wedges, _ = ax.pie(
        values,
        colors=clrs,
        startangle=90,
        wedgeprops=dict(linewidth=1.0, edgecolor="white"),
        radius=0.98,
    )

    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    ax.set_xlim(-1.30, 1.18)
    ax.set_ylim(-1.20, 1.28)
    ax.text(0.5, 1.015, title, transform=ax.transAxes,
            fontsize=TITLE_FONTSIZE, fontweight=TITLE_FONTWEIGHT,
            ha="center", va="bottom", clip_on=False, zorder=10)
    # Implementation detail; see the repository documentation.
    ax.text(-1.17, 1.06, letter,
            fontsize=LETTER_FONTSIZE, fontweight=LETTER_FONTWEIGHT,
            ha="right", va="top", clip_on=False)

    pct_list = normalize_percentages(values, total, decimals=1)

    SMALL_PCT  = 12.0  # Implementation detail; see the repository documentation.
    R_INSIDE   = 0.62
    R_LINE_END = 1.08
    R_STEP     = 0.10
    MIN_GAP    = 24.0

    raw_items = [((wedge.theta1 + wedge.theta2) / 2, wedge, val, pct, lbl)
                 for wedge, val, pct, lbl in zip(wedges, values, pct_list, labels)]

    small_idx = [i for i, (_, _, _, pct, _) in enumerate(raw_items) if pct < SMALL_PCT]
    raw_small_angles = [raw_items[i][0] for i in small_idx]
    spread_small_angles = spread_clustered_angles(raw_small_angles, MIN_GAP)

    label_angles = [None] * len(raw_items)
    small_rank = [None] * len(raw_items)
    for rank, (i, ang) in enumerate(zip(small_idx, spread_small_angles)):
        label_angles[i] = ang
        small_rank[i] = rank
    for i, (angle, _, _, _, _) in enumerate(raw_items):
        if label_angles[i] is None:
            label_angles[i] = angle

    for (angle, wedge, val, pct, lbl), lbl_angle, rank in zip(raw_items, label_angles, small_rank):
        override = label_overrides.get(lbl, {})
        lbl_angle = override.get("angle", lbl_angle)

        rad_w = math.radians(angle)
        rad_l = math.radians(lbl_angle)
        cos_w, sin_w = math.cos(rad_w), math.sin(rad_w)
        cos_l, sin_l = math.cos(rad_l), math.sin(rad_l)
        text  = f"{pct:.1f}%\n(n={val})"

        if pct >= SMALL_PCT:
            ax.text(R_INSIDE * cos_w, R_INSIDE * sin_w, text,
                    ha="center", va="center", fontsize=LABEL_FONTSIZE)
        else:
            r_end = R_LINE_END + R_STEP * (rank % 2)
            r_end = override.get("r_end", r_end)
            x_text = r_end * cos_l
            y_text = r_end * sin_l
            # Implementation detail; see the repository documentation.
            # Implementation detail; see the repository documentation.
            if cos_l > 0.18:
                ha = "left"
                x_text += 0.025
            elif cos_l < -0.18:
                ha = "right"
                x_text -= 0.025
            else:
                ha = "center"
            ax.annotate(
                text,
                xy=(0.99 * cos_w, 0.99 * sin_w),
                xytext=(x_text, y_text),
                ha=ha, va="center", fontsize=LABEL_FONTSIZE,
                arrowprops=dict(arrowstyle="-", color="black", lw=0.7,
                                shrinkA=0, shrinkB=0),
            )

    ax.set_axis_off()


def add_row_legend(legend_ax, colors_dict, ncol=None):
    legend_ax.axis('off')
    handles = [mpatches.Patch(facecolor=c, label=l) for l, c in colors_dict.items()]
    legend_ax.legend(
        handles=handles, loc='center',
        ncol=ncol or len(handles),
        frameon=False, fontsize=LEGEND_FONTSIZE,
        handlelength=1.1, handleheight=1.1,
        columnspacing=1.7, borderaxespad=0.0,
    )


# ---------------------------------------------------------------------------
# Implementation detail; see the repository documentation.
# ---------------------------------------------------------------------------
# Validation.
# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
PANEL_LABEL_OVERRIDES = {}


# ---------------------------------------------------------------------------
# Main entry point.
# ---------------------------------------------------------------------------

def render_combined_figure(data, out_path):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    width_in  = WIDTH_MM / 25.4
    height_in = HEIGHT_MM / 25.4

    fig = plt.figure(figsize=(width_in, height_in))
    gs = GridSpec(4, 3, height_ratios=[10, 1.25, 10, 1.25],
                  left=0.035, right=0.985, top=0.955, bottom=0.025,
                  hspace=0.34, wspace=0.18, figure=fig)

    axes_row_I   = [fig.add_subplot(gs[0, c]) for c in range(3)]
    legend_ax_I  = fig.add_subplot(gs[1, :])
    axes_row_II  = [fig.add_subplot(gs[2, c]) for c in range(3)]
    legend_ax_II = fig.add_subplot(gs[3, :])

    letters = ["A", "B", "C", "D", "E", "F"]

    for c, antigen in enumerate(ANTIGENS):
        counts_i, _ = data[antigen]
        title = f"HLA I {antigen}"
        draw_pie_ax(axes_row_I[c], counts_i, COLORS_I, title, letters[c],
                    label_overrides=PANEL_LABEL_OVERRIDES.get(title))

    for c, antigen in enumerate(ANTIGENS):
        _, counts_ii = data[antigen]
        title = f"HLA II {antigen}"
        draw_pie_ax(axes_row_II[c], counts_ii, COLORS_II, title, letters[3 + c],
                    label_overrides=PANEL_LABEL_OVERRIDES.get(title))

    add_row_legend(legend_ax_I, COLORS_I, ncol=len(COLORS_I))
    add_row_legend(legend_ax_II, COLORS_II, ncol=len(COLORS_II))

    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    fig.savefig(out_path, dpi=DPI, facecolor='white')
    plt.close(fig)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    data = {}
    for ag in ANTIGENS:
        print(f"\n{'='*50}")
        print(f"Processing details{ag}")
        data[ag] = get_counts_for_antigen(ag)

    out_path = os.path.join(OUT_DIR, "pie_RM_all.png")
    render_combined_figure(data, out_path)
    print(f"Saved output{out_path}")
    print("Completed successfully")


if __name__ == "__main__":
    main()
