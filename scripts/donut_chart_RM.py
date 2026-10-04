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
from collections import defaultdict

# ---------------------------------------------------------------------------
# Formatting.
# ---------------------------------------------------------------------------
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]

TITLE_FONTSIZE          = 10
TITLE_FONTWEIGHT        = "bold"
LETTER_FONTSIZE         = 10
LETTER_FONTWEIGHT       = "bold"
LABEL_FONTSIZE          = 6.4  # Implementation detail; see the repository documentation.
CENTER_TOTAL_FONTSIZE   = 12.5 # Implementation detail; see the repository documentation.
CENTER_CAPTION_FONTSIZE = 4.8  # Implementation detail; see the repository documentation.
LEGEND_FONTSIZE         = 7

# ---------------------------------------------------------------------------
# Configuration.
# ---------------------------------------------------------------------------

DEFAULT_ROOT = f"{DATA_ROOT}"
DEFAULT_OUT_DIR = f"{DATA_ROOT}/p24/Visualization"

ROOT = os.environ.get("RM_ROOT", DEFAULT_ROOT)
OUT_DIR = os.environ.get("RM_OUT_DIR", DEFAULT_OUT_DIR)

ANTIGENS = ["p24", "pp65", "PtxS1"]  # Implementation detail; see the repository documentation.

MHC_I_TOOLS = {
    "IEDB": [
        ("Matches IEDB_I_Consensus",        "IEDB_I_Consensus"),
        ("Matches IEDB_I_NetMHCpan_4.1_BA", "IEDB_I_NetMHCpan_4.1_BA"),
        ("Matches IEDB_I_NetMHCpan_4.1_EL", "IEDB_I_NetMHCpan_4.1_EL"),
    ],
    "NetCTL":    [("Matches NetCTL_1.2",    "NetCTL_1.2")],
    "NetMHC":    [("Matches NetMHC_4.0",    "NetMHC_4.0")],
    "NetMHCpan": [("Matches NetMHCpan_4.1", "NetMHCpan_4.1")],
}

MHC_II_TOOLS = {
    "IEDB II": [
        ("Matches IEDB_II_Consensus",           "IEDB_II_Consensus"),
        ("Matches IEDB_II_NetMHCIIpan_4.1_BA",  "IEDB_II_NetMHCIIpan_4.1_BA"),
        ("Matches IEDB_II_NetMHCIIpan_4.1_EL",  "IEDB_II_NetMHCIIpan_4.1_EL"),
    ],
    "NetMHCII":    [("Matches NetMHC_II_2.3",    "NetMHC_II_2.3")],
    "NetMHCIIpan": [("Matches NetMHCIIpan_4.1",  "NetMHCIIpan_4.1")],
}

MHC_I_COLORS = {
    "IEDB":      "#FFB7B2",
    "NetCTL":    "#FFDAC1",
    "NetMHC":    "#B5EAD7",
    "NetMHCpan": "#C7CEEA",
    "Shared (\u22652 tools)": "#A7ADB4",
}

MHC_II_COLORS = {
    "IEDB II":     "#FFB7B2",
    "NetMHCII":    "#B5EAD7",
    "NetMHCIIpan": "#C7CEEA",
    "Shared (\u22652 tools)": "#A7ADB4",
}

# Implementation detail; see the repository documentation.
WIDTH_MM  = 179   # Implementation detail; see the repository documentation.
HEIGHT_MM = 158   # Implementation detail; see the repository documentation.
DPI       = 600

# ---------------------------------------------------------------------------
# Input loading.
# ---------------------------------------------------------------------------

def guess_allele_col(cols):
    for c in cols:
        if str(c).strip().lower() in ("allele", "hla allele"):
            return c
    return None


def guess_peptide_col(cols):
    for c in cols:
        if str(c).strip().lower() == "peptide":
            return c
    return None


def load_pairs_from_file(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Required input or value was not found{path}")

    try:
        xl = pd.ExcelFile(path)
        sheet = "Union_OR" if "Union_OR" in xl.sheet_names else xl.sheet_names[0]
        df = pd.read_excel(path, sheet_name=sheet, header=0)
    except Exception as e:
        raise RuntimeError(f"Processing details{path}: {e}") from e

    cols = list(df.columns)
    a_col = guess_allele_col(cols)
    p_col = guess_peptide_col(cols)

    if a_col is None or p_col is None:
        raise ValueError(f"Required input or value was not found{path}: {cols}")

    pairs = set()
    for _, row in df.iterrows():
        a = str(row[a_col]).strip()
        p = str(row[p_col]).strip()
        if a and p and a.lower() != "nan" and p.lower() != "nan":
            pairs.add((a, p))
    return pairs


def collect_tool_pairs(base_dir, antigen, tools_config):
    tool_pairs = {}
    for tool_name, file_list in tools_config.items():
        pairs = set()
        for folder, suffix in file_list:
            path = os.path.join(base_dir, folder, f"RM_{antigen}_{suffix}.xlsx")
            pairs |= load_pairs_from_file(path)  # Implementation detail; see the repository documentation.
        print(f"  {tool_name}: {len(pairs)}Processing details")
        tool_pairs[tool_name] = pairs
    return tool_pairs


def compute_segments(tool_pairs):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    pair_count = defaultdict(int)
    for pairs in tool_pairs.values():
        for p in pairs:
            pair_count[p] += 1

    total = len(pair_count)
    if total == 0:
        print("Warning")
        return None, 0

    shared = {p for p, c in pair_count.items() if c >= 2}
    segments = {}
    for name, pairs in tool_pairs.items():
        u = len(pairs - shared)
        if u > 0:
            segments[name] = u
    if len(shared) > 0:
        segments["Shared (\u22652 tools)"] = len(shared)

    return segments, total


def get_segments_for_antigen(antigen):
    base = os.path.join(ROOT, f"Processing details{antigen}")

    print(f"Processing details{antigen}")
    try:
        mhc_i_dir = os.path.join(base, "Matches MHC I")
        tool_pairs_I = collect_tool_pairs(mhc_i_dir, antigen, MHC_I_TOOLS)
        segments_I, total_I = compute_segments(tool_pairs_I)
    except Exception as e:
        print(f"Error{antigen}: {e}")
        segments_I, total_I = None, 0

    print(f"Processing details{antigen}")
    try:
        mhc_ii_dir = os.path.join(base, "Matches MHC II")
        tool_pairs_II = collect_tool_pairs(mhc_ii_dir, antigen, MHC_II_TOOLS)
        segments_II, total_II = compute_segments(tool_pairs_II)
    except Exception as e:
        print(f"Error{antigen}: {e}")
        segments_II, total_II = None, 0

    return segments_I, total_I, segments_II, total_II


# ---------------------------------------------------------------------------
# Normalization.
# ---------------------------------------------------------------------------

def radius_clearing_center_box(angle_deg, half_w, half_h, min_r=0.0):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    rad = math.radians(angle_deg)
    cos_a, sin_a = math.cos(rad), math.sin(rad)
    candidates = [min_r]
    if abs(cos_a) > 1e-9:
        r_x = half_w / abs(cos_a)
        if abs(r_x * sin_a) <= half_h:
            candidates.append(r_x)
    if abs(sin_a) > 1e-9:
        r_y = half_h / abs(sin_a)
        if abs(r_y * cos_a) <= half_w:
            candidates.append(r_y)
    return max(candidates)


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

    # Implementation detail; see the repository documentation.
    gaps = [(sorted_angles[(i + 1) % n] - sorted_angles[i]) % 360 for i in range(n)]
    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
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

    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
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

def draw_donut_ax(ax, segments, total, colors, title, letter, label_overrides=None):
    ax.set_aspect("equal", adjustable="box")

    if not segments or total <= 0:
        ax.text(0.5, 0.5, "Processing details", transform=ax.transAxes,
                ha='center', va='center', fontsize=LABEL_FONTSIZE)
        ax.axis('off')
        return

    if label_overrides is None:
        label_overrides = {}

    labels = list(segments.keys())
    sizes  = [segments[l] for l in labels]
    clrs   = [colors.get(l, "#AAAAAA") for l in labels]

    wedges, _ = ax.pie(
        sizes,
        colors=clrs,
        startangle=90,
        radius=1.10,
        wedgeprops=dict(width=0.54, edgecolor='white', linewidth=1.0),
    )

    ax.set_xlim(-1.38, 1.28)
    ax.set_ylim(-1.38, 1.50)
    # Implementation detail; see the repository documentation.
    ax.text(-1.24, 1.36, letter,
            fontsize=LETTER_FONTSIZE, fontweight=LETTER_FONTWEIGHT,
            ha="right", va="top", clip_on=False)

    pct_list = normalize_percentages(sizes, total, decimals=1)

    SMALL_PCT        = 11.0
    R_INSIDE         = 0.83   # Implementation detail; see the repository documentation.
    R_LINE_END_SHORT = 1.20
    R_LINE_END_LONG  = 1.30
    R_STEP           = 0.09
    TINY_PCT         = 3.0
    MIN_GAP          = 24.0

    ax.text(0, 0.075, str(total), ha='center', va='center',
            fontsize=CENTER_TOTAL_FONTSIZE, fontweight='normal')
    ax.text(0, -0.12, 'Recovered matches', ha='center', va='center',
            fontsize=CENTER_CAPTION_FONTSIZE, fontweight='normal')

    items = [(w, v, p, (w.theta1 + w.theta2) / 2, lbl)
             for w, v, p, lbl in zip(wedges, sizes, pct_list, labels)]

    small_idx = [i for i, (_, _, pct, _, _) in enumerate(items) if pct < SMALL_PCT]
    raw_small_angles = [items[i][3] for i in small_idx]
    spread_small_angles = spread_clustered_angles(raw_small_angles, MIN_GAP)

    label_angles = [None] * len(items)
    small_rank = [None] * len(items)
    for rank, (i, ang) in enumerate(zip(small_idx, spread_small_angles)):
        label_angles[i] = ang
        small_rank[i] = rank
    for i, (_, _, _, angle, _) in enumerate(items):
        if label_angles[i] is None:
            label_angles[i] = angle

    for (wedge, val, pct, w_angle, lbl), lbl_angle, rank in zip(items, label_angles, small_rank):
        if val == 0:
            continue
        # Implementation detail; see the repository documentation.
        # Implementation detail; see the repository documentation.
        # Implementation detail; see the repository documentation.
        override = label_overrides.get(lbl, {})
        lbl_angle = override.get("angle", lbl_angle)

        text  = f"{pct:.1f}%\n(n={val})"
        rad_w = math.radians(w_angle)
        rad_l = math.radians(lbl_angle)
        cos_w, sin_w = math.cos(rad_w), math.sin(rad_w)
        cos_l, sin_l = math.cos(rad_l), math.sin(rad_l)

        if pct >= SMALL_PCT:
            # Implementation detail; see the repository documentation.
            # Implementation detail; see the repository documentation.
            ax.text(R_INSIDE * cos_w, R_INSIDE * sin_w, text,
                    ha="center", va="center", fontsize=LABEL_FONTSIZE)
        else:
            r_base = R_LINE_END_LONG if pct < TINY_PCT else R_LINE_END_SHORT
            # Implementation detail; see the repository documentation.
            r_base += R_STEP * (rank % 2)
            r_base = override.get("r_end", r_base)
            x_text = r_base * cos_l
            y_text = r_base * sin_l
            if cos_l > 0.18:
                ha = "left"
                x_text += 0.025
            elif cos_l < -0.18:
                ha = "right"
                x_text -= 0.025
            else:
                ha = "center"

            # Implementation detail; see the repository documentation.
            # Implementation detail; see the repository documentation.
            # Implementation detail; see the repository documentation.
            if sin_l > 0.35:
                va = "bottom"
                y_text += 0.015
            elif sin_l < -0.35:
                va = "top"
                y_text -= 0.015
            else:
                va = "center"
            ax.annotate(
                text,
                xy=(1.11 * cos_w, 1.11 * sin_w),
                xytext=(x_text, y_text),
                ha=ha, va=va, fontsize=LABEL_FONTSIZE,
                arrowprops=dict(arrowstyle="-", color="black", lw=0.7,
                                shrinkA=3, shrinkB=0),
            )

    ax.set_axis_off()


def add_row_legend(legend_ax, colors_dict, ncol=3):
    legend_ax.axis('off')
    handles = [mpatches.Patch(facecolor=c, label=l) for l, c in colors_dict.items()]
    legend_ax.legend(
        handles=handles, loc='center', ncol=ncol,
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
# Implementation detail; see the repository documentation.
# PANEL_LABEL_OVERRIDES = {
#     "HLA II PtxS1": {
#         "IEDB II":   {"angle": 60},
#         "NetMHCII":  {"angle": 100},
#     },
# }
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
                  left=0.035, right=0.985, top=0.90, bottom=0.025,
                  hspace=0.38, wspace=0.18, figure=fig)

    axes_row_I   = [fig.add_subplot(gs[0, c]) for c in range(3)]
    legend_ax_I  = fig.add_subplot(gs[1, :])
    axes_row_II  = [fig.add_subplot(gs[2, c]) for c in range(3)]
    legend_ax_II = fig.add_subplot(gs[3, :])

    letters = ["A", "B", "C", "D", "E", "F"]

    for c, antigen in enumerate(ANTIGENS):
        segments_I, total_I, _, _ = data[antigen]
        title = f"HLA I {antigen}"
        draw_donut_ax(axes_row_I[c], segments_I, total_I, MHC_I_COLORS,
                      title, letters[c],
                      label_overrides=PANEL_LABEL_OVERRIDES.get(title))

    for c, antigen in enumerate(ANTIGENS):
        _, _, segments_II, total_II = data[antigen]
        title = f"HLA II {antigen}"
        draw_donut_ax(axes_row_II[c], segments_II, total_II, MHC_II_COLORS,
                      title, letters[3 + c],
                      label_overrides=PANEL_LABEL_OVERRIDES.get(title))

    # Implementation detail; see the repository documentation.
    # Output generation.
    # Implementation detail; see the repository documentation.
    fig.canvas.draw()
    for c, antigen in enumerate(ANTIGENS):
        pos = axes_row_I[c].get_position()
        fig.text((pos.x0 + pos.x1) / 2, 0.975, f"HLA I {antigen}",
                 fontsize=TITLE_FONTSIZE, fontweight=TITLE_FONTWEIGHT,
                 ha="center", va="top")
    for c, antigen in enumerate(ANTIGENS):
        pos = axes_row_II[c].get_position()
        fig.text((pos.x0 + pos.x1) / 2, pos.y1 + 0.012,
                 f"HLA II {antigen}",
                 fontsize=TITLE_FONTSIZE, fontweight=TITLE_FONTWEIGHT,
                 ha="center", va="bottom")

    add_row_legend(legend_ax_I, MHC_I_COLORS, ncol=len(MHC_I_COLORS))
    add_row_legend(legend_ax_II, MHC_II_COLORS, ncol=len(MHC_II_COLORS))

    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    fig.savefig(out_path, dpi=DPI, facecolor='white', edgecolor='none')
    plt.close(fig)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    data = {}
    for ag in ANTIGENS:
        print(f"\n{'='*50}")
        print(f"Processing details{ag}")
        data[ag] = get_segments_for_antigen(ag)

    out_path = os.path.join(OUT_DIR, "donut_coverage_all.png")
    render_combined_figure(data, out_path)
    print(f"Saved output{out_path}")
    print("Completed successfully")


if __name__ == "__main__":
    main()
