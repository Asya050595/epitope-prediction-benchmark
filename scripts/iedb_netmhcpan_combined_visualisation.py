#!/usr/bin/env python3
"""
Combined visualisation of recovered matches for IEDB and DTU predictors.

The script preserves the data-reading and plotting logic of:
  - iedb_visualisation.py
  - netmhcpan_visualisation.py

It reads the same four comparison workbooks and creates one figure containing
the eight former plots:

  A  IEDB                         B  NetMHC(pan)
     HLA I | HLA II                  HLA I | HLA II

  C  IEDB overlap                 D  NetMHCpan overlap
     HLA I | HLA II                  HLA I | HLA II

Output: IEDB_NetMHCpan_RM_combined.png (600 dpi).

The default project root can be overridden with RM_ROOT. The output directory
can be overridden with RM_OUT_DIR.
"""

from project_paths import DATA_ROOT

import os
import re
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from matplotlib.ticker import MaxNLocator


# -----------------------------------------------------------------------------
# Paths and global figure settings
# -----------------------------------------------------------------------------

DEFAULT_ROOT = f"{DATA_ROOT}"
ROOT = Path(os.environ.get("RM_ROOT", DEFAULT_ROOT))
OUT_DIR = Path(
    os.environ.get(
        "RM_OUT_DIR",
        ROOT / "p24" / "Visualization",
    )
)

ANTIGENS = ["p24", "pp65", "PtxS1"]

IEDB_SUMMARY_I = (
    ROOT / "p24" / "Matches MHC I comparison"
    / "IEDB_I_match_comparison_SUMMARY.xlsx"
)
IEDB_SUMMARY_II = (
    ROOT / "p24" / "Matches MHC II comparison"
    / "IEDB_II_match_comparison_SUMMARY.xlsx"
)
DTU_SUMMARY_I = (
    ROOT / "p24" / "Matches MHC I comparison"
    / "NetMHCpan_4.1_match_comparison_ALL.xlsx"
)
DTU_SUMMARY_II = (
    ROOT / "p24" / "Matches MHC II comparison"
    / "NetMHCIIpan_4.1_match_comparison_ALL.xlsx"
)

OUTPUT_NAME = "IEDB_NetMHCpan_RM_combined.png"
OUTPUT_DPI = 600
# Eight information-dense plots with genuine 10-pt text require a landscape
# canvas wider than a standard 180-mm two-column figure. Reducing this canvas
# to 180 mm would also reduce the effective type size.
FIG_WIDTH_MM = 360
FIG_HEIGHT_MM = 285
FONT_SIZE = 10
# Match the final pie-chart figure exactly.
TITLE_FONTSIZE = 15
LETTER_FONTSIZE = 15

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"],
    "font.size": FONT_SIZE,
    "axes.titlesize": FONT_SIZE,
    "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE,
    "ytick.labelsize": FONT_SIZE,
    "legend.fontsize": FONT_SIZE,
})


# -----------------------------------------------------------------------------
# IEDB: read the comparison summaries exactly as in iedb_visualisation.py
# -----------------------------------------------------------------------------

def _value(value, fallback=0):
    """Safely convert an Excel cell to the requested numeric type."""
    try:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return fallback
        return type(fallback)(value)
    except (TypeError, ValueError):
        return fallback


def read_iedb_summary(path: Path) -> dict:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if not path.is_file():
        raise FileNotFoundError(f"Required input or value was not found{path}")

    df = pd.read_excel(path, sheet_name="Summary", header=None)
    all_stats = {}

    for _, row in df.iterrows():
        antigen = str(row.iloc[0]).strip()
        if antigen not in ANTIGENS:
            continue

        values = row.values
        union = _value(values[7], 0)
        all_stats[antigen] = {
            "union": union,
            "single": {
                "EL": _value(values[1], 0),
                "BA": _value(values[3], 0),
                "Consensus": _value(values[5], 0),
            },
            "single_pct": {
                "EL": _value(values[2], 0.0),
                "BA": _value(values[4], 0.0),
                "Consensus": _value(values[6], 0.0),
            },
            "venn": {
                "Only EL": _value(values[8], 0),
                "Only BA": _value(values[10], 0),
                "Only Consensus": _value(values[12], 0),
                "EL + BA": _value(values[14], 0),
                "EL + Consensus": _value(values[16], 0),
                "BA + Consensus": _value(values[18], 0),
                "EL + BA + Consensus": _value(values[20], 0),
            },
            "venn_pct": {
                "Only EL": _value(values[9], 0.0),
                "Only BA": _value(values[11], 0.0),
                "Only Consensus": _value(values[13], 0.0),
                "EL + BA": _value(values[15], 0.0),
                "EL + Consensus": _value(values[17], 0.0),
                "BA + Consensus": _value(values[19], 0.0),
                "EL + BA + Consensus": _value(values[21], 0.0),
            },
        }
        print(
            f"  {antigen}: union={union}, "
            f"EL={all_stats[antigen]['single']['EL']}, "
            f"BA={all_stats[antigen]['single']['BA']}, "
            f"Consensus={all_stats[antigen]['single']['Consensus']}"
        )

    return all_stats


def iedb_rows(stats: dict) -> list:
    """Return non-empty antigen rows in the fixed publication order."""
    return [
        {"antigen": antigen, **stats[antigen]}
        for antigen in ANTIGENS
        if antigen in stats and stats[antigen]["union"] > 0
    ]


# -----------------------------------------------------------------------------
# DTU: read the comparison summaries exactly as in netmhcpan_visualisation.py
# -----------------------------------------------------------------------------

ANTIGEN_COLUMNS = {
    "p24": (1, 2),
    "pp65": (3, 4),
    "PtxS1": (5, 6),
}

ROWS_I = {
    "total": 4,
    "exact1": 6,
    "ge2": 8,
    "union": 10,
    "el": 12,
    "only_el": 13,
    "ba": 14,
    "only_ba": 15,
    "both": 16,
}

ROWS_II = {
    "total": 4,
    "exact1": 6,
    "ge2": 8,
    "union": 10,
    "el": 12,
    "ba": 13,
    "aff": 14,
    "only_el": 16,
    "only_ba": 17,
    "only_aff": 18,
    "el_ba": 19,
    "el_aff": 20,
    "ba_aff": 21,
    "all3": 22,
}


def _parse_percentage(value) -> float:
    """Convert '92.1%' or 0.921 to a unit fraction."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return 0.0
    if isinstance(value, str):
        match = re.search(r"[\d.]+", value.strip().replace(",", "."))
        return float(match.group()) / 100.0 if match else 0.0
    return float(value)


def _to_int(value) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def read_dtu_mhc_i(path: Path) -> dict:
    """Read NetMHCpan MHC-I EL/BA summary data."""
    if not path.is_file():
        raise FileNotFoundError(f"Required input or value was not found{path}")

    df = pd.read_excel(path, sheet_name=0, header=None)
    result = {}
    for antigen, (count_col, pct_col) in ANTIGEN_COLUMNS.items():
        def count(key):
            return _to_int(df.iloc[ROWS_I[key], count_col])

        def pct(key):
            return _parse_percentage(df.iloc[ROWS_I[key], pct_col])

        result[antigen] = {
            "total": count("total"),
            "exact1": count("exact1"),
            "exact1_pct": pct("exact1"),
            "ge2": count("ge2"),
            "ge2_pct": pct("ge2"),
            "union": count("union"),
            "union_pct": pct("union"),
            "el": count("el"),
            "el_pct": pct("el"),
            "only_el": count("only_el"),
            "only_el_pct": pct("only_el"),
            "ba": count("ba"),
            "ba_pct": pct("ba"),
            "only_ba": count("only_ba"),
            "only_ba_pct": pct("only_ba"),
            "both": count("both"),
            "both_pct": pct("both"),
        }
        print(
            f"  HLA I {antigen}: total={result[antigen]['total']}, "
            f"EL={result[antigen]['el']}, BA={result[antigen]['ba']}, "
            f"EL+BA={result[antigen]['both']}"
        )
    return result


def read_dtu_mhc_ii(path: Path) -> dict:
    """Read NetMHCIIpan MHC-II EL/BA/affinity summary data."""
    if not path.is_file():
        raise FileNotFoundError(f"Required input or value was not found{path}")

    df = pd.read_excel(path, sheet_name=0, header=None)
    result = {}
    for antigen, (count_col, pct_col) in ANTIGEN_COLUMNS.items():
        def count(key):
            return _to_int(df.iloc[ROWS_II[key], count_col])

        def pct(key):
            return _parse_percentage(df.iloc[ROWS_II[key], pct_col])

        result[antigen] = {
            "total": count("total"),
            "exact1": count("exact1"),
            "exact1_pct": pct("exact1"),
            "ge2": count("ge2"),
            "ge2_pct": pct("ge2"),
            "union": count("union"),
            "union_pct": pct("union"),
            "el": count("el"),
            "el_pct": pct("el"),
            "ba": count("ba"),
            "ba_pct": pct("ba"),
            "aff": count("aff"),
            "aff_pct": pct("aff"),
            "only_el": count("only_el"),
            "only_el_pct": pct("only_el"),
            "only_ba": count("only_ba"),
            "only_ba_pct": pct("only_ba"),
            "only_aff": count("only_aff"),
            "only_aff_pct": pct("only_aff"),
            "el_ba": count("el_ba"),
            "el_ba_pct": pct("el_ba"),
            "el_aff": count("el_aff"),
            "el_aff_pct": pct("el_aff"),
            "ba_aff": count("ba_aff"),
            "ba_aff_pct": pct("ba_aff"),
            "all3": count("all3"),
            "all3_pct": pct("all3"),
        }
        print(
            f"  HLA II {antigen}: total={result[antigen]['total']}, "
            f"EL={result[antigen]['el']}, BA={result[antigen]['ba']}, "
            f"Aff={result[antigen]['aff']}, all3={result[antigen]['all3']}"
        )
    return result


# -----------------------------------------------------------------------------
# Original colour schemes
# -----------------------------------------------------------------------------

ANTIGEN_COLORS = {
    "p24": "#E8926B",
    "pp65": "#8FBA8F",
    "PtxS1": "#F0C674",
}

IEDB_VENN_KEYS = [
    "Only EL",
    "Only BA",
    "Only Consensus",
    "EL + BA",
    "EL + Consensus",
    "BA + Consensus",
    "EL + BA + Consensus",
]
IEDB_VENN_COLORS = [
    "#E8A598",
    "#C4A8D4",
    "#9FC4D4",
    "#F2C9A0",
    "#A8D4C2",
    "#B8BDD8",
    "#7A6E94",
]

COLOR_ONLY_EL = "#E8A598"
COLOR_BOTH_I = "#A8D4C2"
COLOR_ONLY_BA = "#C4A8D4"

DTU_VENN_II_ORDER = [
    "only_el",
    "only_ba",
    "only_aff",
    "el_ba",
    "el_aff",
    "ba_aff",
    "all3",
]
DTU_VENN_II_COLORS = {
    "only_el": "#E8A598",
    "only_ba": "#C4A8D4",
    "only_aff": "#9FC4D4",
    "el_ba": "#F2C9A0",
    "el_aff": "#A8D4C2",
    "ba_aff": "#B8BDD8",
    "all3": "#7A6E94",
}
DTU_VENN_II_LABELS = {
    "only_el": "Only EL",
    "only_ba": "Only BA",
    "only_aff": "Only Aff",
    "el_ba": "EL + BA",
    "el_aff": "EL + Aff",
    "ba_aff": "BA + Aff",
    "all3": "EL + BA + Aff",
}


def _text_color(hex_color: str, threshold: float = 0.4) -> str:
    """Choose black or white text from the linearized sRGB luminance."""
    value = hex_color.lstrip("#")
    red, green, blue = (
        int(value[index:index + 2], 16) / 255.0 for index in (0, 2, 4)
    )

    def linear(channel):
        if channel <= 0.04045:
            return channel / 12.92
        return ((channel + 0.055) / 1.055) ** 2.4

    luminance = (
        0.2126 * linear(red)
        + 0.7152 * linear(green)
        + 0.0722 * linear(blue)
    )
    return "black" if luminance >= threshold else "white"


def normalize_percentages(values, total, decimals=1):
    """Round percentages so that the displayed values sum to exactly 100%."""
    if total <= 0 or not values:
        raise ValueError(
            "Processing details"
            "Processing details"
        )
    if any(value < 0 for value in values):
        raise ValueError("Processing details")

    scale = 10 ** decimals
    raw = [value / total * 100 * scale for value in values]
    floored = [math.floor(value) for value in raw]
    remainder = round(100 * scale) - sum(floored)
    order = sorted(
        range(len(values)),
        key=lambda index: raw[index] - floored[index],
        reverse=True,
    )

    for index in range(max(remainder, 0)):
        floored[order[index % len(order)]] += 1
    for index in range(max(-remainder, 0)):
        floored[order[-(index + 1)]] -= 1

    expected_sum = 100 * scale
    if sum(floored) != expected_sum:
        raise ArithmeticError(
            "Processing details"
        )
    return [value / scale for value in floored]


# -----------------------------------------------------------------------------
# Plot helpers: same data and chart types, adapted to supplied Axes
# -----------------------------------------------------------------------------

def _style_axes(ax, horizontal=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    ax.grid(False)
    ax.set_axisbelow(False)
    if horizontal:
        ax.tick_params(axis="y", length=3, left=True)
    else:
        ax.tick_params(axis="x", length=0)


def _antigen_handles(antigens):
    return [
        mpatches.Patch(
            facecolor=ANTIGEN_COLORS[antigen],
            edgecolor="black",
            linewidth=0.7,
            label=antigen,
        )
        for antigen in antigens
    ]


def draw_grouped_bars(ax, data: dict, keys: list, labels: list,
                      show_ylabel: bool = True):
    """Draw the original grouped count bars for three antigens."""
    antigens = [antigen for antigen in ANTIGENS if antigen in data]
    x_values = np.arange(len(keys))
    # Use most of the available category width so that the three antigen bars
    # remain visually substantial in the combined multi-panel figure.
    group_width = 0.90
    bar_width = group_width / max(len(antigens), 1)
    offsets = (
        np.arange(len(antigens)) - (len(antigens) - 1) / 2
    ) * bar_width
    max_count = max(
        (data[antigen].get(key, 0) for antigen in antigens for key in keys),
        default=1,
    )

    for index, antigen in enumerate(antigens):
        counts = [data[antigen].get(key, 0) for key in keys]
        percentages = [data[antigen].get(f"{key}_pct", 0.0) for key in keys]
        bars = ax.bar(
            x_values + offsets[index],
            counts,
            width=bar_width * 0.94,
            color=ANTIGEN_COLORS[antigen],
            edgecolor="black",
            linewidth=0.7,
        )
        for bar, count, percentage in zip(bars, counts, percentages):
            if count == 0:
                continue
            # Shift the pp65 and PtxS1 labels slightly to the right so their
            # rotated text does not visually run into the preceding bar.
            label_dx = {"p24": 0, "pp65": 3, "PtxS1": 5}.get(antigen, 0)
            ax.annotate(
                f"{count} ({percentage * 100:.1f}%)",
                xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                xytext=(label_dx, 4),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=FONT_SIZE,
                rotation=90,
            )

    ax.set_xticks(x_values)
    ax.set_xticklabels(labels)
    if show_ylabel:
        ax.set_ylabel("Number of matches (n)")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    ax.set_ylim(0, max(max_count * 1.60, 1))
    _style_axes(ax)
    return _antigen_handles(antigens)


def draw_iedb_bars(ax, rows: list, show_ylabel: bool = True):
    """Draw IEDB EL/BA/Consensus bars from the original single statistics."""
    data = {
        row["antigen"]: {
            "EL": row["single"].get("EL", 0),
            "EL_pct": row["single_pct"].get("EL", 0.0),
            "BA": row["single"].get("BA", 0),
            "BA_pct": row["single_pct"].get("BA", 0.0),
            "Consensus": row["single"].get("Consensus", 0),
            "Consensus_pct": row["single_pct"].get("Consensus", 0.0),
        }
        for row in rows
    }
    return draw_grouped_bars(
        ax,
        data,
        ["EL", "BA", "Consensus"],
        ["EL", "BA", "Consensus"],
        show_ylabel=show_ylabel,
    )


def _annotate_wide_segment(ax, x_center, y_center, width, count, color):
    """Place a two-line label inside a segment that is wide enough."""
    label = f"{width * 100:.1f}%\n(n={count})"
    ax.text(
        x_center,
        y_center,
        label,
        ha="center",
        va="center",
        color=_text_color(color),
        fontsize=FONT_SIZE,
    )


def _short_callout_layout(
    items,
    x_min_by_key=None,
    x_shift_by_key=None,
    default_x_min=0.09,
    min_text_gap=0.22,
):
    """Return minimally separated local x positions for callouts.

    Labels stay close to the segment they describe.  When neighbouring labels
    would touch, their x positions are spread only as far as required.  This
    retains the clean single annotation level used by the donut charts while
    avoiding the very long connectors produced by fixed global positions.
    """
    if not items:
        return []

    x_min_by_key = x_min_by_key or {}
    x_shift_by_key = x_shift_by_key or {}
    positions = []
    minimums = []
    for key, center, width, count in items:
        x_min = x_min_by_key.get(key, default_x_min)
        minimums.append(x_min)
        preferred_x = center + x_shift_by_key.get(key, 0.0)
        positions.append(float(np.clip(preferred_x, x_min, 0.975)))

    # Forward pass: introduce only the separation required by adjacent text.
    for index in range(1, len(positions)):
        positions[index] = max(
            positions[index], positions[index - 1] + min_text_gap
        )

    # If the group reached the right boundary, compact it back from the right
    # without changing label order or crossing leader lines.
    if positions[-1] > 0.975:
        positions[-1] = 0.975
        for index in range(len(positions) - 2, -1, -1):
            positions[index] = min(
                positions[index], positions[index + 1] - min_text_gap
            )
        for index in range(len(positions)):
            positions[index] = max(positions[index], minimums[index])
            if index:
                positions[index] = max(
                    positions[index], positions[index - 1] + min_text_gap
                )

    return [
        (label_x, key, center, width, count)
        for label_x, (key, center, width, count) in zip(positions, items)
    ]


def draw_stacked_overlap(ax, row_names, segment_rows, segment_definitions,
                         expected_totals, xlabel,
                         callout_directions=None,
                         callout_x_min=None,
                         callout_x_shift=None,
                         callout_y_offset=None,
                         callout_ha=None,
                         callout_row_overrides=None,
                         default_callout_x_min=0.09):
    """Draw the original mutually exclusive overlap categories."""
    if not (
        len(row_names) == len(segment_rows) == len(expected_totals)
    ):
        raise ValueError("Table status")

    # Validate that mutually exclusive categories form a complete partition,
    # then apply largest-remainder rounding. Consequently, both the stacked
    # bar widths and the one-decimal labels total exactly 100.0% in every row.
    segment_keys = [definition[0] for definition in segment_definitions]
    normalized_rows = []
    for row_name, segments, expected_total in zip(
        row_names, segment_rows, expected_totals
    ):
        counts = [int(segments[key][0]) for key in segment_keys]
        category_total = sum(counts)
        expected_total = int(expected_total)
        if category_total != expected_total:
            raise ValueError(
                f"Validation status{row_name}Processing details"
                f"Processing details{category_total}Validation status"
                f"{expected_total}."
            )

        displayed_percentages = normalize_percentages(
            counts, expected_total, decimals=1
        )
        if not math.isclose(
            sum(displayed_percentages), 100.0, abs_tol=1e-9
        ):
            raise ArithmeticError(
                f"Processing details{row_name}Processing details"
            )
        normalized_rows.append({
            key: (count, percentage / 100.0)
            for key, count, percentage in zip(
                segment_keys, counts, displayed_percentages
            )
        })
    segment_rows = normalized_rows

    # Dedicated vertical lanes above and below every bar keep callouts from
    # different antigens separate even when many narrow segments are present.
    y_values = np.arange(len(row_names), dtype=float) * 2.75
    lefts = np.zeros(len(row_names), dtype=float)
    used_handles = []
    narrow_labels = [[] for _ in row_names]

    for segment_index, (key, color, label) in enumerate(segment_definitions):
        widths = np.asarray(
            [segments[key][1] for segments in segment_rows], dtype=float
        )
        counts = np.asarray(
            [segments[key][0] for segments in segment_rows], dtype=int
        )
        ax.barh(
            y_values,
            widths,
            left=lefts,
            height=0.72,
            color=color,
            edgecolor="black",
            linewidth=0.7,
        )

        if np.any(counts > 0):
            used_handles.append(
                mpatches.Patch(
                    facecolor=color,
                    edgecolor="black",
                    linewidth=0.7,
                    label=label,
                )
            )

        for row_index, (width, count) in enumerate(zip(widths, counts)):
            if count == 0:
                continue
            center = lefts[row_index] + width / 2
            # At 10 pt, labels in adjacent segments of about 20% can touch.
            # Keep text inside only when the segment provides a comfortable
            # horizontal margin; narrower labels use collision-free callouts.
            if width >= 0.23:
                _annotate_wide_segment(
                    ax, center, y_values[row_index], width, count, color
                )
            else:
                narrow_labels[row_index].append((key, center, width, count))
        lefts += widths

    # Labels from narrow segments remain complete.  As in the donut charts,
    # each callout begins at the outer edge of the coloured shape and ends at
    # the nearest edge of its label.  Nearby labels use compact vertical lanes
    # instead of distant fixed x positions, so the connectors stay short.
    callout_directions = callout_directions or {}
    callout_x_min = callout_x_min or {}
    callout_x_shift = callout_x_shift or {}
    callout_y_offset = callout_y_offset or {}
    callout_ha = callout_ha or {}
    callout_row_overrides = callout_row_overrides or {}
    for row_index, items in enumerate(narrow_labels):
        row_name = row_names[row_index]
        items.sort(key=lambda item: item[1])
        upper = []
        lower = []
        for item_index, item in enumerate(items):
            key = item[0]
            override = callout_row_overrides.get((row_name, key), {})
            default_direction = 1 if item_index % 2 == 0 else -1
            direction = override.get(
                "direction",
                callout_directions.get(key, default_direction),
            )
            (upper if direction > 0 else lower).append(item)
        for direction, group in ((1, upper), (-1, lower)):
            if not group:
                continue
            row_x_min = dict(callout_x_min)
            row_x_shift = dict(callout_x_shift)
            for key, _, _, _ in group:
                override = callout_row_overrides.get((row_name, key), {})
                if "x_min" in override:
                    row_x_min[key] = override["x_min"]
                if "x_shift" in override:
                    row_x_shift[key] = override["x_shift"]
            for label_x, key, center, width, count in _short_callout_layout(
                group,
                x_min_by_key=row_x_min,
                x_shift_by_key=row_x_shift,
                default_x_min=default_callout_x_min,
            ):
                override = callout_row_overrides.get((row_name, key), {})
                bar_edge_y = y_values[row_index] + direction * 0.36
                y_offset = override.get(
                    "y_offset", callout_y_offset.get(key, 0.90)
                )
                label_y = (
                    y_values[row_index]
                    + direction * y_offset
                )
                label_text = f"{width * 100:.1f}%\n(n={count})"
                if override.get("separate_text_from_leader", False):
                    # Keep the already approved leader line fixed while
                    # allowing a small independent text correction.
                    ax.plot(
                        [center, label_x],
                        [bar_edge_y, label_y],
                        color="black",
                        lw=0.8,
                        clip_on=False,
                    )
                    text_y_offset = override.get("text_y_offset", y_offset)
                    text_y = (
                        y_values[row_index]
                        + direction * text_y_offset
                    )
                    ax.text(
                        label_x + override.get("text_dx", 0.0),
                        text_y,
                        label_text,
                        ha=override.get("text_ha", "center"),
                        va=override.get(
                            "text_va",
                            "bottom" if direction > 0 else "top",
                        ),
                        fontsize=FONT_SIZE,
                        color="black",
                        bbox=(
                            {
                                "facecolor": "white",
                                "edgecolor": "none",
                                "pad": 0.0,
                            }
                            if override.get("mask_leader_behind_text", False)
                            else None
                        ),
                        clip_on=False,
                    )
                else:
                    ax.annotate(
                        label_text,
                        xy=(center, bar_edge_y),
                        xytext=(label_x, label_y),
                        textcoords="data",
                        ha=override.get("ha", callout_ha.get(key, "center")),
                        va="bottom" if direction > 0 else "top",
                        fontsize=FONT_SIZE,
                        color="black",
                        arrowprops={
                            "arrowstyle": "-",
                            "color": "black",
                            "lw": 0.8,
                            "shrinkA": 3,
                            "shrinkB": 0,
                        },
                        clip_on=False,
                    )

    ax.set_yticks(y_values)
    ax.set_yticklabels(row_names)
    # The stacked bars begin exactly at x=0, so they remain visually connected
    # to the y-axis. A small margin is retained only at the right edge.
    ax.set_xlim(0.0, 1.055)
    # The shared x-axis label is drawn once beneath the HLA-I/HLA-II pair.
    # Repeating it on both compact axes competes with the shared legend.
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", "0.5", "1.0"])
    if len(y_values):
        ax.set_ylim(y_values[0] - 1.65, y_values[-1] + 1.65)
    _style_axes(ax, horizontal=True)
    return used_handles


def draw_iedb_overlap(ax, rows: list, mhc_class: str):
    segment_rows = [
        {
            key: (row["venn"].get(key, 0), row["venn_pct"].get(key, 0.0))
            for key in IEDB_VENN_KEYS
        }
        for row in rows
    ]
    definitions = list(zip(IEDB_VENN_KEYS, IEDB_VENN_COLORS, IEDB_VENN_KEYS))
    callout_directions = {}
    callout_x_min = {}
    default_callout_x_min = 0.11
    if mhc_class == "I":
        # Preserve the requested category semantics: Only EL is above a short
        # leader and centred safely clear of the y-axis; EL + BA is below
        # every HLA-I bar. BA + Consensus receives a small rightward shift in
        # all three antigen rows.
        callout_directions = {"Only EL": 1, "EL + BA": -1}
        callout_x_min = {"Only EL": 0.17}

    return draw_stacked_overlap(
        ax,
        [row["antigen"] for row in rows],
        segment_rows,
        definitions,
        [row["union"] for row in rows],
        "Proportion of union pairs",
        callout_directions=callout_directions,
        callout_x_min=callout_x_min,
        callout_x_shift=(
            {"BA + Consensus": 0.14}
            if mhc_class == "I" else {}
        ),
        callout_y_offset=(
            {"Only EL": 0.68} if mhc_class == "I" else {}
        ),
        callout_row_overrides=(
            {
                ("PtxS1", "EL + BA"): {"y_offset": 0.78},
                ("PtxS1", "Only Consensus"): {
                    "separate_text_from_leader": True,
                    "text_dx": 0.04,
                    "text_y_offset": 0.90,
                    "text_ha": "center",
                    "text_va": "bottom",
                },
                ("pp65", "EL + BA"): {"y_offset": 0.78},
                ("pp65", "Only EL"): {
                    "x_min": 0.22,
                    "y_offset": 0.54,
                    "separate_text_from_leader": True,
                    "text_dx": -0.04,
                    "text_y_offset": 0.44,
                    "text_ha": "center",
                    "text_va": "bottom",
                    "mask_leader_behind_text": True,
                },
                ("pp65", "BA + Consensus"): {
                    "x_shift": 0.28,
                    "y_offset": 0.64,
                },
                ("p24", "Only EL"): {
                    "x_min": 0.22,
                    "y_offset": 0.54,
                    "separate_text_from_leader": True,
                    "text_dx": -0.04,
                    "text_y_offset": 0.44,
                    "text_ha": "center",
                    "text_va": "bottom",
                    "mask_leader_behind_text": True,
                },
                ("p24", "BA + Consensus"): {
                    "x_shift": 0.28,
                    "y_offset": 0.64,
                },
            }
            if mhc_class == "I" else {}
        ),
        default_callout_x_min=default_callout_x_min,
    )


def draw_dtu_overlap_i(ax, stats: dict):
    antigens = [antigen for antigen in ANTIGENS if antigen in stats]
    segment_rows = []
    for antigen in antigens:
        values = stats[antigen]
        total = values["total"] if values["total"] > 0 else 1
        segment_rows.append({
            "only_el": (values["only_el"], values["only_el"] / total),
            "both": (values["both"], values["both"] / total),
            "only_ba": (values["only_ba"], values["only_ba"] / total),
        })
    definitions = [
        ("only_el", COLOR_ONLY_EL, "Only EL"),
        ("both", COLOR_BOTH_I, "EL + BA"),
        ("only_ba", COLOR_ONLY_BA, "Only BA"),
    ]
    return draw_stacked_overlap(
        ax,
        antigens,
        segment_rows,
        definitions,
        [stats[antigen]["total"] for antigen in antigens],
        "Proportion of total RM",
        callout_directions={"only_el": 1},
        callout_x_min={"only_el": 0.14},
        default_callout_x_min=0.11,
    )


def draw_dtu_overlap_ii(ax, stats: dict):
    antigens = [antigen for antigen in ANTIGENS if antigen in stats]
    segment_rows = []
    for antigen in antigens:
        values = stats[antigen]
        total = values["total"] if values["total"] > 0 else 1
        segment_rows.append({
            key: (values.get(key, 0), values.get(key, 0) / total)
            for key in DTU_VENN_II_ORDER
        })
    definitions = [
        (key, DTU_VENN_II_COLORS[key], DTU_VENN_II_LABELS[key])
        for key in DTU_VENN_II_ORDER
    ]
    return draw_stacked_overlap(
        ax,
        antigens,
        segment_rows,
        definitions,
        [stats[antigen]["total"] for antigen in antigens],
        "Proportion of total RM",
        # Retain the current HLA-II arrangement for the two requested
        # categories regardless of which other narrow categories are present.
        callout_directions={"el_aff": 1, "ba_aff": 1},
        callout_x_shift={"ba_aff": 0.10},
        default_callout_x_min=0.11,
    )


# -----------------------------------------------------------------------------
# Combined 2 × 2 panel composition
# -----------------------------------------------------------------------------

def _draw_panel_header(fig, spec, letter, title):
    """Draw a bold panel letter at upper left and a centred bold title."""
    header = fig.add_subplot(spec)
    header.axis("off")
    header.text(
        0.00,
        0.55,
        letter,
        ha="left",
        va="center",
        fontsize=LETTER_FONTSIZE,
        fontweight="bold",
    )
    header.text(
        0.54,
        0.55,
        title,
        ha="center",
        va="center",
        fontsize=TITLE_FONTSIZE,
        fontweight="bold",
    )



def _new_chart_axis(fig, spec, title):
    axis = fig.add_subplot(spec)
    axis.set_title(title, fontweight="bold", pad=7)
    return axis


def _new_legend_axis(fig, spec):
    axis = fig.add_subplot(spec)
    axis.axis("off")
    return axis


def _draw_legend(axis, handles, columns=1, anchor_y=1.0):
    if not handles:
        return
    axis.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, anchor_y),
        ncol=columns,
        frameon=False,
        handlelength=1.0,
        handletextpad=0.45,
        columnspacing=0.8,
        borderaxespad=0.0,
    )


def make_combined_figure(iedb_i, iedb_ii, dtu_i, dtu_ii, output_path: Path):
    """Create the single eight-chart figure requested for publication."""
    fig = plt.figure(
        figsize=(FIG_WIDTH_MM / 25.4, FIG_HEIGHT_MM / 25.4),
        facecolor="white",
    )
    # A direct 4-column grid gives every chart substantially more drawing area
    # than a 2x2 grid containing another nested 1x2 grid in each cell.
    grid = fig.add_gridspec(
        6,
        4,
        left=0.040,
        right=0.995,
        top=0.995,
        bottom=0.025,
        wspace=0.44,
        hspace=0.23,
        height_ratios=[0.12, 1.00, 0.15, 0.12, 1.35, 0.38],
    )

    rows_i = iedb_rows(iedb_i)
    rows_ii = iedb_rows(iedb_ii)

    # The overlap panels contain numerous external labels. Give each HLA-I /
    # HLA-II pair its own wider internal gutter so callouts cannot compete in
    # the centre, without altering the spacing of the upper bar charts.
    overlap_grid_c = grid[4, 0:2].subgridspec(1, 2, wspace=0.78)
    overlap_grid_d = grid[4, 2:4].subgridspec(1, 2, wspace=0.78)

    # Shared headers reproduce the requested A-D / title format.
    _draw_panel_header(fig, grid[0, 0:2], "A", "IEDB all matches")
    _draw_panel_header(
        fig, grid[0, 2:4], "B", "DTU NetMHC(II)pan all matches"
    )
    _draw_panel_header(fig, grid[3, 0:2], "C", "IEDB thresholds overlap")
    _draw_panel_header(
        fig, grid[3, 2:4], "D", "DTU NetMHC(II)pan thresholds overlap"
    )

    # A — IEDB totals by mode
    ax_i = _new_chart_axis(fig, grid[1, 0], "HLA I")
    ax_ii = _new_chart_axis(fig, grid[1, 1], "HLA II")
    legend = _new_legend_axis(fig, grid[2, 0:2])
    handles_i = draw_iedb_bars(ax_i, rows_i, show_ylabel=True)
    draw_iedb_bars(ax_ii, rows_ii, show_ylabel=False)
    _draw_legend(legend, handles_i, columns=3)

    # B — standalone NetMHC(pan) totals by mode
    ax_i = _new_chart_axis(fig, grid[1, 2], "HLA I")
    ax_ii = _new_chart_axis(fig, grid[1, 3], "HLA II")
    legend = _new_legend_axis(fig, grid[2, 2:4])
    handles_i = draw_grouped_bars(
        ax_i,
        dtu_i,
        ["el", "ba", "both"],
        ["EL", "BA", "EL + BA"],
        show_ylabel=True,
    )
    handles_ii = draw_grouped_bars(
        ax_ii,
        dtu_ii,
        ["el", "ba", "aff"],
        ["EL", "BA", "Affinity"],
        show_ylabel=False,
    )
    _draw_legend(legend, handles_i, columns=3)

    # C — IEDB mutually exclusive overlap categories
    ax_i = _new_chart_axis(fig, overlap_grid_c[0, 0], "HLA I")
    ax_ii = _new_chart_axis(fig, overlap_grid_c[0, 1], "HLA II")
    legend = _new_legend_axis(fig, grid[5, 0:2])
    handles_i = draw_iedb_overlap(ax_i, rows_i, "I")
    handles_ii = draw_iedb_overlap(ax_ii, rows_ii, "II")
    legend.text(
        0.5, 0.72, "Proportion of union pairs",
        ha="center", va="bottom", fontsize=FONT_SIZE,
    )
    _draw_legend(legend, handles_i, columns=3, anchor_y=0.40)

    # D — standalone NetMHC(pan) mutually exclusive overlap categories
    ax_i = _new_chart_axis(fig, overlap_grid_d[0, 0], "HLA I")
    ax_ii = _new_chart_axis(fig, overlap_grid_d[0, 1], "HLA II")
    legend = _new_legend_axis(fig, grid[5, 2:4])
    handles_i = draw_dtu_overlap_i(ax_i, dtu_i)
    handles_ii = draw_dtu_overlap_ii(ax_ii, dtu_ii)
    legend.text(
        0.5, 0.72, "Proportion of total RM",
        ha="center", va="bottom", fontsize=FONT_SIZE,
    )
    # The HLA-I and HLA-II overlap palettes are identical except that EL + BA
    # has a different colour in the two original plots. Retain both colours and
    # qualify those two legend entries explicitly.
    combined_handles = []
    seen = set()
    for handle in handles_i + handles_ii:
        label = handle.get_label()
        color = tuple(handle.get_facecolor())
        key = (label, color)
        if key in seen:
            continue
        seen.add(key)
        if label == "EL + BA":
            label = (
                "EL + BA (HLA I)"
                if color == tuple(to_rgba(COLOR_BOTH_I))
                else "EL + BA (HLA II)"
            )
        combined_handles.append(
            mpatches.Patch(
                facecolor=handle.get_facecolor(),
                edgecolor="black",
                linewidth=0.7,
                label=label,
            )
        )
    _draw_legend(legend, combined_handles, columns=3, anchor_y=0.40)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path,
        dpi=OUTPUT_DPI,
        facecolor="white",
        edgecolor="none",
    )
    plt.close(fig)
    print(f"Saved output{output_path}")


def main():
    print("=" * 72)
    print("Processing details")
    print("=" * 72)

    print(f"Processing details{IEDB_SUMMARY_I}")
    iedb_i = read_iedb_summary(IEDB_SUMMARY_I)

    print(f"Processing details{IEDB_SUMMARY_II}")
    iedb_ii = read_iedb_summary(IEDB_SUMMARY_II)

    print(f"Processing details{DTU_SUMMARY_I}")
    dtu_i = read_dtu_mhc_i(DTU_SUMMARY_I)

    print(f"Processing details{DTU_SUMMARY_II}")
    dtu_ii = read_dtu_mhc_ii(DTU_SUMMARY_II)

    make_combined_figure(
        iedb_i,
        iedb_ii,
        dtu_i,
        dtu_ii,
        OUT_DIR / OUTPUT_NAME,
    )
    print("Completed successfully")


if __name__ == "__main__":
    main()
