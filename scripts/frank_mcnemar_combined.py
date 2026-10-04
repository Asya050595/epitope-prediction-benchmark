#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Portable implementation of frank_mcnemar_combined. See the repository README and data/README.md for inputs, outputs, and execution instructions."""

from __future__ import annotations

from project_paths import DATA_ROOT

import itertools
import os
import re
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import openpyxl
import pandas as pd
from scipy.stats import binomtest


# ---------------------------------------------------------------------
# Configuration.
# ---------------------------------------------------------------------

DEFAULT_ROOT = Path(f"{DATA_ROOT}")
BASE_DIR = Path(os.environ.get("RM_ROOT", DEFAULT_ROOT))
OUTPUT_DIR = BASE_DIR / "p24" / "Visualization"
OUTPUT_FILE = OUTPUT_DIR / "FRANK_McNemar_HLA_I_II_combined.png"

ANTIGENS = ["p24", "pp65"]
MHC_CLASSES = ("I", "II")

FONT_FAMILY = "DejaVu Sans"
FONT_SIZE = 10
TITLE_SIZE = 11
PANEL_LETTER_SIZE = 11
DPI = 600

# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
MM_TO_INCH = 1.0 / 25.4
FIG_WIDTH_MM = 280
FIG_HEIGHT_MM = 225
FIGSIZE = (FIG_WIDTH_MM * MM_TO_INCH, FIG_HEIGHT_MM * MM_TO_INCH)

PANEL_TITLES = {
    "A": "HLA I",
    "B": "HLA II",
    "C": "HLA I",
    "D": "HLA II",
}

ROW_TITLES = {
    "mcnemar": "Pairwise McNemar Test of Conditional Sensitivity",
    "frank": "FRANK Score by Epitope Prediction Tool",
}


def mcnemar_input(antigen: str, mhc_class: str) -> Path:
    return (
        BASE_DIR
        / f"{antigen}"
        / f"McNemar test MHC {mhc_class}"
        / f"McNemar_{antigen}_MHC_{mhc_class}.xlsx"
    )


def frank_input(antigen: str, mhc_class: str) -> Path:
    return (
        BASE_DIR
        / f"{antigen}"
        / f"FRANK MHC {mhc_class}"
        / f"FRANK_{antigen}_MHC_{mhc_class}.xlsx"
    )


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [FONT_FAMILY, "Arial", "Liberation Sans"],
            "font.size": FONT_SIZE,
            "axes.titlesize": TITLE_SIZE,
            "axes.labelsize": FONT_SIZE,
            "xtick.labelsize": FONT_SIZE,
            "ytick.labelsize": FONT_SIZE,
            "legend.fontsize": FONT_SIZE,
            "axes.unicode_minus": False,
        }
    )


def add_panel_heading(ax, letter: str, title: str, y: float = 1.08) -> None:
    """Portable implementation of frank_mcnemar_combined. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    ax.text(
        -0.15,
        y,
        letter,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=PANEL_LETTER_SIZE,
        fontweight="bold",
        clip_on=False,
    )
    ax.set_title(title, fontsize=TITLE_SIZE, fontweight="bold", y=y, pad=0)


class PipelineError(Exception):
    """Portable implementation of frank_mcnemar_combined. See the repository README and data/README.md for inputs, outputs, and execution instructions."""


# ---------------------------------------------------------------------
# Input loading.
# ---------------------------------------------------------------------

# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
MC_GROUP_WIDTH = 3.60
MC_GROUP_SPACING = 5.20
MC_BAR_WIDTH_FRACTION = 0.80
MC_AXIS_BAR_GAP = 0.45

# Implementation detail; see the repository documentation.
COLOR_IEDB_CONSENSUS = "#E07A5F"
COLOR_IEDB_BA = "#81B29A"
COLOR_IEDB_EL = "#F2CC8F"
COLOR_NETMHCPAN = "#3D5A80"
COLOR_NETMHC = "#98C1D9"
COLOR_NETCTL = "#BC6C25"

MC_BRACKET_TOP_Y = 100.0
MC_BRACKET_STEP_Y = 7.0
MC_BRACKET_TEXT_GAP_Y = 0.01
MC_P_THRESHOLDS = [(0.001, "***"), (0.01, "**"), (0.05, "*")]


def _find_header_row(ws, header_value: str, col: int = 1, max_search: int = 15) -> int:
    for row in range(1, max_search + 1):
        if ws.cell(row=row, column=col).value == header_value:
            return row
    raise PipelineError(
        f"Не найдена строка заголовка '{header_value}' в листе "
        f"'{ws.title}' (проверены первые {max_search} строк)."
    )


def load_mcnemar_overview(path: Path) -> list[dict]:
    if not path.is_file():
        raise PipelineError(f"Не найден входной файл McNemar: {path}")

    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    if "Tool_overview" not in wb.sheetnames:
        raise PipelineError(f"В файле нет листа 'Tool_overview': {path}")
    ws = wb["Tool_overview"]
    header_row = _find_header_row(ws, "Tool")

    required = ("Tool", "Sensitivity (%)", "95% CI lower (%)", "95% CI upper (%)")
    columns = {
        ws.cell(row=header_row, column=c).value: c
        for c in range(1, ws.max_column + 1)
    }
    missing = [name for name in required if name not in columns]
    if missing:
        raise PipelineError(
            f"В листе 'Tool_overview' отсутствуют столбцы {missing}: {path}"
        )

    rows = []
    for row in range(header_row + 1, ws.max_row + 1):
        tool = ws.cell(row=row, column=columns["Tool"]).value
        sensitivity = ws.cell(row=row, column=columns["Sensitivity (%)"]).value
        if tool is None or sensitivity is None:
            continue
        sensitivity = float(sensitivity)
        lower = ws.cell(row=row, column=columns["95% CI lower (%)"]).value
        upper = ws.cell(row=row, column=columns["95% CI upper (%)"]).value
        rows.append(
            {
                "tool": str(tool).strip(),
                "sensitivity": sensitivity,
                "ci_lower": float(lower) if lower is not None else sensitivity,
                "ci_upper": float(upper) if upper is not None else sensitivity,
            }
        )

    rows.sort(key=lambda item: item["sensitivity"], reverse=True)
    return rows


def load_mcnemar_holm(path: Path) -> dict[frozenset, float]:
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    if "Pairwise_results" not in wb.sheetnames:
        raise PipelineError(f"В файле нет листа 'Pairwise_results': {path}")
    ws = wb["Pairwise_results"]
    header_row = _find_header_row(ws, "Comparison ID")
    columns = {
        ws.cell(row=header_row, column=c).value: c
        for c in range(1, ws.max_column + 1)
    }
    required = ("Tool 1", "Tool 2", "Holm-adjusted p")
    missing = [name for name in required if name not in columns]
    if missing:
        raise PipelineError(
            f"В листе 'Pairwise_results' отсутствуют столбцы {missing}: {path}"
        )

    values = {}
    for row in range(header_row + 1, ws.max_row + 1):
        tool_1 = ws.cell(row=row, column=columns["Tool 1"]).value
        tool_2 = ws.cell(row=row, column=columns["Tool 2"]).value
        p_value = ws.cell(row=row, column=columns["Holm-adjusted p"]).value
        if tool_1 is None or tool_2 is None or p_value is None:
            continue
        values[frozenset((str(tool_1).strip(), str(tool_2).strip()))] = float(p_value)
    return values


def mcnemar_tool_color(tool_name: str, context: str = "") -> str:
    is_iedb = re.search(r"\bIEDB\b", tool_name) is not None
    if is_iedb:
        if re.search(r"\bEL\b", tool_name):
            return COLOR_IEDB_EL
        if re.search(r"\bBA\b", tool_name):
            return COLOR_IEDB_BA
        if re.search(r"\bConsensus\b", tool_name, re.IGNORECASE):
            return COLOR_IEDB_CONSENSUS
    else:
        if re.search(r"NetMHC(?:II)?pan", tool_name):
            return COLOR_NETMHCPAN
        if re.search(r"\bNetMHC\b", tool_name):
            return COLOR_NETMHC
        if re.search(r"\bNetCTL\b", tool_name):
            return COLOR_NETCTL
    raise PipelineError(
        f"Для инструмента '{tool_name}' ({context}) не задано цветовое правило."
    )


def stars_for_mcnemar(p_value: float) -> str:
    for threshold, stars in MC_P_THRESHOLDS:
        if p_value < threshold:
            return stars
    return ""


def build_mcnemar_group(antigen: str, mhc_class: str) -> dict:
    path = mcnemar_input(antigen, mhc_class)
    overview = load_mcnemar_overview(path)
    tools = [row["tool"] for row in overview]
    holm = load_mcnemar_holm(path)

    significant = []
    pairs = list(itertools.combinations(range(len(tools)), 2))
    pairs.sort(key=lambda pair: (pair[1] - pair[0], pair[0]))
    for left, right in pairs:
        p_value = holm.get(frozenset((tools[left], tools[right])))
        if p_value is None:
            print(
                f"ВНИМАНИЕ: нет Holm p-value для '{tools[left]}' vs "
                f"'{tools[right]}' ({antigen}, HLA {mhc_class}); скобка пропущена.",
                file=sys.stderr,
            )
            continue
        stars = stars_for_mcnemar(p_value)
        # Implementation detail; see the repository documentation.
        # Implementation detail; see the repository documentation.
        if len(stars) >= 3:
            significant.append((left, right, stars))

    return {
        "antigen": antigen,
        "tools": tools,
        "sensitivities": [row["sensitivity"] for row in overview],
        "ci_lower": [row["ci_lower"] for row in overview],
        "ci_upper": [row["ci_upper"] for row in overview],
        "colors": [
            mcnemar_tool_color(tool, f"{antigen}, HLA {mhc_class}") for tool in tools
        ],
        "brackets": significant,
    }


def build_mcnemar_groups(mhc_class: str) -> list[dict]:
    groups = [build_mcnemar_group(antigen, mhc_class) for antigen in ANTIGENS]
    for group_index, group in enumerate(groups):
        center = group_index * MC_GROUP_SPACING
        n_tools = len(group["tools"])
        if n_tools <= 1:
            offsets = [0.0]
            bar_width = 0.22
        else:
            offsets = np.linspace(-MC_GROUP_WIDTH / 2, MC_GROUP_WIDTH / 2, n_tools)
            bar_width = float(
                (offsets[1] - offsets[0]) * MC_BAR_WIDTH_FRACTION
            )
        group["x_positions"] = [center + offset for offset in offsets]
        group["center"] = center
        group["bar_width"] = bar_width

    return groups


def mcnemar_heading_y(groups: list[dict]) -> float:
    """Portable implementation of frank_mcnemar_combined. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    return 1.08


def draw_mcnemar_panel(
    ax,
    groups: list[dict],
    letter: str,
    title: str,
    heading_y: Optional[float] = None,
) -> float:

    for group in groups:
        for index, (x, sensitivity, color) in enumerate(
            zip(group["x_positions"], group["sensitivities"], group["colors"])
        ):
            ax.bar(
                x,
                sensitivity,
                width=group["bar_width"],
                color=color,
                edgecolor="black",
                linewidth=0.8,
                zorder=3,
            )
            lower = group["ci_lower"][index]
            upper = group["ci_upper"][index]
            ax.errorbar(
                x,
                sensitivity,
                yerr=[
                    [max(0.0, sensitivity - lower)],
                    [max(0.0, upper - sensitivity)],
                ],
                fmt="none",
                ecolor="black",
                elinewidth=0.9,
                capsize=3,
                capthick=0.9,
                zorder=5,
            )
        bracket_count = len(group["brackets"])
        for level, (left, right, stars) in enumerate(group["brackets"]):
            # Implementation detail; see the repository documentation.
            # Implementation detail; see the repository documentation.
            # Implementation detail; see the repository documentation.
            # Output generation.
            y = MC_BRACKET_TOP_Y - (bracket_count - 1 - level) * MC_BRACKET_STEP_Y
            x_left = group["x_positions"][left]
            x_right = group["x_positions"][right]
            ax.plot(
                [x_left, x_right],
                [y, y],
                color="black",
                linewidth=0.9,
                clip_on=False,
                zorder=4,
            )
            ax.text(
                (x_left + x_right) / 2,
                y + MC_BRACKET_TEXT_GAP_Y,
                stars,
                ha="center",
                va="bottom",
                fontsize=FONT_SIZE,
                fontweight="bold",
                clip_on=False,
                zorder=4,
            )

    # Implementation detail; see the repository documentation.
    # Metric calculation.
    ax.set_ylim(0, 100)
    ax.set_yticks(range(0, 101, 20))
    ax.set_yticklabels([f"{value}%" for value in range(0, 101, 20)])
    ax.set_xticks([group["center"] for group in groups])
    ax.set_xticklabels([group["antigen"] for group in groups])
    ax.set_ylabel("Sensitivity")
    widest_bar = max(group["bar_width"] for group in groups)
    x_margin = MC_GROUP_WIDTH / 2 + widest_bar / 2 + MC_AXIS_BAR_GAP
    ax.set_xlim(
        groups[0]["center"] - x_margin,
        groups[-1]["center"] + x_margin,
    )
    ax.tick_params(axis="x", length=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if heading_y is None:
        heading_y = mcnemar_heading_y(groups)
    add_panel_heading(ax, letter, title, heading_y)
    return heading_y


# ---------------------------------------------------------------------
# Input loading.
# ---------------------------------------------------------------------

METHOD_SELECTION = {
    "I": [
        ("IEDB_I_Consensus", "Consensus_percentile"),
        ("IEDB_I_NetMHCpan_4.1_BA", "Median_percentile"),
        ("IEDB_I_NetMHCpan_4.1_EL", "Median_percentile"),
        ("NetMHCpan_4.1", "Rank_EL"),
        ("NetMHC_4.0", "Rank_pct"),
        ("NetCTL_1.2", "COMB"),
    ],
    "II": [
        ("IEDB_II_Consensus", "Consensus_percentile"),
        ("IEDB_II_NetMHCIIpan_4.1_BA", "Median_percentile"),
        ("IEDB_II_NetMHCIIpan_4.1_EL", "Median_percentile"),
        ("NetMHCIIpan_4.1", "Rank_EL"),
        ("NetMHC_II_2.3", "Rank_pct"),
    ],
}

TOOL_DISPLAY = {
    "I": {
        "IEDB_I_Consensus": "IEDB I Consensus",
        "IEDB_I_NetMHCpan_4.1_BA": "IEDB I NetMHCpan 4.1 BA",
        "IEDB_I_NetMHCpan_4.1_EL": "IEDB I NetMHCpan 4.1 EL",
        "NetMHCpan_4.1": "NetMHCpan 4.1",
        "NetMHC_4.0": "NetMHC 4.0",
        "NetCTL_1.2": "NetCTL 1.2",
    },
    "II": {
        "IEDB_II_Consensus": "IEDB II Consensus",
        "IEDB_II_NetMHCIIpan_4.1_BA": "IEDB II NetMHCIIpan 4.1 BA",
        "IEDB_II_NetMHCIIpan_4.1_EL": "IEDB II NetMHCIIpan 4.1 EL",
        "NetMHCIIpan_4.1": "NetMHCIIpan 4.1",
        "NetMHC_II_2.3": "NetMHC II 2.3",
    },
}

FRANK_PALETTE = {
    "I": {
        "IEDB_I_Consensus": "#E07A5F",
        "IEDB_I_NetMHCpan_4.1_BA": "#81B29A",
        "IEDB_I_NetMHCpan_4.1_EL": "#F2CC8F",
        "NetMHCpan_4.1": "#3D5A80",
        "NetMHC_4.0": "#98C1D9",
        "NetCTL_1.2": "#BC6C25",
    },
    "II": {
        "IEDB_II_Consensus": "#E07A5F",
        "IEDB_II_NetMHCIIpan_4.1_BA": "#81B29A",
        "IEDB_II_NetMHCIIpan_4.1_EL": "#F2CC8F",
        "NetMHCIIpan_4.1": "#3D5A80",
        "NetMHC_II_2.3": "#98C1D9",
    },
}

REQUIRED_SUMMARY_COLUMNS = ["Tool", "Score", "Sheet", "Status_flag"]
REQUIRED_DETAIL_COLUMNS = [
    "Allele",
    "Peptide",
    "Length",
    "Epitope_score",
    "N_better",
    "N_total_peptides",
    "N_missing_windows",
    "FRANK_%",
    "Status",
]

FRANK_ALPHA = 0.05
FRANK_STAR_THRESHOLDS = [
    (0.0001, "****"),
    (0.001, "***"),
    (0.01, "**"),
    (0.05, "*"),
]
# Implementation detail; see the repository documentation.
FRANK_BOX_SPAN = 3.15
FRANK_GROUP_SPACING = 4.75
FRANK_BOX_WIDTH_FRACTION = 0.79
FRANK_AXIS_BOX_GAP = 0.45
FRANK_POINT_JITTER = 0.55
FRANK_POINT_ALPHA = 0.35
FRANK_POINT_SIZE = 5
FRANK_Y_MIN = 0
FRANK_Y_MAX = 100
# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
FRANK_BRACKET_GAP = 3.0
FRANK_BRACKET_STEP = 6.0
FRANK_BRACKET_TICK = 1.2
FRANK_BRACKET_TEXT_GAP = 1.8
FRANK_BRACKET_MAX_Y = 90.0


@dataclass
class ToolData:
    tool: str
    score: str
    sheet: str
    antigen: str
    mhc_class: str
    df_full: pd.DataFrame
    df_ok: pd.DataFrame


@dataclass
class PairwiseResult:
    mhc_class: str
    antigen: str
    tool_a: str
    tool_b: str
    n_common_ok_pairs: int
    n_denominator_mismatch: int
    n_valid_pairs: int
    n_wins_a: int
    n_wins_b: int
    n_ties: int
    n_nontie: int
    p_raw: Optional[float]
    skip_reason: Optional[str]
    p_holm: Optional[float] = None
    stars: str = ""
    significant: bool = False


def load_frank_tool_data(
    antigen: str,
    mhc_class: str,
    path: Path,
    method_selection: list[tuple[str, str]],
) -> dict[str, ToolData]:
    if not path.is_file():
        raise PipelineError(f"Не найден входной файл FRANK: {path}")

    try:
        excel = pd.ExcelFile(path, engine="openpyxl")
    except Exception as exc:
        raise PipelineError(f"Не удалось открыть файл {path}: {exc}") from exc
    if "Summary" not in excel.sheet_names:
        raise PipelineError(f"В файле отсутствует лист 'Summary': {path}")

    summary = excel.parse("Summary")
    missing = [name for name in REQUIRED_SUMMARY_COLUMNS if name not in summary.columns]
    if missing:
        raise PipelineError(f"В листе 'Summary' отсутствуют столбцы {missing}: {path}")

    result = {}
    for tool, score in method_selection:
        matches = summary[(summary["Tool"] == tool) & (summary["Score"] == score)]
        if len(matches) != 1:
            available = summary.loc[summary["Tool"] == tool, "Score"].tolist()
            raise PipelineError(
                f"В {path} комбинация Tool='{tool}', Score='{score}' найдена "
                f"{len(matches)} раз. Доступные Score: {available}"
            )
        sheet = matches.iloc[0]["Sheet"]
        if not isinstance(sheet, str) or sheet not in excel.sheet_names:
            raise PipelineError(
                f"Для Tool='{tool}', Score='{score}' указан отсутствующий лист "
                f"'{sheet}' в файле {path}."
            )

        frame = excel.parse(sheet)
        missing = [name for name in REQUIRED_DETAIL_COLUMNS if name not in frame.columns]
        if missing:
            raise PipelineError(
                f"В листе '{sheet}' отсутствуют столбцы {missing}: {path}"
            )
        duplicates = frame.duplicated(
            subset=["Allele", "Peptide", "Length"], keep=False
        )
        if duplicates.any():
            raise PipelineError(
                f"В листе '{sheet}' найдены дубликаты Allele+Peptide+Length "
                f"({int(duplicates.sum())} строк): {path}"
            )
        ok = frame[frame["Status"] == "ok"].copy()
        if ok.empty:
            raise PipelineError(
                f"В листе '{sheet}' нет строк со Status == 'ok': {path}"
            )
        result[tool] = ToolData(
            tool=tool,
            score=score,
            sheet=sheet,
            antigen=antigen,
            mhc_class=mhc_class,
            df_full=frame,
            df_ok=ok,
        )
    return result


def pairwise_frank(
    data_a: pd.DataFrame,
    data_b: pd.DataFrame,
    mhc_class: str,
    antigen: str,
    tool_a: str,
    tool_b: str,
) -> PairwiseResult:
    key = ["Allele", "Peptide", "Length"]
    a = data_a[key + ["N_better", "N_total_peptides"]].copy()
    b = data_b[key + ["N_better", "N_total_peptides"]].copy()
    merged = a.merge(b, on=key, how="inner", suffixes=("_A", "_B"))
    n_common = len(merged)

    if n_common == 0:
        return PairwiseResult(
            mhc_class,
            antigen,
            tool_a,
            tool_b,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            None,
            "нет общих валидных пар",
        )

    mismatch = merged["N_total_peptides_A"] != merged["N_total_peptides_B"]
    valid = merged.loc[~mismatch]
    wins_a = int((valid["N_better_A"] < valid["N_better_B"]).sum())
    wins_b = int((valid["N_better_A"] > valid["N_better_B"]).sum())
    ties = int((valid["N_better_A"] == valid["N_better_B"]).sum())
    n_nontie = wins_a + wins_b

    if n_nontie == 0:
        return PairwiseResult(
            mhc_class,
            antigen,
            tool_a,
            tool_b,
            n_common,
            int(mismatch.sum()),
            len(valid),
            wins_a,
            wins_b,
            ties,
            0,
            None,
            "нет пар без ничьей",
        )

    p_value = binomtest(wins_a, n_nontie, p=0.5, alternative="two-sided").pvalue
    return PairwiseResult(
        mhc_class,
        antigen,
        tool_a,
        tool_b,
        n_common,
        int(mismatch.sum()),
        len(valid),
        wins_a,
        wins_b,
        ties,
        n_nontie,
        float(p_value),
        None,
    )


def holm_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values, kind="mergesort")
    adjusted = np.empty(len(p_values), dtype=float)
    running_max = 0.0
    for rank, index in enumerate(order):
        value = min(1.0, max(running_max, (len(p_values) - rank) * p_values[index]))
        adjusted[index] = value
        running_max = value
    return adjusted


def frank_stars(p_value: Optional[float]) -> str:
    if p_value is None or np.isnan(p_value):
        return ""
    for threshold, stars in FRANK_STAR_THRESHOLDS:
        if p_value < threshold:
            return stars
    return ""


def build_frank_context(mhc_class: str) -> dict:
    selection = METHOD_SELECTION[mhc_class]
    all_data = {}
    for antigen in ANTIGENS:
        path = frank_input(antigen, mhc_class)
        all_data[antigen] = load_frank_tool_data(
            antigen, mhc_class, path, selection
        )

    tools = [tool for tool, _ in selection]
    results = []
    for antigen in ANTIGENS:
        for tool_a, tool_b in itertools.combinations(tools, 2):
            results.append(
                pairwise_frank(
                    all_data[antigen][tool_a].df_ok,
                    all_data[antigen][tool_b].df_ok,
                    mhc_class,
                    antigen,
                    tool_a,
                    tool_b,
                )
            )

    valid_indices = [
        index for index, result in enumerate(results) if result.p_raw is not None
    ]
    if valid_indices:
        adjusted = holm_adjust(
            np.array([results[index].p_raw for index in valid_indices], dtype=float)
        )
        for index, p_holm in zip(valid_indices, adjusted):
            results[index].p_holm = float(p_holm)
            results[index].stars = frank_stars(float(p_holm))
            results[index].significant = bool(p_holm < FRANK_ALPHA)

    return {
        "tools": tools,
        "tool_display": TOOL_DISPLAY[mhc_class],
        "palette": FRANK_PALETTE[mhc_class],
        "all_data": all_data,
        "results": results,
    }


def _frank_offsets(n_tools: int) -> np.ndarray:
    if n_tools == 1:
        return np.array([0.0])
    return np.linspace(-FRANK_BOX_SPAN / 2, FRANK_BOX_SPAN / 2, n_tools)


def _upper_boxplot_whisker(values: np.ndarray) -> float:
    """Portable implementation of frank_mcnemar_combined. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    clean = np.asarray(values, dtype=float)
    clean = clean[np.isfinite(clean)]
    if clean.size == 0:
        return 0.0
    q1, q3 = np.percentile(clean, [25, 75])
    upper_fence = q3 + 1.5 * (q3 - q1)
    inside = clean[clean <= upper_fence]
    return float(np.max(inside)) if inside.size else float(q3)


def draw_frank_panel(
    ax,
    context: dict,
    letter: str,
    title: str,
) -> None:
    tools = context["tools"]
    all_data = context["all_data"]
    results = context["results"]
    palette = context["palette"]
    tool_display = context["tool_display"]

    offsets = _frank_offsets(len(tools))
    box_width = (
        (offsets[1] - offsets[0]) * FRANK_BOX_WIDTH_FRACTION
        if len(tools) > 1
        else 0.5
    )
    centers = np.arange(len(ANTIGENS)) * FRANK_GROUP_SPACING
    positions = {}
    whisker_tops = {}

    for antigen_index, antigen in enumerate(ANTIGENS):
        for tool_index, tool in enumerate(tools):
            x = centers[antigen_index] + offsets[tool_index]
            positions[(antigen, tool)] = x
            values = all_data[antigen][tool].df_ok["FRANK_%"].dropna().values
            if len(values) == 0:
                continue
            whisker_tops[(antigen, tool)] = _upper_boxplot_whisker(values)
            color = palette[tool]
            ax.boxplot(
                [values],
                positions=[x],
                widths=box_width,
                patch_artist=True,
                showfliers=False,
                zorder=3,
                boxprops={"facecolor": color, "edgecolor": "black", "linewidth": 0.8},
                medianprops={"color": "black", "linewidth": 1.1},
                whiskerprops={"color": "black", "linewidth": 0.8},
                capprops={"color": "black", "linewidth": 0.8},
            )
            # Implementation detail; see the repository documentation.
            seed = sum(ord(char) for char in f"{antigen}|{tool}")
            rng = np.random.default_rng(seed)
            jitter = rng.uniform(
                -box_width * FRANK_POINT_JITTER / 2,
                box_width * FRANK_POINT_JITTER / 2,
                size=len(values),
            )
            ax.scatter(
                x + jitter,
                values,
                s=FRANK_POINT_SIZE,
                color=color,
                edgecolor="black",
                linewidth=0.15,
                alpha=FRANK_POINT_ALPHA,
                zorder=2,
            )

    for antigen in ANTIGENS:
        significant = [
            result
            for result in results
            if result.antigen == antigen and result.significant
        ]
        significant.sort(
            key=lambda result: abs(
                positions[(antigen, result.tool_a)]
                - positions[(antigen, result.tool_b)]
            )
        )
        if not significant:
            continue

        # Implementation detail; see the repository documentation.
        # Output generation.
        # Implementation detail; see the repository documentation.
        bracket_levels = []
        for result in significant:
            compared_top = max(
                whisker_tops.get((antigen, result.tool_a), 0.0),
                whisker_tops.get((antigen, result.tool_b), 0.0),
            )
            y = compared_top + FRANK_BRACKET_GAP
            if bracket_levels:
                y = max(y, bracket_levels[-1] + FRANK_BRACKET_STEP)
            bracket_levels.append(y)

        overflow = max(0.0, bracket_levels[-1] - FRANK_BRACKET_MAX_Y)
        bracket_levels = [y - overflow for y in bracket_levels]

        for result, y in zip(significant, bracket_levels):
            x_left = positions[(antigen, result.tool_a)]
            x_right = positions[(antigen, result.tool_b)]
            if x_left > x_right:
                x_left, x_right = x_right, x_left

            ax.plot(
                [x_left, x_left, x_right, x_right],
                [
                    y,
                    y + FRANK_BRACKET_TICK,
                    y + FRANK_BRACKET_TICK,
                    y,
                ],
                color="black",
                linewidth=0.8,
                zorder=4,
                clip_on=True,
            )
            ax.text(
                (x_left + x_right) / 2,
                y + FRANK_BRACKET_TICK + FRANK_BRACKET_TEXT_GAP,
                result.stars,
                ha="center",
                va="bottom",
                fontsize=FONT_SIZE,
                fontweight="bold",
                clip_on=True,
                zorder=5,
                bbox={
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.85,
                    "pad": 0.10,
                },
            )

    ax.set_xticks(centers)
    ax.set_xticklabels(ANTIGENS)
    ax.set_ylabel("FRANK, %")
    ax.set_ylim(FRANK_Y_MIN, FRANK_Y_MAX)
    ax.set_yticks(range(0, 101, 20))
    x_margin = FRANK_BOX_SPAN / 2 + box_width / 2 + FRANK_AXIS_BOX_GAP
    ax.set_xlim(centers[0] - x_margin, centers[-1] + x_margin)
    ax.axhline(0, color="grey", linewidth=0.6, zorder=0)
    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    for boundary in (centers[:-1] + centers[1:]) / 2:
        ax.axvline(boundary, color="#D9D9D9", linewidth=0.7, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    add_panel_heading(ax, letter, title, 1.055)


# ---------------------------------------------------------------------
# Implementation detail; see the repository documentation.
# ---------------------------------------------------------------------

def add_shared_legend(fig: plt.Figure) -> None:
    """Portable implementation of frank_mcnemar_combined. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    entries = [
        ("IEDB Consensus", COLOR_IEDB_CONSENSUS),
        ("IEDB NetMHCpan / NetMHCIIpan BA", COLOR_IEDB_BA),
        ("IEDB NetMHCpan / NetMHCIIpan EL", COLOR_IEDB_EL),
        ("DTU NetMHCpan / NetMHCIIpan", COLOR_NETMHCPAN),
        ("NetMHC / NetMHCII", COLOR_NETMHC),
        ("NetCTL", COLOR_NETCTL),
    ]
    handles = [
        Patch(facecolor=color, edgecolor="black", linewidth=0.8, label=label)
        for label, color in entries
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.025),
        ncol=3,
        frameon=False,
        fontsize=FONT_SIZE,
        handlelength=1.2,
        columnspacing=1.4,
        labelspacing=0.45,
    )


def add_row_title(
    fig: plt.Figure,
    axes_row,
    panel_heading_ys,
    title: str,
) -> None:
    """Portable implementation of frank_mcnemar_combined. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    boxes = [ax.get_position() for ax in axes_row]
    x_center = (min(box.x0 for box in boxes) + max(box.x1 for box in boxes)) / 2
    heading_tops = [
        box.y0 + heading_y * box.height
        for box, heading_y in zip(boxes, panel_heading_ys)
    ]
    fig.text(
        x_center,
        max(heading_tops) + 0.035,
        title,
        ha="center",
        va="bottom",
        fontsize=TITLE_SIZE,
        fontweight="bold",
    )


def build_final_figure() -> plt.Figure:
    configure_style()
    print("Загрузка данных McNemar...")
    mcnemar_i = build_mcnemar_groups("I")
    mcnemar_ii = build_mcnemar_groups("II")

    print("Загрузка данных FRANK и расчёт попарных тестов + Holm...")
    frank_i = build_frank_context("I")
    frank_ii = build_frank_context("II")

    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE)
    shared_mcnemar_heading_y = max(
        mcnemar_heading_y(mcnemar_i),
        mcnemar_heading_y(mcnemar_ii),
    )
    mcnemar_heading_i = draw_mcnemar_panel(
        axes[0, 0],
        mcnemar_i,
        "A",
        PANEL_TITLES["A"],
        shared_mcnemar_heading_y,
    )
    mcnemar_heading_ii = draw_mcnemar_panel(
        axes[0, 1],
        mcnemar_ii,
        "B",
        PANEL_TITLES["B"],
        shared_mcnemar_heading_y,
    )
    draw_frank_panel(axes[1, 0], frank_i, "C", PANEL_TITLES["C"])
    draw_frank_panel(axes[1, 1], frank_ii, "D", PANEL_TITLES["D"])
    add_shared_legend(fig)

    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    fig.subplots_adjust(
        left=0.085,
        right=0.985,
        top=0.85,
        bottom=0.13,
        wspace=0.18,
        hspace=0.72,
    )

    add_row_title(
        fig,
        axes[0, :],
        [mcnemar_heading_i, mcnemar_heading_ii],
        ROW_TITLES["mcnemar"],
    )
    add_row_title(
        fig,
        axes[1, :],
        [1.055, 1.055],
        ROW_TITLES["frank"],
    )
    return fig


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figure = build_final_figure()
    figure.savefig(
        OUTPUT_FILE,
        dpi=DPI,
        facecolor="white",
    )
    plt.close(figure)
    print(f"Сохранено: {OUTPUT_FILE}")


if __name__ == "__main__":
    try:
        main()
    except PipelineError as exc:
        print(f"\n[ОШИБКА] {exc}", file=sys.stderr)
        sys.exit(1)
    except Exception:
        print("\n[НЕОЖИДАННАЯ ОШИБКА]", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)
