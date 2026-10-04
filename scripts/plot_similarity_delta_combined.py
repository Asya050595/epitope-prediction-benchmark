#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combined visualization of changes in pairwise similarity between epitope
prediction tools after adding unrecovered predictions (UP):

    Delta J = J(RM + UP) - J(RM)
    Delta O = O(RM + UP) - O(RM)

The script creates one publication-ready 2 x 2 figure:
    A - HLA I Change in Jaccard Similarity
    B - HLA II Change in Jaccard Similarity
    C - HLA I Change in Overlap Coefficient
    D - HLA II Change in Overlap Coefficient

Each panel contains separate heatmaps for p24 and pp65. Only tools listed in
LABELS_I or LABELS_II are retained. The output width is 300 mm and the image
is saved at 600 dpi using a 17 pt sans-serif font.

Expected input files, relative to the Control directory:
    Runs for p24/Jaccard Index MHC I/Jaccardp24_MHC_I.xlsx
    Runs for pp65/Jaccard Index MHC I/Jaccardpp65_MHC_I.xlsx
    ...and the corresponding HLA II and Overlap Coefficient files.

In the real project the Russian directory names are used; see input_path().
The root directory can be overridden with the RM_ROOT environment variable.
"""

from __future__ import annotations

from project_paths import DATA_ROOT

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import openpyxl
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch, Rectangle


# -----------------------------------------------------------------------------
# Paths and output settings
# -----------------------------------------------------------------------------

DEFAULT_ROOT = f"{DATA_ROOT}"
BASE_DIR = Path(os.environ.get("RM_ROOT", DEFAULT_ROOT))
OUTPUT_DIR = BASE_DIR / "p24" / "Visualization"
OUTPUT_FILE = "Similarity_delta_HLA_I_II_combined.png"

ANTIGENS = ("p24", "pp65")

FIGURE_WIDTH_MM = 300.0
# This height retains the previously enlarged square heatmaps; the extra width
# is used to separate the left colorbars from the long HLA II tool labels.
FIGURE_HEIGHT_MM = 405.0
MM_PER_INCH = 25.4
DPI = 600
FONT_SIZE = 17


# -----------------------------------------------------------------------------
# Tool labels. These dictionaries are also explicit allowlists: any tool not
# listed here is omitted from both the calculation and the visualization.
# -----------------------------------------------------------------------------

LABELS_I = {
    "IEDB_I_Consensus": "IEDB-C",
    "IEDB_I_NetMHCpan_4.1_BA": "IEDB-BA",
    "IEDB_I_NetMHCpan_4.1_EL": "IEDB-EL",
    "NetCTL_1.2": "NetCTL",
    "NetMHC_4.0": "NetMHC",
    "NetMHCpan_4.1": "NetMHCpan",
}

LABELS_II = {
    "IEDB_II_Consensus": "IEDB-C",
    "IEDB_II_NetMHCIIpan_4.1_BA": "IEDB-BA",
    "IEDB_II_NetMHCIIpan_4.1_EL": "IEDB-EL",
    "NetMHC_II_2.3": "NetMHCII",
    "NetMHCIIpan_4.1": "NetMHCIIpan",
}

LABELS_BY_CLASS = {"I": LABELS_I, "II": LABELS_II}


METRICS = {
    "jaccard": {
        "folder": "Jaccard Index MHC {mhc_class}",
        "filename": "Jaccard{antigen}_MHC_{mhc_class}.xlsx",
        "sheet_rm": "Jaccard_RM",
        "sheet_all": "Jaccard_RM_plus_UP",
    },
    "overlap": {
        "folder": "Overlap Coefficient MHC {mhc_class}",
        "filename": "Overlap{antigen}_MHC_{mhc_class}.xlsx",
        "sheet_rm": "Overlap_RM",
        "sheet_all": "Overlap_RM_plus_UP",
    },
}

PANEL_CONFIGS = (
    ("A", "I", "jaccard", "HLA I ΔJ"),
    ("B", "II", "jaccard", "HLA II ΔJ"),
    ("C", "I", "overlap", "HLA I ΔO"),
    ("D", "II", "overlap", "HLA II ΔO"),
)


# -----------------------------------------------------------------------------
# Pastel diverging palettes. The first color represents a negative change, the
# last a positive change, and zero is pure white.
# -----------------------------------------------------------------------------

PALETTES = {
    "pastel_blue_peach": (
        "#6F9FC1", "#BED7E6", "#FFFFFF", "#F1C5B4", "#D98E72"
    ),
    "pastel_teal_lavender": (
        "#6FAEA6", "#BEDCD7", "#FFFFFF", "#DDD1EA", "#A88CC7"
    ),
    "pastel_sage_rose": (
        "#7FA58B", "#C7D8CC", "#FFFFFF", "#E8C7CD", "#C98F99"
    ),
    "pastel_blue_gold": (
        "#7EA6C2", "#C8DCE8", "#FFFFFF", "#F2DEB3", "#D8AD62"
    ),
}

# Recommended default: restrained, readable, and less pink than the other
# options. Change only this key to recolor the entire figure.
ACTIVE_PALETTE = "pastel_blue_peach"

MISSING_COLOR = "#D0D3D6"
NORM_GAMMA = 0.45


def setup_style() -> None:
    """Apply one consistent publication style to every text element."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": FONT_SIZE,
            "axes.titlesize": FONT_SIZE,
            "axes.labelsize": FONT_SIZE,
            "xtick.labelsize": FONT_SIZE,
            "ytick.labelsize": FONT_SIZE,
            "legend.fontsize": FONT_SIZE,
            "figure.titlesize": FONT_SIZE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def make_colormap(name: str) -> LinearSegmentedColormap:
    if name not in PALETTES:
        available = ", ".join(sorted(PALETTES))
        raise KeyError(f"Unknown palette '{name}'. Available: {available}")
    return LinearSegmentedColormap.from_list(name, PALETTES[name], N=256)


class PowerTwoSlopeNorm(Normalize):
    """Symmetric diverging normalization with extra contrast near zero."""

    def __init__(self, vmin: float, vcenter: float, vmax: float, gamma: float = 0.45):
        super().__init__(vmin=vmin, vmax=vmax, clip=False)
        self.vcenter = vcenter
        self.gamma = gamma

    def __call__(self, value, clip=None):
        result, is_scalar = self.process_value(value)
        data = np.ma.getdata(result)
        neg = data <= self.vcenter
        negative_distance = np.zeros_like(data, dtype=float)
        positive_distance = np.zeros_like(data, dtype=float)

        if self.vcenter != self.vmin:
            negative_distance[neg] = np.clip(
                (self.vcenter - data[neg]) / (self.vcenter - self.vmin), 0, 1
            )
        if self.vmax != self.vcenter:
            positive_distance[~neg] = np.clip(
                (data[~neg] - self.vcenter) / (self.vmax - self.vcenter), 0, 1
            )

        mapped = np.where(
            neg,
            0.5 - 0.5 * negative_distance**self.gamma,
            0.5 + 0.5 * positive_distance**self.gamma,
        )
        output = np.ma.masked_array(mapped, mask=np.ma.getmask(result))
        return output.item() if is_scalar else output

    def inverse(self, value):
        result, is_scalar = self.process_value(value)
        position = np.ma.getdata(result).astype(float)
        neg = position <= 0.5
        output = np.zeros_like(position)

        negative_distance = np.clip(1 - 2 * position[neg], 0, 1) ** (
            1.0 / self.gamma
        )
        output[neg] = self.vcenter - negative_distance * (
            self.vcenter - self.vmin
        )

        positive_distance = np.clip(2 * position[~neg] - 1, 0, 1) ** (
            1.0 / self.gamma
        )
        output[~neg] = self.vcenter + positive_distance * (
            self.vmax - self.vcenter
        )

        output = np.ma.masked_array(output, mask=np.ma.getmask(result))
        return output.item() if is_scalar else output


def input_path(metric: str, antigen: str, mhc_class: str) -> Path:
    config = METRICS[metric]
    folder = config["folder"].format(mhc_class=mhc_class)
    filename = config["filename"].format(
        antigen=antigen, mhc_class=mhc_class
    )
    return BASE_DIR / f"Processing details{antigen}" / folder / filename


def read_matrix(path: Path, sheet_name: str, labels_map: dict[str, str]):
    """Read a square matrix and retain only explicitly supported tools."""
    workbook = openpyxl.load_workbook(path, data_only=True, read_only=True)
    try:
        if sheet_name not in workbook.sheetnames:
            raise KeyError(
                f"Sheet '{sheet_name}' not found in {path}. "
                f"Available sheets: {workbook.sheetnames}"
            )
        worksheet = workbook[sheet_name]
        rows = list(worksheet.iter_rows(values_only=True))
    finally:
        workbook.close()

    if len(rows) < 2:
        raise ValueError(f"Empty matrix in {path} / {sheet_name}")

    column_names = list(rows[0][1:])
    row_names = [row[0] for row in rows[1:]]
    values = [list(row[1:]) for row in rows[1:]]

    if row_names != column_names:
        raise ValueError(
            f"Rows and columns differ in {path} / {sheet_name}: "
            f"{row_names} vs {column_names}"
        )

    keep_indices = [
        index for index, tool_name in enumerate(column_names) if tool_name in labels_map
    ]
    if not keep_indices:
        raise ValueError(f"No supported tools found in {path} / {sheet_name}")

    retained = [[values[i][j] for j in keep_indices] for i in keep_indices]
    retained_names = [column_names[index] for index in keep_indices]
    matrix = np.asarray(retained, dtype=float)
    short_labels = [labels_map[name] for name in retained_names]
    return matrix, short_labels


def load_delta(metric: str, antigen: str, mhc_class: str):
    path = input_path(metric, antigen, mhc_class)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    config = METRICS[metric]
    labels_map = LABELS_BY_CLASS[mhc_class]
    rm_matrix, rm_labels = read_matrix(path, config["sheet_rm"], labels_map)
    all_matrix, all_labels = read_matrix(path, config["sheet_all"], labels_map)

    if rm_labels != all_labels:
        raise ValueError(f"Tool order differs between sheets in {path}")

    missing = np.isnan(rm_matrix) | np.isnan(all_matrix)
    delta = all_matrix - rm_matrix
    return delta, missing, rm_labels


def load_panel_data(metric: str, mhc_class: str):
    panels = []
    maximum = 0.0

    for antigen in ANTIGENS:
        delta, missing, labels = load_delta(metric, antigen, mhc_class)
        panels.append((antigen, delta, missing, labels))
        finite = delta[np.isfinite(delta)]
        if finite.size:
            maximum = max(maximum, float(np.max(np.abs(finite))))

    maximum = max(0.05, np.ceil(maximum / 0.05) * 0.05)
    norm = PowerTwoSlopeNorm(
        vmin=-maximum,
        vcenter=0.0,
        vmax=maximum,
        gamma=NORM_GAMMA,
    )
    return panels, norm


def draw_heatmap(
    ax, antigen, delta, missing, labels, norm, cmap, *, show_xlabels: bool
) -> None:
    count = len(labels)
    plot_data = delta.copy()
    np.fill_diagonal(plot_data, np.nan)
    masked = np.ma.masked_invalid(plot_data)

    ax.imshow(masked, cmap=cmap, norm=norm, aspect="equal", interpolation="none")

    for row in range(count):
        for column in range(count):
            if row != column and missing[row, column]:
                ax.add_patch(
                    Rectangle(
                        (column - 0.5, row - 0.5),
                        1,
                        1,
                        facecolor=MISSING_COLOR,
                        edgecolor="none",
                        zorder=1.5,
                    )
                )

    ax.set_xticks(range(count))
    if show_xlabels:
        ax.set_xticklabels(labels, rotation=45, ha="right")
    else:
        ax.set_xticklabels([])
    ax.set_yticks(range(count), labels=labels)
    ax.set_xticks(np.arange(-0.5, count, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, count, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0, zorder=3)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(which="major", bottom=False, left=False, pad=2)

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_title(antigen, fontweight="bold", pad=14)


def draw_panel(fig, outer_spec, letter, mhc_class, metric, title, cmap) -> bool:
    panel_grid = outer_spec.subgridspec(
        4,
        2,
        width_ratios=(1.0, 0.065),
        height_ratios=(0.38, 1.0, 0.30, 1.0),
        wspace=0.18,
        hspace=0.0,
    )
    heatmap_axes = [
        fig.add_subplot(panel_grid[1, 0]),
        fig.add_subplot(panel_grid[3, 0]),
    ]
    colorbar_axis = fig.add_subplot(panel_grid[1:, 1])

    panels, norm = load_panel_data(metric, mhc_class)
    for index, (axis, panel_data) in enumerate(zip(heatmap_axes, panels)):
        draw_heatmap(
            axis,
            *panel_data,
            norm,
            cmap,
            show_xlabels=index == len(heatmap_axes) - 1,
        )

    scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array([])
    colorbar = fig.colorbar(scalar_mappable, cax=colorbar_axis)
    colorbar.ax.tick_params(labelsize=FONT_SIZE, length=4, pad=6)
    colorbar.ax.axhline(0.0, color="#404040", linewidth=0.8, zorder=5)

    has_undefined = False
    for _, _, missing, _ in panels:
        off_diagonal_missing = missing.copy()
        np.fill_diagonal(off_diagonal_missing, False)
        if np.any(off_diagonal_missing):
            has_undefined = True
            break

    panel_position = outer_spec.get_position(fig)
    title_y = panel_position.y1 - 0.008
    fig.text(
        (panel_position.x0 + panel_position.x1) / 2,
        title_y,
        title,
        ha="center",
        va="top",
        fontsize=FONT_SIZE,
        fontweight="bold",
        linespacing=1.0,
    )
    fig.text(
        panel_position.x0 - 0.055,
        title_y,
        letter,
        ha="left",
        va="top",
        fontsize=FONT_SIZE,
        fontweight="bold",
    )
    return has_undefined


def build_figure():
    setup_style()
    cmap = make_colormap(ACTIVE_PALETTE)

    figure_size = (
        FIGURE_WIDTH_MM / MM_PER_INCH,
        FIGURE_HEIGHT_MM / MM_PER_INCH,
    )
    fig = plt.figure(figsize=figure_size, facecolor="white")
    outer_grid = fig.add_gridspec(
        2,
        2,
        left=0.115,
        right=0.905,
        bottom=0.125,
        top=0.955,
        wspace=0.70,
        hspace=0.32,
    )

    any_undefined = False
    for index, config in enumerate(PANEL_CONFIGS):
        row, column = divmod(index, 2)
        any_undefined |= draw_panel(
            fig, outer_grid[row, column], *config, cmap
        )

    if any_undefined:
        fig.legend(
            handles=[
                Patch(
                    facecolor=MISSING_COLOR,
                    edgecolor="none",
                    label="Undefined similarity",
                )
            ],
            loc="lower center",
            bbox_to_anchor=(0.5, 0.018),
            frameon=False,
            ncol=1,
            handlelength=1.3,
        )
    return fig


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / OUTPUT_FILE
    figure = build_figure()
    figure.savefig(output_path, dpi=DPI, facecolor="white")
    plt.close(figure)
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
