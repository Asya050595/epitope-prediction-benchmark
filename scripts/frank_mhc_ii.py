#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""English documentation for this module or helper is provided in the repository README and in the surrounding code."""

from project_paths import DATA_ROOT

import os
import re
import sys
import glob
import statistics
from collections import Counter

import pandas as pd
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
from allele_sets import get_allele_sets

# ═══════════════════════════════════════════════════════════════════════════
# Configuration.
# ═══════════════════════════════════════════════════════════════════════════

DEFAULT_BASE = f"{DATA_ROOT}"
BASE = os.path.abspath(os.path.expanduser(os.environ.get("RM_ROOT", DEFAULT_BASE)))

ANTIGEN_META = {
    "p24": {
        "xlsx": f"{BASE}/p24/p24.xlsx",
        "txt":  f"{BASE}/p24/p24.txt",
        "out_dir": f"{BASE}/p24/FRANK MHC II",
    },
    "pp65": {
        "xlsx": f"{BASE}/pp65/pp65.xlsx",
        "txt":  f"{BASE}/pp65/pp65.txt",
        "out_dir": f"{BASE}/pp65/FRANK MHC II",
    },
}

RAW_PATHS = {
    "IEDB_II_Consensus": {
        "p24":   f"{BASE}/p24/MHC II/IEDB II/IEDB_II_Consensus_p24_11-18_pep1.tsv",
        "pp65":  f"{BASE}/pp65/MHC II/IEDB II/IEDB_II_Consensus_pp65_11-18_pep1.tsv",
    },
    "IEDB_II_NetMHCIIpan_4.1_BA": {
        "p24":   f"{BASE}/p24/MHC II/IEDB II/IEDB_II_NetMHCIIpan_4.1_BA_p24_11-18_pep1.tsv",
        "pp65":  f"{BASE}/pp65/MHC II/IEDB II/IEDB_II_NetMHCIIpan_4.1_BA_pp65_11-18_pep1.tsv",
    },
    "IEDB_II_NetMHCIIpan_4.1_EL": {
        "p24":   f"{BASE}/p24/MHC II/IEDB II/IEDB_II_NetMHCIIpan_4.1_EL_p24_11-18_pep1.tsv",
        "pp65":  f"{BASE}/pp65/MHC II/IEDB II/IEDB_II_NetMHCIIpan_4.1_EL_pp65_11-18_pep1.tsv",
    },
    "NetMHCIIpan_4.1": {
        "p24": [
            f"{BASE}/p24/MHC II/NetMHCIIpan/NetMHCIIpan_p24_part_1_11-18.xls",
            f"{BASE}/p24/MHC II/NetMHCIIpan/NetMHCIIpan_p24_part_2_11-18.xls",
        ],
        "pp65": [
            f"{BASE}/pp65/MHC II/NetMHCIIpan/NetMHCIIpan_pp65_11-18.xls",
        ],
    },
    "NetMHC_II_2.3": {
        "p24":   os.path.join(BASE, "p24",   "MHC II", "NetMHC II", "NetMHCII_p24_part_*_*.html"),
        "pp65":  os.path.join(BASE, "pp65",  "MHC II", "NetMHC II", "NetMHCII_pp65_*.html"),
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# Input loading.
# ═══════════════════════════════════════════════════════════════════════════

class ToolLoadError(Exception):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    def __init__(self, kind: str, message: str):
        assert kind in ("input_file_missing", "parse_failed")
        self.kind = kind
        self.message = message
        super().__init__(f"[{kind}] {message}")


def _check_expected_alleles(found: set, expected: set, label: str,
                            fatal: bool = True, allow_extra: bool = False):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    missing = expected - found
    extra = found - expected
    if not missing and (allow_extra or not extra):
        return
    parts = []
    if missing:
        parts.append(f"Required input or value was not found{sorted(missing)}")
    if extra and not allow_extra:
        parts.append(f"Processing details{sorted(extra)}")
    message = f"{label}Allele status" + "; ".join(parts)
    if fatal:
        raise ToolLoadError("parse_failed", message)
    print(f"    [!] {message}Processing details")


def _check_no_cross_part_duplicates(alleles_per_file: dict, label: str):
    all_seen = []
    for alleles in alleles_per_file.values():
        all_seen.extend(alleles)
    counts = Counter(all_seen)
    dup = {a: [fn for fn, al in alleles_per_file.items() if a in al]
           for a, c in counts.items() if c > 1}
    if dup:
        raise ToolLoadError("parse_failed", f"{label}Allele status{dup}")


# ═══════════════════════════════════════════════════════════════════════════
# Normalization.
# ═══════════════════════════════════════════════════════════════════════════

_DASH_LOOKALIKES = {
    "\u2010": "-", "\u2011": "-", "\u2012": "-", "\u2013": "-",
    "\u2014": "-", "\u2212": "-", "\uFE63": "-", "\uFF0D": "-",
    "\u00AD": "",
}


def _normalize_dashes(s: str) -> str:
    fixed = s
    for bad, repl in _DASH_LOOKALIKES.items():
        if bad in fixed:
            fixed = fixed.replace(bad, repl)
    return fixed


def normalize_allele(raw) -> str:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if raw is None:
        raise ValueError("Required input or value was not found")
    s = str(raw).strip()
    if not s:
        raise ValueError("Allele status")
    s = _normalize_dashes(s)
    s = re.sub(r'^HLA-', '', s, flags=re.IGNORECASE)
    m = re.match(r'^([A-Z]+\d)[\*_]?(\d{2})[:_]?(\d{2})$', s, flags=re.IGNORECASE)
    if not m:
        raise ValueError(f"Allele status{raw!r}")
    gene, group, protein = m.groups()
    return f"{gene.upper()}*{group}:{protein}"


def norm_peptide(p) -> str:
    return str(p).strip().upper()


def mhc2_raw_to_standard(raw: str) -> str:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    raw = raw.strip()
    m = re.match(r'^((?:HLA-)?[A-Z]+\d*)[_*](\d{2})(\d{2})$', raw)
    if m:
        locus, g1, g2 = m.groups()
        locus = re.sub(r'^HLA-', '', locus)
        return f"{locus}*{g1}:{g2}"
    return raw


# ═══════════════════════════════════════════════════════════════════════════
# Allele handling.
# Validation.
# Allele handling.
# ═══════════════════════════════════════════════════════════════════════════

_IEDB_II_EXPECTED_RAW = {
    "p24": {
        "DRB1*01:01", "DRB1*03:01", "DRB1*04:01", "DRB1*04:05", "DRB1*07:01",
        "DRB1*09:01", "DRB1*11:01", "DRB1*13:02", "DRB1*15:01", "DRB3*01:01",
        "DRB4*01:01", "DRB5*01:01",
    },
    "pp65": {
        "DRB1*01:01", "DRB1*03:01", "DRB1*04:01", "DRB1*07:01", "DRB1*11:01",
        "DRB1*12:01", "DRB1*13:02", "DRB1*15:01", "DRB3*01:01", "DRB3*02:02",
    },
}

_NETMHCIIPAN_EXPECTED_RAW = {
    "p24": [
        "DRB1_0101", "DRB1_0301", "DRB1_0302", "DRB1_0401", "DRB1_0404", "DRB1_0405",
        "DRB1_0701", "DRB1_0801", "DRB1_0901", "DRB1_1001", "DRB1_1101", "DRB1_1301",
        "DRB1_1302", "DRB1_1303", "DRB1_1304", "DRB1_1401", "DRB1_1501", "DRB1_1502",
        "DRB3_0101", "DRB3_0301", "DRB3_0303", "DRB4_0101", "DRB5_0101",
    ],
    "pp65": [
        "DRB1_0101", "DRB1_0301", "DRB1_0302", "DRB1_0401", "DRB1_0701", "DRB1_1001",
        "DRB1_1101", "DRB1_1104", "DRB1_1201", "DRB1_1301", "DRB1_1302", "DRB1_1307",
        "DRB1_1401", "DRB1_1413", "DRB1_1501", "DRB1_1601", "DRB3_0101", "DRB3_0202",
        "DRB3_0301",
    ],
}

_NETMHCII23_EXPECTED_RAW = {
    "p24": [
        "DRB1_0101", "DRB1_0301", "DRB1_0401", "DRB1_0404", "DRB1_0405",
        "DRB1_0701", "DRB1_0801", "DRB1_0901", "DRB1_1001", "DRB1_1101",
        "DRB1_1301", "DRB1_1302", "DRB1_1501", "DRB3_0101", "DRB3_0301",
        "DRB4_0101", "DRB5_0101",
    ],
    "pp65": [
        "DRB1_0101", "DRB1_0301", "DRB1_0401", "DRB1_0701", "DRB1_1001",
        "DRB1_1101", "DRB1_1201", "DRB1_1301", "DRB1_1302", "DRB1_1501",
        "DRB3_0101", "DRB3_0202", "DRB3_0301",
    ],
}

EXPECTED_ALLELES = {
    "IEDB_II_Consensus":          {ag: set(a) for ag, a in _IEDB_II_EXPECTED_RAW.items()},
    "IEDB_II_NetMHCIIpan_4.1_BA": {ag: set(a) for ag, a in _IEDB_II_EXPECTED_RAW.items()},
    "IEDB_II_NetMHCIIpan_4.1_EL": {ag: set(a) for ag, a in _IEDB_II_EXPECTED_RAW.items()},
    "NetMHCIIpan_4.1": {
        ag: {mhc2_raw_to_standard(a) for a in alist} for ag, alist in _NETMHCIIPAN_EXPECTED_RAW.items()
    },
    "NetMHC_II_2.3": {
        ag: {mhc2_raw_to_standard(a) for a in alist} for ag, alist in _NETMHCII23_EXPECTED_RAW.items()
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# Input loading.
# ═══════════════════════════════════════════════════════════════════════════

def read_sequence(txt_path: str) -> str:
    with open(txt_path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    seq_lines = [
        ln.strip() for ln in lines
        if ln.strip() and not ln.lstrip().startswith((">", ";"))
    ]
    seq = "".join(seq_lines).upper()
    seq = re.sub(r"\s+", "", seq)
    return seq


# ═══════════════════════════════════════════════════════════════════════════
# Input loading.
#
# Allele handling.
# ═══════════════════════════════════════════════════════════════════════════

EXCEL_ALLELE_COL = "HLA allele.1"
EXCEL_PEPTIDE_COL = "HTL epitopes"


def load_validated_pairs(xlsx_path: str) -> list:
    try:
        df = pd.read_excel(xlsx_path, header=0)
    except Exception as exc:
        raise ValueError(f"Processing details{xlsx_path}: {exc}") from exc
    if EXCEL_ALLELE_COL not in df.columns or EXCEL_PEPTIDE_COL not in df.columns:
        raise ValueError(
            f"Processing details{xlsx_path}Required input or value was not found"
            f"'{EXCEL_ALLELE_COL}Processing details{EXCEL_PEPTIDE_COL}'. "
            f"Table status{list(df.columns)}"
        )
    sub = df[[EXCEL_ALLELE_COL, EXCEL_PEPTIDE_COL]].dropna()

    pairs, seen, errors = [], set(), []
    for idx, row in sub.iterrows():
        try:
            a = normalize_allele(row[EXCEL_ALLELE_COL])
        except ValueError as exc:
            errors.append(f"Table status{idx}: {exc}")
            continue
        p = norm_peptide(row[EXCEL_PEPTIDE_COL])
        if a and p and (a, p) not in seen:
            seen.add((a, p))
            pairs.append((a, p))

    if errors:
        sample = "; ".join(errors[:10])
        more = f" (+{len(errors) - 10}Processing details" if len(errors) > 10 else ""
        raise ValueError(
            f"Processing details{xlsx_path}Allele status"
            f"{len(errors)}Error"
            f"Table status"
            f"Processing details{sample}{more}"
        )
    if not pairs:
        raise ValueError(
            f"Processing details{xlsx_path}Required input or value was not found"
            f"Allele status"
        )
    return pairs


def select_evaluable_pairs(pairs: list, antigen: str) -> tuple[list, set, int]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    evaluable_alleles, _ = get_allele_sets("II", antigen)
    present_alleles = {allele for allele, _ in pairs}
    missing_alleles = evaluable_alleles - present_alleles
    if missing_alleles:
        raise ValueError(
            f"Processing details{antigen}Required input or value was not found"
            f"{', '.join(sorted(missing_alleles))}"
        )

    filtered = [(allele, peptide) for allele, peptide in pairs if allele in evaluable_alleles]
    if not filtered:
        raise ValueError(f"Processing details{antigen}Processing details")
    return filtered, evaluable_alleles, len(pairs) - len(filtered)


# ═══════════════════════════════════════════════════════════════════════════
# IEDB TSV (Consensus / NetMHCIIpan BA / NetMHCIIpan EL)
#
# Allele handling.
# Error handling.
# Validation.
# Validation.
# Implementation detail; see the repository documentation.
# ═══════════════════════════════════════════════════════════════════════════

def load_iedb_tsv_multi(tsv_path: str, score_cols: list, expected_alleles: set) -> dict:
    fname = os.path.basename(tsv_path)
    if not os.path.isfile(tsv_path):
        raise ToolLoadError("input_file_missing", f"Required input or value was not found{tsv_path}")

    try:
        df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    except (OSError, pd.errors.ParserError, UnicodeError) as exc:
        raise ToolLoadError("parse_failed", f"{fname}Processing details{exc}") from exc
    df.columns = [c.strip() for c in df.columns]
    required = {"allele", "peptide", *score_cols}
    missing_cols = required - set(df.columns)
    if missing_cols:
        raise ToolLoadError("parse_failed", f"{fname}Required input or value was not found{sorted(missing_cols)}")

    alleles_norm = []
    errors = []
    for idx, raw in df["allele"].items():
        try:
            alleles_norm.append(normalize_allele(raw))
        except ValueError as exc:
            errors.append(f"Table status{idx}: {exc}")
            alleles_norm.append(None)
    if errors:
        sample = "; ".join(errors[:10])
        more = f" (+{len(errors) - 10}Processing details" if len(errors) > 10 else ""
        raise ToolLoadError("parse_failed",
            f"{fname}Allele status{len(errors)}Table status"
            f"Processing details{sample}{more}")

    df = df.assign(allele_norm=alleles_norm)
    df["peptide_norm"] = df["peptide"].apply(norm_peptide)
    for col in score_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    found_alleles = set(df["allele_norm"].unique())
    if expected_alleles is not None:
        _check_expected_alleles(
            found_alleles, expected_alleles, fname, allow_extra=True
        )

    result = {}
    cols_arr = {col: df[col].to_numpy() for col in score_cols}
    alleles = df["allele_norm"].to_numpy()
    peptides = df["peptide_norm"].to_numpy()
    for i in range(len(df)):
        allele, peptide = alleles[i], peptides[i]
        d = result.setdefault(allele, {})
        fields = d.setdefault(peptide, {})
        for col in score_cols:
            v = cols_arr[col][i]
            if pd.isna(v):
                continue
            v = float(v)
            if col not in fields or v < fields[col]:
                fields[col] = v

    if not result:
        raise ToolLoadError("parse_failed", f"{fname}Required input or value was not found")
    return result


# ═══════════════════════════════════════════════════════════════════════════
# Implementation detail; see the repository documentation.
# ═══════════════════════════════════════════════════════════════════════════

def _empty_parse_stats() -> dict:
    return {
        "corrupted_rows": 0,
        "truncated_blocks": 0,
        "rows_per_allele": {},
        "nan_counts": {"rank_el": 0, "rank_ba": 0, "affinity_nM": 0},
        "peptide_rows": [],
    }


def parse_netmhciipan_file(path: str):
    fname = os.path.basename(path)
    if not os.path.isfile(path):
        raise ToolLoadError("input_file_missing", f"Required input or value was not found{path}")
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()
    except OSError as exc:
        raise ToolLoadError("parse_failed", f"{fname}Processing details{exc}") from exc

    header_idx = None
    for i, line in enumerate(lines[:10]):
        if line.split("\t", 1)[0].strip() == "Pos":
            header_idx = i
            break
    if header_idx is None or header_idx == 0:
        raise ToolLoadError("parse_failed", f"{fname}Required input or value was not found")

    allele_row = lines[header_idx - 1].rstrip("\n").split("\t")
    header_row = [c.strip() for c in lines[header_idx].rstrip("\n").split("\t")]

    block_starts = [i for i, v in enumerate(header_row) if v == "Core"]
    if not block_starts:
        raise ToolLoadError("parse_failed", f"{fname}Required input or value was not found")
    try:
        ave_idx = header_row.index("Ave")
    except ValueError:
        ave_idx = len(header_row)
    block_ends = block_starts[1:] + [ave_idx]

    blocks = []
    for start, end in zip(block_starts, block_ends):
        allele_raw = allele_row[start].strip() if start < len(allele_row) else ""
        if not allele_raw:
            raise ToolLoadError("parse_failed", f"{fname}Allele status{start}")
        sub_headers = header_row[start:end]
        try:
            rank_el_idx = start + sub_headers.index("Rank")
            rank_ba_idx = start + sub_headers.index("Rank_BA")
            afin_idx = start + sub_headers.index("nM")
        except ValueError:
            raise ToolLoadError("parse_failed",
                f"{fname}Validation status{allele_raw} ({sub_headers})")
        blocks.append((mhc2_raw_to_standard(allele_raw), rank_el_idx, rank_ba_idx, afin_idx))

    allele_names = [a for a, _, _, _ in blocks]
    if len(allele_names) != len(set(allele_names)):
        dup = [a for a, c in Counter(allele_names).items() if c > 1]
        raise ToolLoadError("parse_failed", f"{fname}Allele status{dup}")

    stats = _empty_parse_stats()
    stats["rows_per_allele"] = {a: 0 for a in allele_names}
    result = {a: {} for a in allele_names}

    def to_float(s):
        s = s.strip()
        if s in ("", "NA"):
            return None
        try:
            return float(s)
        except ValueError:
            return None

    for line in lines[header_idx + 1:]:
        raw = line.rstrip("\n")
        if not raw.strip():
            continue
        parts = raw.split("\t")
        try:
            int(parts[0])
        except (ValueError, IndexError):
            stats["corrupted_rows"] += 1
            continue
        if len(parts) < 2:
            stats["corrupted_rows"] += 1
            continue
        peptide = parts[1].strip().upper()
        stats["peptide_rows"].append(peptide)

        for allele_star, rank_el_idx, rank_ba_idx, afin_idx in blocks:
            needed_max = max(rank_el_idx, rank_ba_idx, afin_idx)
            if needed_max >= len(parts):
                stats["truncated_blocks"] += 1
                continue
            rank_el = to_float(parts[rank_el_idx])
            rank_ba = to_float(parts[rank_ba_idx])
            afin = to_float(parts[afin_idx])
            if rank_el is None:
                stats["nan_counts"]["rank_el"] += 1
            if rank_ba is None:
                stats["nan_counts"]["rank_ba"] += 1
            if afin is None:
                stats["nan_counts"]["affinity_nM"] += 1
            stats["rows_per_allele"][allele_star] += 1
            fields = result[allele_star].setdefault(peptide, {})
            if rank_el is not None and ("rank_el" not in fields or rank_el < fields["rank_el"]):
                fields["rank_el"] = rank_el
            if rank_ba is not None and ("rank_ba" not in fields or rank_ba < fields["rank_ba"]):
                fields["rank_ba"] = rank_ba
            if afin is not None and ("affinity_nM" not in fields or afin < fields["affinity_nM"]):
                fields["affinity_nM"] = afin

    if not stats["peptide_rows"]:
        raise ToolLoadError("parse_failed", f"{fname}Required input or value was not found")

    return result, stats


def _check_netmhciipan_integrity(file_stats: dict, tool_label: str):
    problems = []
    peptide_lists = {}
    for fname, stats in file_stats.items():
        if stats["corrupted_rows"] > 0:
            problems.append(f"{fname}: {stats['corrupted_rows']}Table status")
        if stats["truncated_blocks"] > 0:
            problems.append(f"{fname}: {stats['truncated_blocks']}Processing details")
        nan_total = sum(stats["nan_counts"].values())
        if nan_total > 0:
            problems.append(f"{fname}Processing details{stats['nan_counts']}")
        rpa = stats["rows_per_allele"]
        if len(set(rpa.values())) > 1:
            problems.append(f"{fname}Allele status{rpa}")
        peptide_lists[fname] = stats["peptide_rows"]

    if len(peptide_lists) > 1:
        names = list(peptide_lists.keys())
        base_name, base_list = names[0], peptide_lists[names[0]]
        for other_name in names[1:]:
            other_list = peptide_lists[other_name]
            if other_list != base_list:
                problems.append(
                    f"{base_name}Processing details{other_name}Peptide status"
                    f"Validation status{len(base_list)} vs {len(other_list)}Table status"
                )

    if problems:
        raise ToolLoadError("parse_failed", f"{tool_label}: {'; '.join(problems)}")


def load_tool_netmhciipan_4_1(antigen: str) -> dict:
    merged = {}
    file_stats = {}
    alleles_per_file = {}
    for p in RAW_PATHS["NetMHCIIpan_4.1"][antigen]:
        part, stats = parse_netmhciipan_file(p)
        fname = os.path.basename(p)
        file_stats[fname] = stats
        alleles_per_file[fname] = set(part.keys())
        for allele, d in part.items():
            merged.setdefault(allele, {}).update(d)

    _check_netmhciipan_integrity(file_stats, "NetMHCIIpan_4.1")
    _check_no_cross_part_duplicates(alleles_per_file, "NetMHCIIpan_4.1")

    if not merged:
        raise ToolLoadError("parse_failed", "Required input or value was not found")

    expected, _ = get_allele_sets("II", antigen)
    _check_expected_alleles(
        set(merged.keys()), expected, f"NetMHCIIpan_4.1 [{antigen}]", allow_extra=True
    )
    return merged


# ═══════════════════════════════════════════════════════════════════════════
# Metric calculation.
# Allele handling.
# Implementation detail; see the repository documentation.
# ═══════════════════════════════════════════════════════════════════════════

_NETMHCII23_LINE_RE = re.compile(
    r'^\s+'
    r'((?:DRB[0-9]|DP[AB][0-9]|DQ[AB][0-9])_\d+)'
    r'\s+(\d+)'
    r'\s+([A-Z]+)'
    r'\s+([A-Z]+)'
    r'\s+(-?\d+)'
    r'\s+(\S+)'
    r'\s+(\S+)'
    r'\s+(\S+)'
    r'\s+(\S+)'
    r'\s+(\S+)'
    r'(?:\s+(SB|WB))?'
    r'\s*$'
)


def parse_netmhcii23_html(path: str):
    result = {}
    n_matched = 0
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()
    except OSError as exc:
        raise ToolLoadError("parse_failed",
            f"{os.path.basename(path)}Processing details{exc}") from exc
    for line in lines:
        m = _NETMHCII23_LINE_RE.match(line)
        if not m:
            continue
        try:
            affinity = float(m.group(7))
            rank = float(m.group(8))
        except ValueError:
            continue
        allele_std = mhc2_raw_to_standard(m.group(1))
        peptide = m.group(3).upper()
        fields = result.setdefault(allele_std, {}).setdefault(peptide, {})
        if "rank" not in fields or rank < fields["rank"]:
            fields["rank"] = rank
        if "affinity" not in fields or affinity < fields["affinity"]:
            fields["affinity"] = affinity
        n_matched += 1
    return result, n_matched


def load_tool_netmhc_ii_2_3(antigen: str) -> dict:
    pattern = RAW_PATHS["NetMHC_II_2.3"][antigen]
    files = sorted(glob.glob(pattern))
    if not files:
        raise ToolLoadError("input_file_missing", f"Required input or value was not found{pattern}")
    merged = {}
    total_matched = 0
    for f in files:
        part, n_matched = parse_netmhcii23_html(f)
        total_matched += n_matched
        for allele, d in part.items():
            merged.setdefault(allele, {}).update(d)
    if total_matched == 0:
        raise ToolLoadError("parse_failed",
            f"Table status{len(files)}Processing details")

    expected, _ = get_allele_sets("II", antigen)
    _check_expected_alleles(
        set(merged.keys()), expected, f"NetMHC_II_2.3 [{antigen}]",
        fatal=False, allow_extra=True
    )
    return merged


# ═══════════════════════════════════════════════════════════════════════════
# Implementation detail; see the repository documentation.
# ═══════════════════════════════════════════════════════════════════════════

TOOL_SPECS = [
    {
        "name": "IEDB_II_Consensus",
        "loader": lambda antigen: load_iedb_tsv_multi(
            RAW_PATHS["IEDB_II_Consensus"][antigen],
            ["median binding percentile", "consensus percentile", "smm_align ic50", "nn_align ic50"],
            get_allele_sets("II", antigen)[0],
        ),
        "length_range": (11, 18),
        "scores": [
            {"key": "median binding percentile", "label": "Median_percentile", "higher_is_better": False},
            {"key": "consensus percentile",       "label": "Consensus_percentile", "higher_is_better": False},
            {"key": "smm_align ic50",             "label": "SMM_IC50", "higher_is_better": False},
            {"key": "nn_align ic50",              "label": "NN_IC50", "higher_is_better": False},
        ],
    },
    {
        "name": "IEDB_II_NetMHCIIpan_4.1_BA",
        "loader": lambda antigen: load_iedb_tsv_multi(
            RAW_PATHS["IEDB_II_NetMHCIIpan_4.1_BA"][antigen],
            ["median binding percentile", "netmhciipan_ba ic50"],
            get_allele_sets("II", antigen)[0],
        ),
        "length_range": (11, 18),
        "scores": [
            {"key": "median binding percentile", "label": "Median_percentile", "higher_is_better": False},
            {"key": "netmhciipan_ba ic50",        "label": "BA_IC50", "higher_is_better": False},
        ],
    },
    {
        "name": "IEDB_II_NetMHCIIpan_4.1_EL",
        "loader": lambda antigen: load_iedb_tsv_multi(
            RAW_PATHS["IEDB_II_NetMHCIIpan_4.1_EL"][antigen],
            ["median binding percentile"],
            get_allele_sets("II", antigen)[0],
        ),
        "length_range": (11, 18),
        "scores": [
            {"key": "median binding percentile", "label": "Median_percentile", "higher_is_better": False},
        ],
    },
    {
        "name": "NetMHCIIpan_4.1",
        "loader": load_tool_netmhciipan_4_1,
        "length_range": (11, 18),
        "scores": [
            {"key": "rank_el",     "label": "Rank_EL", "higher_is_better": False},
            {"key": "rank_ba",     "label": "Rank_BA", "higher_is_better": False},
            {"key": "affinity_nM", "label": "Affinity_nM", "higher_is_better": False},
        ],
    },
    {
        "name": "NetMHC_II_2.3",
        "loader": load_tool_netmhc_ii_2_3,
        "length_range": (11, 18),
        "scores": [
            {"key": "rank",     "label": "Rank_pct", "higher_is_better": False},
            {"key": "affinity", "label": "Affinity_nM", "higher_is_better": False},
        ],
    },
]

FATAL_STATUSES = {"input_file_missing", "parse_failed"}
REVIEW_STATUSES = {"epitope_score_missing", "incomplete_coverage",
                    "score_not_available_for_allele", "epitope_not_in_sequence"}
UNSUPPORTED_STATUSES = {"allele_not_supported", "unsupported_length"}


# ═══════════════════════════════════════════════════════════════════════════
# Allele handling.
# ═══════════════════════════════════════════════════════════════════════════

def evaluate_pair(peptide: str, allele: str, sequence: str,
                   length_range, higher_is_better: bool,
                   scalar_data: dict) -> dict:
    k = len(peptide)
    empty = {"status": None, "frank_pct": None, "n_better": None,
             "n_total": None, "epitope_score": None, "n_missing_windows": None}

    if length_range and not (length_range[0] <= k <= length_range[1]):
        empty["status"] = "unsupported_length"
        return empty

    if allele not in scalar_data:
        empty["status"] = "allele_not_supported"
        return empty
    score_lookup = scalar_data[allele]
    if not score_lookup:
        empty["status"] = "score_not_available_for_allele"
        return empty

    if peptide not in sequence:
        empty["status"] = "epitope_not_in_sequence"
        return empty

    total = len(sequence) - k + 1
    if total <= 0:
        empty["status"] = "epitope_not_in_sequence"
        return empty

    if peptide not in score_lookup:
        empty["status"] = "epitope_score_missing"
        empty["n_total"] = total
        return empty
    epitope_score = score_lookup[peptide]

    scores = []
    missing = 0
    for i in range(total):
        s = score_lookup.get(sequence[i:i + k])
        if s is None:
            missing += 1
        else:
            scores.append(s)

    if missing > 0:
        return {
            "status": "incomplete_coverage", "frank_pct": None, "n_better": None,
            "n_total": total, "epitope_score": epitope_score, "n_missing_windows": missing,
        }

    if higher_is_better:
        n_better = sum(1 for s in scores if s > epitope_score)
    else:
        n_better = sum(1 for s in scores if s < epitope_score)

    return {
        "status": "ok",
        "frank_pct": 100.0 * n_better / total,
        "n_better": n_better,
        "n_total": total,
        "epitope_score": epitope_score,
        "n_missing_windows": 0,
    }


def extract_scalar_dict(tool_data: dict, key: str) -> dict:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    out = {}
    for allele, pep_dict in tool_data.items():
        d = {}
        for peptide, fields in pep_dict.items():
            v = fields.get(key)
            if v is not None:
                d[peptide] = v
        out[allele] = d
    return out


# ═══════════════════════════════════════════════════════════════════════════
# Implementation detail; see the repository documentation.
# ═══════════════════════════════════════════════════════════════════════════

def make_sheet_name(name: str, used: set) -> str:
    safe = re.sub(r'[\\/*\[\]:?]', "_", name)[:31]
    base = safe
    i = 1
    while safe in used:
        suffix = f"_{i}"
        safe = base[:31 - len(suffix)] + suffix
        i += 1
    used.add(safe)
    return safe


# ═══════════════════════════════════════════════════════════════════════════
# Implementation detail; see the repository documentation.
# ═══════════════════════════════════════════════════════════════════════════

def process_antigen(antigen: str, meta: dict):
    print(f"\n{'=' * 70}")
    print(f"Processing details{antigen}  (MHC II)")
    print(f"{'=' * 70}")

    if not os.path.isfile(meta["xlsx"]):
        print(f"Required input or value was not found{meta['xlsx']}")
        return None
    if not os.path.isfile(meta["txt"]):
        print(f"Required input or value was not found{meta['txt']}")
        return None

    sequence = read_sequence(meta["txt"])
    if not sequence:
        print(f"Error{meta['txt']}")
        return None
    if not re.fullmatch(r"[A-Z]+", sequence):
        bad_chars = sorted(set(re.findall(r"[^A-Z]", sequence)))
        print(f"Error{bad_chars}")
        return None
    print(f"Processing details{len(sequence)}Processing details")

    try:
        all_validated_pairs = load_validated_pairs(meta["xlsx"])
        validated_pairs, evaluable_alleles, n_excluded = select_evaluable_pairs(
            all_validated_pairs, antigen
        )
    except ValueError as exc:
        print(f"Error{exc}")
        return None
    print(
        f"  Evaluable set: {len(evaluable_alleles)}Allele status"
        f"{len(validated_pairs)}Processing details{n_excluded}Processing details"
    )
    n_reference = len(validated_pairs)

    summary_rows = []
    detail_sheets = {}
    used_sheet_names = set()
    fatal_issues = []

    for spec in TOOL_SPECS:
        print(f"\n  --- {spec['name']} ---")
        tool_error = None
        try:
            tool_data = spec["loader"](antigen)
        except ToolLoadError as exc:
            tool_error = exc
            tool_data = {}
            fatal_issues.append((spec["name"], exc.kind, exc.message))
            print(f"Error{exc.kind}] {exc.message}")

        for score_spec in spec["scores"]:
            if tool_error is not None:
                rows = [{
                    "Allele": a, "Peptide": p, "Length": len(p),
                    "Epitope_score": None, "N_better": None, "N_total_peptides": None,
                    "N_missing_windows": None, "FRANK_%": None, "Status": tool_error.kind,
                } for a, p in validated_pairs]
            else:
                scalar_data = extract_scalar_dict(tool_data, score_spec["key"])
                rows = []
                for allele, peptide in validated_pairs:
                    res = evaluate_pair(
                        peptide, allele, sequence,
                        spec["length_range"], score_spec["higher_is_better"], scalar_data,
                    )
                    rows.append({
                        "Allele": allele, "Peptide": peptide, "Length": len(peptide),
                        "Epitope_score": res["epitope_score"],
                        "N_better": res["n_better"],
                        "N_total_peptides": res["n_total"],
                        "N_missing_windows": res["n_missing_windows"],
                        "FRANK_%": res["frank_pct"],
                        "Status": res["status"],
                    })

            df = pd.DataFrame(rows)
            sheet_name = make_sheet_name(f"{spec['name']}_{score_spec['label']}", used_sheet_names)
            detail_sheets[sheet_name] = df

            status_counts = df["Status"].value_counts().to_dict()
            n_evaluated = status_counts.get("ok", 0)
            n_input_failed = sum(v for k, v in status_counts.items() if k in FATAL_STATUSES)
            n_review = sum(v for k, v in status_counts.items() if k in REVIEW_STATUSES)
            n_unsupported = sum(v for k, v in status_counts.items() if k in UNSUPPORTED_STATUSES)

            ok = df[df["Status"] == "ok"]
            if n_evaluated > 0:
                mean_v = float(ok["FRANK_%"].mean())
                median_v = float(statistics.median(ok["FRANK_%"]))
                min_v = float(ok["FRANK_%"].min())
                max_v = float(ok["FRANK_%"].max())
            else:
                mean_v = median_v = min_v = max_v = None

            if n_input_failed > 0:
                status_flag = "INPUT_ERROR"
            elif n_review > 0:
                status_flag = "REVIEW_NEEDED"
            else:
                status_flag = "OK"

            skip_str = "; ".join(f"{k}: {v}" for k, v in sorted(status_counts.items()) if k != "ok")

            summary_rows.append({
                "Tool": spec["name"],
                "Score": score_spec["label"],
                "Sheet": sheet_name,
                "Reference_set": "Evaluable set",
                "N_evaluable_alleles": len(evaluable_alleles),
                "Status_flag": status_flag,
                "N_reference": n_reference,
                "N_evaluated": n_evaluated,
                "Coverage_%": round(100.0 * n_evaluated / n_reference, 1) if n_reference else None,
                "N_unsupported": n_unsupported,
                "N_review_needed": n_review,
                "N_input_failed": n_input_failed,
                "Mean_FRANK_%": mean_v,
                "Median_FRANK_%": median_v,
                "Min_FRANK_%": min_v,
                "Max_FRANK_%": max_v,
                "Skip_reasons": skip_str,
            })

            tag = f"[{status_flag}]"
            if n_evaluated:
                print(f"    {tag} [{score_spec['label']}Processing details{n_evaluated}/{n_reference}   "
                      f"Mean={mean_v:.2f}%   Median={median_v:.2f}%")
            else:
                print(f"    {tag} [{score_spec['label']}Processing details{n_reference}")
            if skip_str:
                print(f"      {skip_str}")

    summary_df = pd.DataFrame(summary_rows)
    return summary_df, detail_sheets, fatal_issues


# ═══════════════════════════════════════════════════════════════════════════
# Output generation.
# ═══════════════════════════════════════════════════════════════════════════

HEADER_FILL = PatternFill("solid", start_color="2F5496", end_color="2F5496")
SUMMARY_FILL = PatternFill("solid", start_color="1F3864", end_color="1F3864")
FATAL_FILL = PatternFill("solid", start_color="C00000", end_color="C00000")
REVIEW_ROW_FILL = PatternFill("solid", start_color="FFF2CC", end_color="FFF2CC")
ERROR_ROW_FILL = PatternFill("solid", start_color="F8CBAD", end_color="F8CBAD")
HEADER_FONT = Font(name="Arial", bold=True, color="FFFFFF", size=11)
BODY_FONT = Font(name="Arial", size=10)
ALIGN_CENTER = Alignment(horizontal="center", vertical="center", wrap_text=True)
ALIGN_LEFT = Alignment(horizontal="left", vertical="center", wrap_text=True)


def _autofit(ws):
    for col in ws.columns:
        max_len = 0
        col_letter = get_column_letter(col[0].column)
        for cell in col:
            try:
                val = str(cell.value) if cell.value is not None else ""
                max_len = max(max_len, len(val))
            except Exception:
                pass
        ws.column_dimensions[col_letter].width = min(max_len + 4, 60)


def _format_sheet(ws, fill):
    for cell in ws[1]:
        cell.fill = fill
        cell.font = HEADER_FONT
        cell.alignment = ALIGN_CENTER
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.font = BODY_FONT
            cell.alignment = ALIGN_LEFT
    ws.freeze_panes = "A2"
    _autofit(ws)


def _highlight_summary_rows(ws, summary_df: pd.DataFrame):
    try:
        list(summary_df.columns).index("Status_flag")
    except ValueError:
        return
    for r, flag in enumerate(summary_df["Status_flag"], start=2):
        fill = ERROR_ROW_FILL if flag == "INPUT_ERROR" else (REVIEW_ROW_FILL if flag == "REVIEW_NEEDED" else None)
        if fill is None:
            continue
        for cell in ws[r]:
            cell.fill = fill


def save_frank_excel(summary_df: pd.DataFrame, detail_sheets: dict, fatal_issues: list, out_path: str):
    sheets = {}
    if fatal_issues:
        sheets["FATAL_ERRORS"] = pd.DataFrame(
            [{"Tool": t, "Kind": k, "Message": m} for t, k, m in fatal_issues]
        )
    sheets["Summary"] = summary_df
    sheets.update(detail_sheets)

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, index=False, sheet_name=sheet_name)

    wb = openpyxl.load_workbook(out_path)
    for sheet_name in sheets.keys():
        if sheet_name == "FATAL_ERRORS":
            _format_sheet(wb[sheet_name], FATAL_FILL)
        elif sheet_name == "Summary":
            _format_sheet(wb[sheet_name], SUMMARY_FILL)
            _highlight_summary_rows(wb[sheet_name], summary_df)
        else:
            _format_sheet(wb[sheet_name], HEADER_FILL)
    wb.save(out_path)


# ═══════════════════════════════════════════════════════════════════════════
# main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    print("Processing details")
    antigen_failed = []
    all_fatal_issues = []

    for antigen, meta in ANTIGEN_META.items():
        result = process_antigen(antigen, meta)
        if result is None:
            antigen_failed.append(antigen)
            continue
        summary_df, detail_sheets, fatal_issues = result
        for tool, kind, message in fatal_issues:
            all_fatal_issues.append((antigen, tool, kind, message))

        out_dir = meta["out_dir"]
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f"FRANK_{antigen}_MHC_II.xlsx")
        save_frank_excel(summary_df, detail_sheets, fatal_issues, out_path)
        print(f"Saved output{out_path}")
        if fatal_issues:
            print(f"Warning{len(fatal_issues)}Loaded input"
                  f"Processing details")

    print(f"\n{'=' * 70}")
    print("Processing details")
    print(f"{'=' * 70}")
    if antigen_failed:
        print(f"Required input or value was not found"
              f"Allele status{', '.join(antigen_failed)}")
    if all_fatal_issues:
        print(f"Error{len(all_fatal_issues)}):")
        for antigen, tool, kind, message in all_fatal_issues:
            print(f"  [{antigen}] {tool} ({kind}): {message}")
        print("Warning"
              "Processing details"
              "Processing details"
              "Allele status"
              "Table status"
              "Processing details")

    if antigen_failed or all_fatal_issues:
        sys.exit(1)
    else:
        print("Error")


if __name__ == "__main__":
    main()
