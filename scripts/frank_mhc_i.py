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
from bs4 import BeautifulSoup
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
        "out_dir": f"{BASE}/p24/FRANK MHC I",
    },
    "pp65": {
        "xlsx": f"{BASE}/pp65/pp65.xlsx",
        "txt":  f"{BASE}/pp65/pp65.txt",
        "out_dir": f"{BASE}/pp65/FRANK MHC I",
    },
}

RAW_PATHS = {
    "IEDB_I_Consensus": {
        "p24":   f"{BASE}/p24/MHC I/IEDB I/IEDB_I_Consensus_p24_8-14.tsv",
        "pp65":  f"{BASE}/pp65/MHC I/IEDB I/IEDB_I_Consensus_pp65_8-14.tsv",
    },
    "IEDB_I_NetMHCpan_4.1_BA": {
        "p24":   f"{BASE}/p24/MHC I/IEDB I/IEDB_I_NetMHCpan_4.1_BA_p24_8-14.tsv",
        "pp65":  f"{BASE}/pp65/MHC I/IEDB I/IEDB_I_NetMHCpan_4.1_BA_pp65_8-14.tsv",
    },
    "IEDB_I_NetMHCpan_4.1_EL": {
        "p24":   f"{BASE}/p24/MHC I/IEDB I/IEDB_I_NetMHCpan_4.1_EL_p24_8-14.tsv",
        "pp65":  f"{BASE}/pp65/MHC I/IEDB I/IEDB_I_NetMHCpan_4.1_EL_pp65_8-14.tsv",
    },
    "NetMHCpan_4.1": {
        "p24": [
            f"{BASE}/p24/MHC I/NetMHCpan/NetMHCpan_p24_part_1_8-14.xls",
            f"{BASE}/p24/MHC I/NetMHCpan/NetMHCpan_p24_part_2_8-14.xls",
            f"{BASE}/p24/MHC I/NetMHCpan/NetMHCpan_p24_part_3_8-14.xls",
            f"{BASE}/p24/MHC I/NetMHCpan/NetMHCpan_p24_part_4_8-14.xls",
        ],
        "pp65": [
            f"{BASE}/pp65/MHC I/NetMHCpan/NetMHCpan_pp65_part_1_8-14.xls",
            f"{BASE}/pp65/MHC I/NetMHCpan/NetMHCpan_pp65_part_2_8-14.xls",
        ],
    },
    "NetMHC_4.0": {
        "p24":   os.path.join(BASE, "p24",   "MHC I", "NetMHC", "NetMHC_p24*8-14.xls"),
        "pp65":  os.path.join(BASE, "pp65",  "MHC I", "NetMHC", "NetMHC_pp65*8-14.xls"),
    },
    "NetCTL_1.2": {
        "p24":   os.path.join(BASE, "p24",   "MHC I", "NetCTL"),
        "pp65":  os.path.join(BASE, "pp65",  "MHC I", "NetCTL"),
    },
}

# Allele handling.
ALLELE_TO_SUPERTYPE = {
    "A*01:01": "A1",  "A*30:01": "A1",  "A*30:02": "A1",
    "A*02:01": "A2",  "A*02:07": "A2",  "A*68:02": "A2",
    "A*11:01": "A3",  "A*11:03": "A3",  "A*32:01": "A3",
    "A*33:03": "A3",  "A*68:01": "A3",  "A*74:01": "A3",
    "A*24:02": "A24", "A*24:07": "A24",
    "A*25:01": "A26", "A*26:01": "A26", "A*26:02": "A26",
    "A*26:03": "A26", "A*26:38": "A26",
    "B*07:02": "B7",  "B*35:01": "B7",  "B*35:02": "B7",
    "B*35:03": "B7",  "B*35:05": "B7",  "B*35:08": "B7",
    "B*35:11": "B7",  "B*42:01": "B7",  "B*42:02": "B7",
    "B*51:01": "B7",  "B*53:01": "B7",  "B*67:01": "B7",
    "B*08:01": "B8",
    "B*14:01": "B27", "B*14:02": "B27", "B*14:03": "B27",
    "B*27:05": "B27", "B*48:01": "B27",
    "B*39:01": "B39", "B*39:10": "B39",
    "B*13:02": "B44", "B*40:01": "B44", "B*40:02": "B44",
    "B*40:06": "B44", "B*44:02": "B44", "B*44:03": "B44",
    "B*44:15": "B44", "B*45:01": "B44", "B*50:01": "B44",
    "B*57:01": "B58", "B*57:02": "B58", "B*57:03": "B58",
    "B*58:01": "B58", "B*58:02": "B58",
    "B*15:01": "B62", "B*15:02": "B62", "B*15:03": "B62",
    "B*15:10": "B62", "B*15:16": "B62", "B*15:17": "B62",
    "B*15:24": "B62", "B*15:40": "B62", "B*52:01": "B62",
}

NETCTL_THRESHOLD = 0.75
NETCTL_EXPECTED_FILE_COUNTS = {"p24": 12, "pp65": 9}
NETCTL_EXPECTED_SUPERTYPES = {
    "p24":   {"A1", "A2", "A3", "A24", "A26", "B7", "B8", "B27", "B39", "B44", "B58", "B62"},
    "pp65":  {"A1", "A2", "A3", "A24", "A26", "B7", "B44", "B58", "B62"},
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
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
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

def strip_hla_prefix(allele) -> str:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    a = str(allele).strip()
    return a[4:] if a.upper().startswith("HLA-") else a


def nostar_to_star(allele_nostar: str) -> str:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    a = allele_nostar.strip()
    if len(a) >= 2 and a[0].isalpha():
        return f"{a[0]}*{a[1:]}"
    return a


def netmhc40_allele_to_star(raw: str) -> str:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    core = raw.strip()
    if core.upper().startswith("HLA-"):
        core = core[4:]
    m = re.match(r'^([A-Za-z]+)(\d{4,})$', core)
    if not m:
        return core
    gene, digits = m.group(1), m.group(2)
    return f"{gene}*{digits[:2]}:{digits[2:]}"


def norm_peptide(p) -> str:
    return str(p).strip().upper()


# ═══════════════════════════════════════════════════════════════════════════
# Allele handling.
# Validation.
# Allele handling.
# ═══════════════════════════════════════════════════════════════════════════

_IEDB_I_EXPECTED_RAW = {
    "p24": {
        "A*01:01", "A*02:01", "A*11:01", "A*24:02", "A*26:01", "A*30:02", "A*68:01", "A*68:02",
        "B*07:02", "B*08:01", "B*15:01", "B*35:01", "B*40:01", "B*44:02", "B*44:03",
        "B*51:01", "B*53:01", "B*57:01", "B*58:01",
    },
    "pp65": {
        "A*01:01", "A*02:01", "A*11:01", "A*24:02", "A*26:01", "A*30:01", "A*32:01", "A*68:01",
        "B*07:02", "B*15:01", "B*35:01", "B*40:01", "B*44:02", "B*44:03",
        "B*51:01", "B*53:01", "B*57:01", "B*58:01",
    },
}

_NETMHCPAN_EXPECTED_RAW = {
    "p24": [
        "A01:01", "A02:01", "A02:07", "A11:01", "A11:03",
        "A24:02", "A25:01", "A26:01", "A26:02", "A26:03",
        "A26:38", "A30:02", "A33:03", "A68:01", "A68:02",
        "A74:01", "B07:02", "B08:01", "B13:02", "B14:01",
        "B14:02", "B14:03", "B15:01", "B15:02", "B15:03",
        "B15:10", "B15:16", "B15:17", "B15:24", "B15:40",
        "B27:05", "B35:01", "B35:02", "B35:03", "B35:05",
        "B35:08", "B39:01", "B39:10", "B40:01", "B40:02",
        "B40:06", "B42:01", "B42:02", "B44:02", "B44:03",
        "B44:15", "B45:01", "B48:01", "B50:01", "B51:01",
        "B53:01", "B57:01", "B57:02", "B57:03", "B58:01",
        "B58:02", "B67:01", "C01:02", "C03:03", "C03:04",
        "C04:01", "C06:02", "C07:01", "C08:02", "C18:01",
        "E01:01", "E01:03",
    ],
    "pp65": [
        "A01:01", "A02:01", "A02:07", "A11:01", "A24:02",
        "A24:07", "A26:01", "A30:01", "A32:01", "A68:01",
        "B07:02", "B15:01", "B35:01", "B35:02", "B35:03",
        "B35:08", "B35:11", "B40:01", "B40:02", "B40:06",
        "B42:01", "B42:02", "B44:02", "B44:03", "B51:01",
        "B52:01", "B53:01", "B57:01", "B58:01", "C01:02",
        "C04:01", "C08:01", "C12:02", "C15:02",
    ],
}

_NETMHC40_EXPECTED_RAW = {
    "p24": {
        "HLA-A0101", "HLA-A0201", "HLA-A0207", "HLA-A1101", "HLA-A2402",
        "HLA-A2501", "HLA-A2601", "HLA-A2602", "HLA-A2603", "HLA-A3002",
        "HLA-A6801", "HLA-A6802", "HLA-B0702", "HLA-B0801", "HLA-B1402",
        "HLA-B1501", "HLA-B1502", "HLA-B1503", "HLA-B1517", "HLA-B2705",
        "HLA-B3501", "HLA-B3503", "HLA-B3901", "HLA-B4001", "HLA-B4002",
        "HLA-B4201", "HLA-B4402", "HLA-B4403", "HLA-B4501", "HLA-B4801",
        "HLA-B5101", "HLA-B5301", "HLA-B5701", "HLA-B5801", "HLA-B5802",
        "HLA-C0303", "HLA-C0401", "HLA-C0602", "HLA-C0701", "HLA-C0802",
        "HLA-E0101",
    },
    "pp65": {
        "HLA-A0101", "HLA-A0201", "HLA-A0207", "HLA-A1101", "HLA-A2402",
        "HLA-A2601", "HLA-A3001", "HLA-A3201", "HLA-A6801", "HLA-B0702",
        "HLA-B1501", "HLA-B3501", "HLA-B3503", "HLA-B4001", "HLA-B4002",
        "HLA-B4201", "HLA-B4402", "HLA-B4403", "HLA-B5101", "HLA-B5301",
        "HLA-B5701", "HLA-B5801", "HLA-C0401", "HLA-C1502",
    },
}

EXPECTED_ALLELES = {
    "IEDB_I_Consensus":         {ag: set(a) for ag, a in _IEDB_I_EXPECTED_RAW.items()},
    "IEDB_I_NetMHCpan_4.1_BA":  {ag: set(a) for ag, a in _IEDB_I_EXPECTED_RAW.items()},
    "IEDB_I_NetMHCpan_4.1_EL":  {ag: set(a) for ag, a in _IEDB_I_EXPECTED_RAW.items()},
    "NetMHCpan_4.1": {
        ag: {nostar_to_star(a) for a in alist} for ag, alist in _NETMHCPAN_EXPECTED_RAW.items()
    },
    "NetMHC_4.0": {
        ag: {netmhc40_allele_to_star(a) for a in aset} for ag, aset in _NETMHC40_EXPECTED_RAW.items()
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
# ═══════════════════════════════════════════════════════════════════════════

def load_validated_pairs(xlsx_path: str) -> list:
    try:
        df = pd.read_excel(xlsx_path)
    except Exception as exc:
        raise ValueError(f"Processing details{xlsx_path}: {exc}") from exc
    required = {"HLA allele", "CTL epitopes"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Processing details{xlsx_path}Required input or value was not found{missing}")
    pairs, seen = [], set()
    for _, row in df[["HLA allele", "CTL epitopes"]].dropna().iterrows():
        a = strip_hla_prefix(row["HLA allele"])
        p = norm_peptide(row["CTL epitopes"])
        if a and p and (a, p) not in seen:
            seen.add((a, p))
            pairs.append((a, p))
    if not pairs:
        raise ValueError(
            f"Processing details{xlsx_path}Required input or value was not found"
            f"Allele status"
        )
    return pairs


def select_evaluable_pairs(pairs: list, antigen: str) -> tuple[list, set, int]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    evaluable_alleles, _ = get_allele_sets("I", antigen)
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
# IEDB TSV (Consensus / NetMHCpan BA / NetMHCpan EL)
#
# Validation.
# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
# Allele handling.
# Implementation detail; see the repository documentation.
# Allele handling.
# ═══════════════════════════════════════════════════════════════════════════

def load_iedb_tsv_multi(tsv_path: str, score_cols: list, length_range: tuple,
                         allowed_nan_cols: set, expected_alleles: set) -> dict:
    fname = os.path.basename(tsv_path)
    if not os.path.isfile(tsv_path):
        raise ToolLoadError("input_file_missing", f"Required input or value was not found{tsv_path}")

    try:
        df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    except (OSError, pd.errors.ParserError, UnicodeError) as exc:
        raise ToolLoadError("parse_failed", f"{fname}Processing details{exc}") from exc
    df.columns = [c.strip() for c in df.columns]
    required = {"allele", "peptide", "peptide length", *score_cols}
    missing_cols = required - set(df.columns)
    if missing_cols:
        raise ToolLoadError("parse_failed", f"{fname}Required input or value was not found{sorted(missing_cols)}")

    lengths = pd.to_numeric(df["peptide length"], errors="coerce")
    n_nan_len = lengths.isna().sum()
    if n_nan_len > 0:
        raise ToolLoadError("parse_failed",
            f"{fname}Processing details{n_nan_len}Processing details")

    actual_lens = df["peptide"].astype(str).str.strip().str.len()
    mismatch = int((actual_lens != lengths).sum())
    if mismatch > 0:
        raise ToolLoadError("parse_failed",
            f"{fname}Validation status"
            f"{mismatch}Table status")

    out_of_range = int(((lengths < length_range[0]) | (lengths > length_range[1])).sum())
    if out_of_range > 0:
        raise ToolLoadError("parse_failed",
            f"{fname}: {out_of_range}Peptide status"
            f"{length_range[0]}-{length_range[1]}")

    # Implementation detail; see the repository documentation.
    for col in score_cols:
        converted = pd.to_numeric(df[col], errors="coerce")
        n_missing = int(converted.isna().sum())
        if n_missing > 0 and col not in allowed_nan_cols:
            raise ToolLoadError("parse_failed",
                f"{fname}Table status{col}Processing details{n_missing}Processing details"
                f"Processing details{sorted(allowed_nan_cols)})")

    df["allele_norm"] = df["allele"].apply(strip_hla_prefix)
    df["peptide_norm"] = df["peptide"].apply(norm_peptide)
    for col in score_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Allele handling.
    peptides_by_allele = {a: set(g["peptide_norm"]) for a, g in df.groupby("allele_norm")}
    if len(peptides_by_allele) > 1:
        ref_allele = next(iter(peptides_by_allele))
        ref_set = peptides_by_allele[ref_allele]
        differing = [a for a in peptides_by_allele if peptides_by_allele[a] != ref_set]
        if differing:
            raise ToolLoadError("parse_failed",
                f"{fname}Allele status{ref_allele}): "
                f"{sorted(differing)[:5]}")

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

_NETMHCPAN_BLOCK_FIELDS = ["core", "icore", "EL-score", "EL_Rank", "BA-score", "BA_Rank"]
_NETMHCPAN_FIXED_COLS = 3
_NETMHCPAN_TAIL_COLS = 2


def _validate_netmhcpan_header(lines: list, fname: str):
    if len(lines) < 3:
        raise ToolLoadError("parse_failed", f"{fname}Table status")

    allele_row = lines[0].rstrip("\n").split("\t")
    raw_alleles = [tok.strip() for tok in allele_row if tok.strip()]
    col_starts = [i for i, tok in enumerate(allele_row) if tok.strip()]
    if not raw_alleles:
        raise ToolLoadError("parse_failed", f"{fname}Allele status")

    alleles_norm = [nostar_to_star(a.replace("HLA-", "").replace("*", "")) for a in raw_alleles]
    if len(alleles_norm) != len(set(alleles_norm)):
        dup = [a for a, c in Counter(alleles_norm).items() if c > 1]
        raise ToolLoadError("parse_failed", f"{fname}Allele status{dup}")

    n = len(raw_alleles)
    col_row = lines[1].rstrip("\n").split("\t")
    expected_total = _NETMHCPAN_FIXED_COLS + 6 * n + _NETMHCPAN_TAIL_COLS
    if len(col_row) != expected_total:
        raise ToolLoadError("parse_failed",
            f"{fname}Validation status{expected_total}Table status"
            f"({_NETMHCPAN_FIXED_COLS}+6×{n}+{_NETMHCPAN_TAIL_COLS}Processing details{len(col_row)}")

    for allele_raw, cs in zip(raw_alleles, col_starts):
        for offset, exp_name in enumerate(_NETMHCPAN_BLOCK_FIELDS):
            idx = cs + offset
            actual = col_row[idx].strip() if idx < len(col_row) else "Processing details"
            if actual != exp_name:
                raise ToolLoadError("parse_failed",
                    f"{fname}Allele status{allele_raw!r} (+{offset}Validation status"
                    f"'{exp_name}Processing details{actual}'")

    return alleles_norm, col_starts


def parse_netmhcpan_wide(path: str, length_range: tuple):
    fname = os.path.basename(path)
    if not os.path.isfile(path):
        raise ToolLoadError("input_file_missing", f"Required input or value was not found{path}")
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()
    except OSError as exc:
        raise ToolLoadError("parse_failed", f"{fname}Processing details{exc}") from exc

    alleles_norm, col_starts = _validate_netmhcpan_header(lines, fname)
    blocks = [(a, cs + 3, cs + 5) for a, cs in zip(alleles_norm, col_starts)]  # (allele, el_idx, ba_idx)

    result = {a: {} for a, _, _ in blocks}
    signature = []
    seen_sig = set()
    errors = []
    min_len = col_starts[-1] + 6 if col_starts else 0

    for lineno, line in enumerate(lines[2:], start=3):
        raw_line = line.rstrip("\n")
        if not raw_line.strip():
            continue
        cols = raw_line.split("\t")
        if len(cols) < 3:
            errors.append(f"{fname}:{lineno}Table status")
            continue

        peptide_raw = cols[1].strip()
        if not peptide_raw or not peptide_raw.replace("-", "").isalpha():
            errors.append(f"{fname}:{lineno}Processing details{peptide_raw!r}")
            continue
        peptide = peptide_raw.upper()
        if not (length_range[0] <= len(peptide) <= length_range[1]):
            errors.append(f"{fname}:{lineno}Peptide status{len(peptide)}Processing details{length_range}")
            continue

        sig_key = (cols[0].strip(), peptide, cols[2].strip())
        if sig_key in seen_sig:
            errors.append(f"{fname}:{lineno}Table status{sig_key}")
            continue
        seen_sig.add(sig_key)
        signature.append(sig_key)

        if len(cols) < min_len:
            errors.append(f"{fname}:{lineno}Validation status{len(cols)} < {min_len})")
            continue

        for allele_star, el_idx, ba_idx in blocks:
            try:
                el_val = float(cols[el_idx])
                ba_val = float(cols[ba_idx])
            except (ValueError, IndexError):
                errors.append(f"{fname}:{lineno}Allele status{allele_star}Processing details")
                continue
            if not (0.0 <= el_val <= 100.0) or not (0.0 <= ba_val <= 100.0):
                errors.append(
                    f"{fname}:{lineno}Allele status{allele_star}: "
                    f"Processing details{el_val}, BA={ba_val}"
                )
                continue
            fields = result[allele_star].setdefault(peptide, {})
            if "rank_el" not in fields or el_val < fields["rank_el"]:
                fields["rank_el"] = el_val
            if "rank_ba" not in fields or ba_val < fields["rank_ba"]:
                fields["rank_ba"] = ba_val

    if errors:
        sample = "; ".join(errors[:5])
        more = f" (+{len(errors) - 5}Processing details" if len(errors) > 5 else ""
        raise ToolLoadError("parse_failed", f"{fname}: {len(errors)}Table status{sample}{more}")
    if not signature:
        raise ToolLoadError("parse_failed", f"{fname}Required input or value was not found")

    return result, signature


def load_tool_netmhcpan_4_1(antigen: str, length_range: tuple = (8, 14)) -> dict:
    merged = {}
    alleles_per_file = {}
    reference_sig, reference_fname = None, None
    for p in RAW_PATHS["NetMHCpan_4.1"][antigen]:
        part, sig = parse_netmhcpan_wide(p, length_range)
        fname = os.path.basename(p)
        alleles_per_file[fname] = set(part.keys())
        if reference_sig is None:
            reference_sig, reference_fname = sig, fname
        elif sig != reference_sig:
            raise ToolLoadError("parse_failed",
                f"{fname}Validation status{reference_fname} "
                f"({len(sig)} vs {len(reference_sig)}Table status"
                f"Processing details")
        for allele, d in part.items():
            merged.setdefault(allele, {}).update(d)

    _check_no_cross_part_duplicates(alleles_per_file, "NetMHCpan_4.1")
    expected, _ = get_allele_sets("I", antigen)
    _check_expected_alleles(
        set(merged.keys()), expected, f"NetMHCpan_4.1 [{antigen}]", allow_extra=True
    )
    return merged


# ═══════════════════════════════════════════════════════════════════════════
# Implementation detail; see the repository documentation.
#
# Allele handling.
# Validation.
# Allele handling.
# Allele handling.
# Allele handling.
# ═══════════════════════════════════════════════════════════════════════════

def parse_netmhc40_wide(path: str) -> dict:
    fname = os.path.basename(path)
    if not os.path.isfile(path):
        raise ToolLoadError("input_file_missing", f"Required input or value was not found{path}")
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()
    except OSError as exc:
        raise ToolLoadError("parse_failed", f"{fname}Processing details{exc}") from exc
    if len(lines) < 3:
        raise ToolLoadError("parse_failed", f"{fname}Table status")

    row0 = lines[0].rstrip("\n").split("\t")
    row1 = lines[1].rstrip("\n").split("\t")

    blocks = []
    for i, v in enumerate(row0):
        if isinstance(v, str) and v.startswith("HLA-"):
            nm_col, rank_col = i, i + 1
            h_nm = row1[nm_col] if nm_col < len(row1) else ""
            h_rank = row1[rank_col] if rank_col < len(row1) else ""
            if h_nm != "nM" or h_rank != "Rank":
                print(f"    [!] {fname}Allele status{v} (col {i}Processing details"
                      f"('{h_nm}'/'{h_rank}Processing details")
                continue
            blocks.append((netmhc40_allele_to_star(v), nm_col, rank_col))

    if not blocks:
        raise ToolLoadError("parse_failed", f"{fname}Required input or value was not found")

    result = {a: {} for a, _, _ in blocks}
    pep_re = re.compile(r"^[A-Z]{5,}$")
    n_rows = 0
    for line in lines[2:]:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 3:
            continue
        peptide = cols[1].strip().upper()
        if not peptide or not pep_re.match(peptide):
            continue
        n_rows += 1
        for allele_star, nm_col, rank_col in blocks:
            if rank_col >= len(cols):
                continue
            try:
                nm_val = float(cols[nm_col])
                rank_val = float(cols[rank_col])
            except ValueError:
                continue
            fields = result[allele_star].setdefault(peptide, {})
            if "affinity_nM" not in fields or nm_val < fields["affinity_nM"]:
                fields["affinity_nM"] = nm_val
            if "rank_pct" not in fields or rank_val < fields["rank_pct"]:
                fields["rank_pct"] = rank_val

    if n_rows == 0:
        raise ToolLoadError("parse_failed", f"{fname}Required input or value was not found")

    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    result = {a: d for a, d in result.items() if d}
    return result


def load_tool_netmhc_4_0(antigen: str) -> dict:
    pattern = RAW_PATHS["NetMHC_4.0"][antigen]
    files = sorted(glob.glob(pattern))
    if not files:
        raise ToolLoadError("input_file_missing", f"Required input or value was not found{pattern}")
    merged = {}
    for p in files:
        part = parse_netmhc40_wide(p)
        for allele, d in part.items():
            merged.setdefault(allele, {}).update(d)

    expected, _ = get_allele_sets("I", antigen)
    _check_expected_alleles(
        set(merged.keys()), expected, f"NetMHC_4.0 [{antigen}]", allow_extra=True
    )
    return merged


# ═══════════════════════════════════════════════════════════════════════════
# Implementation detail; see the repository documentation.
# Validation.
# ═══════════════════════════════════════════════════════════════════════════

_NETCTL_LINE_RE = re.compile(
    r'^\s*(\d+)\s+ID\s+(\S+)\s+pep\s+(\S+)'
    r'\s+aff\s+([\d.\-]+)\s+aff_rescale\s+([\d.\-]+)'
    r'\s+cle\s+([\d.\-]+)\s+tap\s+([\d.\-]+)'
    r'\s+COMB\s+([\d.\-]+)'
    r'(\s+<-E)?',
    re.MULTILINE,
)


def _extract_netctl_supertype(raw_html: str, filepath: str) -> str:
    basename = os.path.basename(filepath)
    m_html = re.search(r'using\s+MHC\s+supertype\s+(\w+)', raw_html, re.IGNORECASE)
    m_fname = re.search(
        r'NetCTL[_\-]1[._]2[_\-]([A-Za-z0-9]+)(?:\(\d+\))?\.html?$', basename, re.IGNORECASE,
    )
    if m_html and m_fname:
        st_html, st_fname = m_html.group(1).upper(), m_fname.group(1).upper()
        if st_html != st_fname:
            raise ToolLoadError("parse_failed",
                f"{basename}Processing details{st_html}Validation status"
                f"Processing details{st_fname}')")
        return st_html
    if m_html and not m_fname:
        raise ToolLoadError("parse_failed",
            f"{basename}Processing details")
    if m_fname:
        return m_fname.group(1).upper()
    raise ToolLoadError("parse_failed", f"{basename}Processing details")


def parse_netctl_html(path: str):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    basename = os.path.basename(path)
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            raw = fh.read()
    except OSError as exc:
        raise ToolLoadError("parse_failed", f"{basename}Processing details{exc}") from exc

    supertype = _extract_netctl_supertype(raw, path)

    tm = re.search(r'Threshold\s+([0-9.]+)', raw, re.IGNORECASE)
    if not tm:
        raise ToolLoadError("parse_failed", f"{basename}Required input or value was not found")
    html_threshold = float(tm.group(1))
    if abs(html_threshold - NETCTL_THRESHOLD) > 1e-9:
        raise ToolLoadError("parse_failed",
            f"{basename}: Threshold={html_threshold}Validation status{NETCTL_THRESHOLD}")

    soup = BeautifulSoup(raw, "html.parser")
    pre = soup.find("pre")
    text = pre.get_text() if pre else raw

    records = []
    for m in _NETCTL_LINE_RE.finditer(text):
        records.append({
            "residue_num": int(m.group(1)),
            "protein_id": m.group(2),
            "peptide": m.group(3).upper(),
            "comb": float(m.group(8)),
            "is_ligand": bool(m.group(9) and m.group(9).strip()),
        })
    if not records:
        raise ToolLoadError("parse_failed", f"{basename}Required input or value was not found")

    positions = [r["residue_num"] for r in records]
    if positions != list(range(1, max(positions) + 1)):
        raise ToolLoadError("parse_failed",
            f"{basename}Processing details"
            f"({len(positions)}Table status{max(positions)})")

    bad_len = [r["peptide"] for r in records if len(r["peptide"]) != 9]
    if bad_len:
        raise ToolLoadError("parse_failed", f"{basename}Peptide status{bad_len[:5]}")

    for r in records:
        expected_ligand = r["comb"] >= NETCTL_THRESHOLD
        if r["is_ligand"] != expected_ligand and abs(r["comb"] - NETCTL_THRESHOLD) >= 5e-4:
            raise ToolLoadError("parse_failed",
                f"{basename}: COMB={r['comb']}Processing details"
                f"({'Processing details' if r['is_ligand'] else 'Processing details'}Validation status"
                f"{r['peptide']}")

    protein_ids = set(r["protein_id"] for r in records)
    if len(protein_ids) != 1:
        raise ToolLoadError("parse_failed", f"{basename}Processing details{sorted(protein_ids)}")

    scores = {}
    for r in records:
        pep = r["peptide"]
        if pep not in scores or r["comb"] > scores[pep]:
            scores[pep] = r["comb"]

    signature = [(r["residue_num"], r["peptide"]) for r in records]
    return supertype, scores, next(iter(protein_ids)), signature


def load_tool_netctl_1_2(antigen: str) -> dict:
    directory = RAW_PATHS["NetCTL_1.2"][antigen]
    files = sorted(
        set(glob.glob(os.path.join(directory, "*.html")) +
            glob.glob(os.path.join(directory, "*.HTML")))
    )
    if not files:
        raise ToolLoadError("input_file_missing", f"Required input or value was not found{directory}")

    expected_count = NETCTL_EXPECTED_FILE_COUNTS.get(antigen)
    if expected_count is not None and len(files) != expected_count:
        raise ToolLoadError("parse_failed",
            f"NetCTL [{antigen}Processing details{len(files)}Validation status"
            f"{expected_count}Processing details")

    merged = {}
    supertype_to_file = {}
    reference_sig = reference_protein_id = reference_file = None

    for fp in files:
        supertype, scores, protein_id, signature = parse_netctl_html(fp)
        fname = os.path.basename(fp)

        if supertype in supertype_to_file:
            raise ToolLoadError("parse_failed",
                f"NetCTL [{antigen}Processing details{supertype}Processing details"
                f"{supertype_to_file[supertype]}Processing details{fname}")
        supertype_to_file[supertype] = fname

        if reference_sig is None:
            reference_sig, reference_protein_id, reference_file = signature, protein_id, fname
        else:
            if protein_id != reference_protein_id:
                raise ToolLoadError("parse_failed",
                    f"NetCTL [{antigen}]: {fname} — protein ID '{protein_id}Validation status"
                    f"Processing details{reference_protein_id}Processing details{reference_file}")
            if signature != reference_sig:
                raise ToolLoadError("parse_failed",
                    f"NetCTL [{antigen}]: {fname}Validation status"
                    f"Processing details{reference_file}Processing details")

        merged[supertype] = scores

    expected_supertypes = NETCTL_EXPECTED_SUPERTYPES.get(antigen)
    if expected_supertypes is not None:
        found = set(merged.keys())
        if found != expected_supertypes:
            missing, extra = expected_supertypes - found, found - expected_supertypes
            parts = []
            if missing:
                parts.append(f"Required input or value was not found{sorted(missing)}")
            if extra:
                parts.append(f"Processing details{sorted(extra)}")
            raise ToolLoadError("parse_failed", f"NetCTL [{antigen}Processing details" + "; ".join(parts))

    return merged


# ═══════════════════════════════════════════════════════════════════════════
# Implementation detail; see the repository documentation.
# ═══════════════════════════════════════════════════════════════════════════

TOOL_SPECS = [
    {
        "name": "IEDB_I_Consensus",
        "kind": "standard",
        "loader": lambda antigen: load_iedb_tsv_multi(
            RAW_PATHS["IEDB_I_Consensus"][antigen],
            ["median binding percentile", "consensus percentile", "smm ic50", "ann ic50"],
            (8, 14),
            {"smm ic50"},
            get_allele_sets("I", antigen)[0],
        ),
        "length_range": (8, 14),
        "scores": [
            {"key": "median binding percentile", "label": "Median_percentile", "higher_is_better": False},
            {"key": "consensus percentile",       "label": "Consensus_percentile", "higher_is_better": False},
            {"key": "smm ic50",                   "label": "SMM_IC50", "higher_is_better": False},
            {"key": "ann ic50",                   "label": "ANN_IC50", "higher_is_better": False},
        ],
    },
    {
        "name": "IEDB_I_NetMHCpan_4.1_BA",
        "kind": "standard",
        "loader": lambda antigen: load_iedb_tsv_multi(
            RAW_PATHS["IEDB_I_NetMHCpan_4.1_BA"][antigen],
            ["median binding percentile", "netmhcpan_ba ic50"],
            (8, 14),
            set(),
            get_allele_sets("I", antigen)[0],
        ),
        "length_range": (8, 14),
        "scores": [
            {"key": "median binding percentile", "label": "Median_percentile", "higher_is_better": False},
            {"key": "netmhcpan_ba ic50",          "label": "BA_IC50", "higher_is_better": False},
        ],
    },
    {
        "name": "IEDB_I_NetMHCpan_4.1_EL",
        "kind": "standard",
        "loader": lambda antigen: load_iedb_tsv_multi(
            RAW_PATHS["IEDB_I_NetMHCpan_4.1_EL"][antigen],
            ["median binding percentile"],
            (8, 14),
            set(),
            get_allele_sets("I", antigen)[0],
        ),
        "length_range": (8, 14),
        "scores": [
            {"key": "median binding percentile", "label": "Median_percentile", "higher_is_better": False},
        ],
    },
    {
        "name": "NetMHCpan_4.1",
        "kind": "standard",
        "loader": load_tool_netmhcpan_4_1,
        "length_range": (8, 14),
        "scores": [
            {"key": "rank_el", "label": "Rank_EL", "higher_is_better": False},
            {"key": "rank_ba", "label": "Rank_BA", "higher_is_better": False},
        ],
    },
    {
        "name": "NetMHC_4.0",
        "kind": "standard",
        "loader": load_tool_netmhc_4_0,
        "length_range": (8, 14),
        "scores": [
            {"key": "affinity_nM", "label": "Affinity_nM", "higher_is_better": False},
            {"key": "rank_pct",    "label": "Rank_pct", "higher_is_better": False},
        ],
    },
    {
        "name": "NetCTL_1.2",
        "kind": "netctl",
        "loader": load_tool_netctl_1_2,
        "length_range": (9, 9),
        "scores": [
            {"key": None, "label": "COMB", "higher_is_better": True},
        ],
    },
]

FATAL_STATUSES = {"input_file_missing", "parse_failed"}
REVIEW_STATUSES = {"epitope_score_missing", "incomplete_coverage",
                    "score_not_available_for_allele", "epitope_not_in_sequence"}
UNSUPPORTED_STATUSES = {"allele_not_supported", "unsupported_length", "no_supertype_mapping"}


# ═══════════════════════════════════════════════════════════════════════════
# Allele handling.
# ═══════════════════════════════════════════════════════════════════════════

def evaluate_pair(peptide: str, allele: str, sequence: str,
                   length_range, kind: str, higher_is_better: bool,
                   scalar_data: dict) -> dict:
    k = len(peptide)
    empty = {"status": None, "frank_pct": None, "n_better": None,
             "n_total": None, "epitope_score": None, "n_missing_windows": None}

    if length_range and not (length_range[0] <= k <= length_range[1]):
        empty["status"] = "unsupported_length"
        return empty

    if kind == "netctl":
        supertype = ALLELE_TO_SUPERTYPE.get(allele)
        if supertype is None:
            empty["status"] = "no_supertype_mapping"
            return empty
        lookup_key = supertype
    else:
        lookup_key = allele

    # Allele handling.
    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    if lookup_key not in scalar_data:
        empty["status"] = "allele_not_supported"
        return empty
    score_lookup = scalar_data[lookup_key]
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
    print(f"Processing details{antigen}  (MHC I)")
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
                if spec["kind"] == "netctl":
                    scalar_data = tool_data
                else:
                    scalar_data = extract_scalar_dict(tool_data, score_spec["key"])
                rows = []
                for allele, peptide in validated_pairs:
                    res = evaluate_pair(
                        peptide, allele, sequence,
                        spec["length_range"], spec["kind"],
                        score_spec["higher_is_better"], scalar_data,
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
        out_path = os.path.join(out_dir, f"FRANK_{antigen}_MHC_I.xlsx")
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
              f"Validation status{', '.join(antigen_failed)}")
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
