#!/usr/bin/env python3
"""
NetCTL 1.2: exact allele-peptide RM/UP from supertype-level predictions.

RM = passed COMB >= 0.75 prediction whose expanded exact allele-peptide pair
is present in the positive reference dataset.
UP = passed prediction whose expanded exact pair is absent from the reference.
No FN calculation or export is performed.

NetCTL predictions are generated per supertype; each passing supertype-peptide
is expanded only to exact HLA-A/B alleles represented in the reference dataset
and mapped to that supertype.
"""

from project_paths import DATA_ROOT

import os
import re
import glob
import sys
import pandas as pd
from bs4 import BeautifulSoup
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

# =============================================================================
# Configuration.
# =============================================================================

THRESHOLD = 0.75
BASE = f"{DATA_ROOT}"
TOOL = "NetCTL_1.2"

# Implementation detail; see the repository documentation.
EXPECTED_FILE_COUNTS = {
    "p24":   12,
    "pp65":  9,
    "PtxS1": 4,
}

ANTIGEN_CONFIG = {
    'p24': {
        "netctl_dir": os.path.join(BASE, "p24/MHC I/NetCTL"),
        "excel": os.path.join(BASE, "p24/p24.xlsx"),
        "out": os.path.join(BASE, "p24/Matches MHC I/Matches NetCTL_1.2"),
        "expected_supertypes": {'A1', 'A2', 'A3', 'A24', 'A26', 'B7', 'B8', 'B27', 'B39', 'B44', 'B58', 'B62'},
    },
    'pp65': {
        "netctl_dir": os.path.join(BASE, "pp65/MHC I/NetCTL"),
        "excel": os.path.join(BASE, "pp65/pp65.xlsx"),
        "out": os.path.join(BASE, "pp65/Matches MHC I/Matches NetCTL_1.2"),
        "expected_supertypes": {'A1', 'A2', 'A3', 'A24', 'A26', 'B7', 'B44', 'B58', 'B62'},
    },
    'PtxS1': {
        "netctl_dir": os.path.join(BASE, "PtxS1/MHC I/NetCTL"),
        "excel": os.path.join(BASE, "PtxS1/PtxS1.xlsx"),
        "out": os.path.join(BASE, "PtxS1/Matches MHC I/Matches NetCTL_1.2"),
        "expected_supertypes": {'A2', 'B7', 'B44', 'B58'},
    },
}

# =============================================================================
# Allele handling.
# =============================================================================

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

# Implementation detail; see the repository documentation.
SUPERTYPE_TO_ALLELES: dict = {}
for _a, _st in ALLELE_TO_SUPERTYPE.items():
    SUPERTYPE_TO_ALLELES.setdefault(_st, []).append(_a)

# Normalization.
_CE_PATTERN = re.compile(r'^[CE]\*')

# =============================================================================
# Implementation detail; see the repository documentation.
# =============================================================================

def extract_supertype_from_html(raw: str, filepath: str) -> tuple:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    basename = os.path.basename(filepath)
    m_html = re.search(r'using\s+MHC\s+supertype\s+(\w+)', raw, re.IGNORECASE)

    m_fname = re.search(
        r'NetCTL[_\-]1[._]2[_\-]([A-Za-z0-9]+)(?:\(\d+\))?\.html?$',
        basename, re.IGNORECASE
    )

    if m_html:
        st_html  = m_html.group(1).upper()
        if not m_fname:
            print(f"Error"
                  f"{basename}")
            return st_html, False
        st_fname = m_fname.group(1)
        if st_html != st_fname.upper():
            print(f"Error"
                  f"Processing details{st_html}Processing details{st_fname}' "
                  f"Processing details{basename})")
            return st_html, False
        return st_html, True

    # Implementation detail; see the repository documentation.
    if m_fname:
        print(f"Required input or value was not found{basename}'; "
              f"Processing details{m_fname.group(1)}")
        return m_fname.group(1), True

    print(f"Error{basename}")
    return "unknown", False


def validate_html_records(records: list, filepath: str) -> None:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if not records:
        return

    basename = os.path.basename(filepath)

    # 1. residue numbers
    positions = [r["residue_num"] for r in records]
    expected_pos = list(range(1, max(positions) + 1))
    if positions != expected_pos:
        raise ValueError(
            f"{basename}Processing details"
            f"Processing details{len(positions)}Validation status{len(expected_pos)})"
        )

    # Peptide handling.
    bad_len = [r["peptide"] for r in records if len(r["peptide"]) != 9]
    if bad_len:
        raise ValueError(
            f"{basename}Peptide status"
            f"{bad_len[:5]}{'...' if len(bad_len) > 5 else ''}"
        )

    # Threshold handling.
    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    # Output generation.
    for r in records:
        expected_ligand = r["comb"] >= THRESHOLD
        if r["is_ligand"] != expected_ligand:
            if abs(r["comb"] - THRESHOLD) < 5e-4:
                print(f"Warning{basename}Processing details{r['comb']:.4f} — "
                      f"Processing details{'Processing details' if r['is_ligand'] else 'Processing details'} "
                      f"Peptide status{r['peptide']}")
            else:
                raise ValueError(
                    f"{basename}: COMB={r['comb']}Processing details"
                    f"({'Processing details' if r['is_ligand'] else 'Processing details'}Validation status"
                    f"Peptide status{r['peptide']}"
                )

    # Implementation detail; see the repository documentation.
    # Validation.


def parse_netctl_html(filepath: str) -> tuple:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
        raw = fh.read()

    supertype, fname_ok = extract_supertype_from_html(raw, filepath)

    # Validation.
    tm = re.search(r'Threshold\s+([0-9.]+)', raw, re.IGNORECASE)
    if not tm:
        raise ValueError(f"{os.path.basename(filepath)}Required input or value was not found")
    html_threshold = float(tm.group(1))
    if abs(html_threshold - THRESHOLD) > 1e-9:
        raise ValueError(
            f"{os.path.basename(filepath)}Processing details{html_threshold}, "
            f"Validation status{THRESHOLD}"
        )

    soup = BeautifulSoup(raw, "html.parser")
    pre  = soup.find("pre")
    text = pre.get_text() if pre else raw

    pat = re.compile(
        r'^\s*(\d+)\s+ID\s+(\S+)\s+pep\s+(\S+)'
        r'\s+aff\s+([\d.\-]+)\s+aff_rescale\s+([\d.\-]+)'
        r'\s+cle\s+([\d.\-]+)\s+tap\s+([\d.\-]+)'
        r'\s+COMB\s+([\d.\-]+)'
        r'(\s+<-E)?',
        re.MULTILINE
    )
    records = []
    for m in pat.finditer(text):
        records.append({
            "residue_num": int(m.group(1)),
            "protein_id":  m.group(2),
            "peptide":     m.group(3).upper(),
            "aff":         float(m.group(4)),
            "aff_rescale": float(m.group(5)),
            "cle":         float(m.group(6)),
            "tap":         float(m.group(7)),
            "comb":        float(m.group(8)),
            "is_ligand":   bool(m.group(9) and m.group(9).strip()),
            "supertype":   supertype,
            "source_file": os.path.basename(filepath),
        })

    validate_html_records(records, filepath)

    return pd.DataFrame(records), fname_ok


def load_all_netctl(netctl_dir: str, antigen: str,
                    expected_supertypes: set) -> pd.DataFrame:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    # Implementation detail; see the repository documentation.
    files = sorted(set(
        glob.glob(os.path.join(netctl_dir, "*.html")) +
        glob.glob(os.path.join(netctl_dir, "*.HTML"))
    ))
    if not files:
        raise FileNotFoundError(f"Processing details{netctl_dir}")

    expected_count = EXPECTED_FILE_COUNTS.get(antigen, len(expected_supertypes))
    print(f"Processing details{len(files)}Validation status{expected_count})")
    if len(files) != expected_count:
        raise ValueError(
            f"Processing details{len(files)}, "
            f"Validation status{expected_count}"
        )

    frames = []
    supertype_to_files: dict = {}
    all_fname_ok = True
    reference_signature  = None   # Implementation detail; see the repository documentation.
    reference_protein_id = None
    reference_file       = None

    for fp in files:
        df, fname_ok = parse_netctl_html(fp)
        if not fname_ok:
            all_fname_ok = False
        if not df.empty:
            # Validation.
            protein_ids = set(df["protein_id"].astype(str))
            if len(protein_ids) != 1:
                raise ValueError(
                    f"{os.path.basename(fp)}Processing details"
                    f"{', '.join(sorted(protein_ids))}"
                )
            protein_id = next(iter(protein_ids))

            # Peptide handling.
            signature = list(zip(
                df["residue_num"].astype(int),
                df["peptide"].astype(str),
            ))

            if reference_signature is None:
                reference_signature  = signature
                reference_protein_id = protein_id
                reference_file       = os.path.basename(fp)
            else:
                if protein_id != reference_protein_id:
                    raise ValueError(
                        f"{os.path.basename(fp)}: protein ID '{protein_id}' "
                        f"Validation status{reference_protein_id}' "
                        f"Processing details{reference_file}. "
                        f"Processing details"
                    )
                if signature != reference_signature:
                    raise ValueError(
                        f"{os.path.basename(fp)}Processing details"
                        f"Validation status{reference_file}. "
                        f"Processing details"
                        f"Processing details"
                    )

            frames.append(df)
            st = df["supertype"].iloc[0]
            supertype_to_files.setdefault(st, []).append(os.path.basename(fp))
            print(f"Processing details{os.path.basename(fp)}: {len(df)}Peptide status"
                  f"Processing details{st}]")
        else:
            raise ValueError(f"Table status"
                             f"{os.path.basename(fp)}")

    if not frames:
        raise ValueError(f"Processing details{netctl_dir}Processing details")

    # Deduplication.
    duplicates = {st: fn for st, fn in supertype_to_files.items() if len(fn) > 1}
    if duplicates:
        msgs = [f"{st}: {', '.join(fn)}" for st, fn in duplicates.items()]
        raise ValueError(f"Processing details{'; '.join(msgs)}")

    pred = pd.concat(frames, ignore_index=True)
    found_supertypes = set(pred["supertype"].dropna().unique())

    extra   = found_supertypes - expected_supertypes
    missing = expected_supertypes - found_supertypes
    if extra or missing:
        parts = []
        if extra:
            parts.append(f"Processing details{', '.join(sorted(extra))}]")
        if missing:
            parts.append(f"Required input or value was not found{', '.join(sorted(missing))}]")
        raise ValueError(f"Processing details{'; '.join(parts)}")

    if not all_fname_ok:
        raise ValueError(
            "Processing details"
            "Error"
        )

    print(f"Processing details{len(files)}Validation status"
          f"{', '.join(sorted(expected_supertypes))}")
    return pred


# =============================================================================
# Implementation detail; see the repository documentation.
# =============================================================================



def load_experimental_data(excel_path: str) -> pd.DataFrame:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    import openpyxl as ox
    wb = ox.load_workbook(excel_path, read_only=True, data_only=True)
    ws = wb.active
    rows_iter = ws.iter_rows(values_only=True)
    header = next(rows_iter, None)
    if header is None:
        wb.close(); raise ValueError(f"Processing details{excel_path}")
    headers = [str(x).strip() if x is not None else "" for x in header]
    try:
        allele_i = headers.index("HLA allele")
        peptide_i = headers.index("CTL epitopes")
    except ValueError:
        wb.close(); raise ValueError(f"Processing details{excel_path}Table status{headers}")
    rows=[]
    for row in rows_iter:
        allele = row[allele_i] if allele_i < len(row) else None
        peptide = row[peptide_i] if peptide_i < len(row) else None
        if not allele or not peptide: continue
        allele=str(allele).strip().replace("HLA-","")
        peptide=str(peptide).strip().upper()
        if allele and peptide: rows.append({"allele_norm":allele,"peptide":peptide})
    wb.close()
    df=pd.DataFrame(rows).drop_duplicates(subset=["allele_norm","peptide"])
    if df.empty: raise ValueError(f"Processing details{excel_path}")
    df["supertype"]=df["allele_norm"].map(ALLELE_TO_SUPERTYPE)
    return df.reset_index(drop=True)


# =============================================================================
# Output generation.
# =============================================================================

HEADER_FILL = PatternFill("solid", start_color="1F4E79", end_color="1F4E79")
HEADER_FONT = Font(name="Arial", bold=True, color="FFFFFF", size=11)
RM_FILL     = PatternFill("solid", start_color="C6EFCE", end_color="C6EFCE")
UP_FILL     = PatternFill("solid", start_color="FFEB9C", end_color="FFEB9C")
ROW_FONT    = Font(name="Arial", size=10)
CENTER      = Alignment(horizontal="center", vertical="center")


def write_excel(df: pd.DataFrame, filepath: str,
                sheet_name: str, row_fill: PatternFill) -> None:
    wb = Workbook()
    ws = wb.active
    ws.title = sheet_name

    cols = list(df.columns)

    # Implementation detail; see the repository documentation.
    for ci, col in enumerate(cols, start=1):
        cell = ws.cell(row=1, column=ci, value=col)
        cell.fill      = HEADER_FILL
        cell.font      = HEADER_FONT
        cell.alignment = CENTER

    if df.empty:
        for ci, col in enumerate(cols, start=1):
            ws.column_dimensions[get_column_letter(ci)].width = min(len(str(col)) + 4, 50)
        ws.freeze_panes = "A2"
        wb.save(filepath)
        print(f"Saved output{os.path.basename(filepath)}Table status")
        return

    for ri, (_, row) in enumerate(df.iterrows(), start=2):
        for ci, col in enumerate(cols, start=1):
            cell = ws.cell(row=ri, column=ci, value=row[col])
            cell.fill      = row_fill
            cell.font      = ROW_FONT
            cell.alignment = CENTER

    for ci, col in enumerate(cols, start=1):
        max_len = max(
            len(str(col)),
            *(len(str(df[col].iloc[i])) for i in range(len(df)))
        )
        ws.column_dimensions[get_column_letter(ci)].width = min(max_len + 4, 50)

    ws.freeze_panes = "A2"
    wb.save(filepath)
    print(f"Saved output{os.path.basename(filepath)}  ({len(df)}Table status")


# =============================================================================
# Main entry point.
# =============================================================================



def process_antigen(antigen: str, cfg: dict) -> None:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    print(f"\n{'='*60}")
    print(f"Processing details{antigen}")
    print(f"{'='*60}")

    print(f"Loaded input{cfg['excel']}")
    exp = load_experimental_data(cfg["excel"])
    exp_exact = set(zip(exp["allele_norm"], exp["peptide"]))
    run_supertypes = set(cfg["expected_supertypes"])
    evaluable = exp[exp["supertype"].isin(run_supertypes)].copy()
    evaluable_alleles = set(evaluable["allele_norm"].unique())
    print(f"Processing details{len(exp_exact)}")
    print(f"Processing details{len(evaluable_alleles)}")

    print(f"Loaded input{cfg['netctl_dir']}")
    pred = load_all_netctl(cfg["netctl_dir"], antigen, run_supertypes)
    pred["peptide_up"] = pred["peptide"].astype(str).str.upper()
    pred_pass = pred[pred["comb"] >= THRESHOLD].copy()
    print(f"Processing details{len(pred)}")
    print(f"Processing details{THRESHOLD}: {len(pred_pass)}")

    best_idx = pred_pass.groupby(["supertype", "peptide_up"])["comb"].idxmax()
    pred_best = pred_pass.loc[best_idx].copy()

    rm_rows, up_rows = [], []
    seen_rm, seen_up = set(), set()
    for _, r in pred_best.iterrows():
        st = r["supertype"]
        peptide = r["peptide_up"]
        alleles = sorted(a for a in evaluable_alleles if ALLELE_TO_SUPERTYPE.get(a) == st)
        for allele in alleles:
            key = (allele, peptide)
            row = {
                "Allele": allele,
                "Supertype": st,
                "Peptide": peptide,
                "COMB_score": round(float(r["comb"]), 4),
                "Source File": r["source_file"],
            }
            if key in exp_exact:
                if key not in seen_rm:
                    rm_rows.append(row); seen_rm.add(key)
            else:
                if key not in seen_up:
                    up_rows.append(row); seen_up.add(key)

    columns = ["Allele", "Supertype", "Peptide", "COMB_score", "Source File"]
    rm_df = pd.DataFrame(rm_rows, columns=columns)
    up_df = pd.DataFrame(up_rows, columns=columns)
    print(f"Processing details{len(rm_df)}  UP={len(up_df)}")

    os.makedirs(cfg["out"], exist_ok=True)
    write_excel(rm_df, os.path.join(cfg["out"], f"RM_{antigen}_{TOOL}.xlsx"), "Recovered Matches", RM_FILL)
    write_excel(up_df, os.path.join(cfg["out"], f"UP_{antigen}_{TOOL}.xlsx"), "Unvalidated Predictions", UP_FILL)
    print(f"Saved output{cfg['out']}")


def validate_config() -> None:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    for antigen, cfg in ANTIGEN_CONFIG.items():
        expected = set(cfg["expected_supertypes"])
        expected_count = EXPECTED_FILE_COUNTS[antigen]
        if len(expected) != expected_count:
            raise ValueError(
                f"{antigen}: {len(expected)}Processing details{expected_count}"
            )
        unknown = expected - set(SUPERTYPE_TO_ALLELES)
        if unknown:
            raise ValueError(f"{antigen}Processing details{sorted(unknown)}")
    print("Validation status")


def main():
    print(f"Validation status")
    print(f"Processing details{THRESHOLD}")

    try:
        validate_config()
    except ValueError as e:
        print(f"Error{e}")
        sys.exit(1)

    failed_antigens = []

    for antigen, cfg in ANTIGEN_CONFIG.items():
        try:
            process_antigen(antigen, cfg)
        except (FileNotFoundError, ValueError) as e:
            print(f"Error{antigen}: {e}")
            print(f"Processing details{antigen}Processing details"
                  f"Processing details")
            failed_antigens.append(antigen)
        except Exception as e:
            import traceback
            print(f"Error{antigen}: {e}")
            traceback.print_exc()
            failed_antigens.append(antigen)

    print(f"\n{'='*60}")
    if failed_antigens:
        print(f"Error"
              f"{', '.join(failed_antigens)}")
        sys.exit(1)
    else:
        print(f"Processing status")
    print("Completed successfully")


if __name__ == "__main__":
    main()
