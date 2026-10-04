#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NetMHCIIpan 4.1 — exact allele-peptide RM/UP.

Thresholds: %Rank_EL <= 10, %Rank_BA <= 10, Affinity < 5000 nM, Union_OR.
Includes blocking parser-integrity checks and exact validation of the prediction
allele set against the configured run for p24, pp65, and PtxS1.
RM = passed prediction present in the positive reference dataset.
UP = passed prediction absent from the positive reference dataset.
No FN calculation or export is performed.
"""

from project_paths import DATA_ROOT

import os
import re
import unicodedata
from collections import Counter
from typing import Optional

import pandas as pd

BASE      = f"{DATA_ROOT}"
TOOL_NAME = "NetMHCIIpan_4.1"

ANTIGENS = {
    "p24": {
        "excel": f"{BASE}/p24/p24.xlsx",
        "xls_files": [
            f"{BASE}/p24/MHC II/NetMHCIIpan/NetMHCIIpan_p24_part_1_11-18.xls",
            f"{BASE}/p24/MHC II/NetMHCIIpan/NetMHCIIpan_p24_part_2_11-18.xls",
        ],
    },
    "pp65": {
        "excel": f"{BASE}/pp65/pp65.xlsx",
        "xls_files": [
            f"{BASE}/pp65/MHC II/NetMHCIIpan/NetMHCIIpan_pp65_11-18.xls",
        ],
    },
    "PtxS1": {
        "excel": f"{BASE}/PtxS1/PtxS1.xlsx",
        "xls_files": [
            f"{BASE}/PtxS1/MHC II/NetMHCIIpan/NetMHCIIpan_PtxS1_11-18.xls",
        ],
    },
}

THRESHOLDS = {
    "rank_el":  ("Rank_EL",     10,    "<="),
    "rank_ba":  ("Rank_BA",     10,    "<="),
    "affinity": ("Affinity_nM", 5000,  "<"),
}

MODES = ["rank_el", "rank_ba", "affinity", "union_or"]

MODE_SHEET_NAMES = {
    "rank_el":  "Rank_EL_lte10",
    "rank_ba":  "Rank_BA_lte10",
    "affinity": "Affinity_lt5000",
    "union_or": "Union_OR",
}

MODE_LABELS = {
    "rank_el":  "%Rank_EL <= 10",
    "rank_ba":  "%Rank_BA <= 10",
    "affinity": "Affinity < 5000",
    "union_or": "Union OR",
}

MODE_METRIC_COLS = {
    "rank_el":  ["Rank_EL"],
    "rank_ba":  ["Rank_BA"],
    "affinity": ["Affinity_nM"],
    "union_or": ["Rank_EL", "Rank_BA", "Affinity_nM"],
}

# Allele handling.
# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
EXPECTED_ALLELES = {
    'p24': [
        'DRB1_0101', 'DRB1_0301', 'DRB1_0302', 'DRB1_0401', 'DRB1_0404', 'DRB1_0405',
        'DRB1_0701', 'DRB1_0801', 'DRB1_0901', 'DRB1_1001', 'DRB1_1101', 'DRB1_1301',
        'DRB1_1302', 'DRB1_1303', 'DRB1_1304', 'DRB1_1401', 'DRB1_1501', 'DRB1_1502',
        'DRB3_0101', 'DRB3_0301', 'DRB3_0303', 'DRB4_0101', 'DRB5_0101',
    ],
    'pp65': [
        'DRB1_0301', 'DRB1_0401', 'DRB1_0402', 'DRB1_0404', 'DRB1_0701', 'DRB1_1101',
        'DRB1_1104', 'DRB1_1501', 'DRB3_0101', 'DRB3_0202',
    ],
    'PtxS1': [
        'DRB1_0101', 'DRB1_1101',
    ],
}


# ──────────────────────────────────────────────────────────────────────────────
# Normalization.
# ──────────────────────────────────────────────────────────────────────────────

def get_output_dir(antigen: str) -> str:
    return os.path.join(
        BASE, f"{antigen}",
        "Matches MHC II", f"Matches {TOOL_NAME}",
    )


def normalize_str(s: str) -> str:
    if not isinstance(s, str):
        s = str(s)
    return unicodedata.normalize("NFKC", s).strip()


def raw_allele_to_standard(raw: str) -> str:
    """Portable implementation of NetMHCIIpan_4_1_search_for_matches. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    raw = raw.strip()
    m = re.match(r'^((?:HLA-)?[A-Z]+\d*)[_*](\d{2})(\d{2})$', raw)
    if m:
        locus, g1, g2 = m.group(1), m.group(2), m.group(3)
        locus = locus.replace("HLA-", "")
        return f"{locus}*{g1}:{g2}"
    return raw


def normalize_allele(allele: str) -> str:
    allele = normalize_str(allele)
    allele = re.sub(r'[\u2012\u2013\u2014\u2015\u2212]', '-', allele)
    allele = re.sub(r'^HLA-', '', allele)
    return allele.strip()


# ──────────────────────────────────────────────────────────────────────────────
# Implementation detail; see the repository documentation.
# ──────────────────────────────────────────────────────────────────────────────

def is_binary_excel(fpath: str) -> bool:
    """Portable implementation of NetMHCIIpan_4_1_search_for_matches. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    with open(fpath, "rb") as fh:
        head = fh.read(8)
    return head.startswith(b"\xD0\xCF\x11\xE0") or head.startswith(b"PK\x03\x04")


def find_header_row_index(lines: list[str], max_scan: int = 10) -> Optional[int]:
    """Portable implementation of NetMHCIIpan_4_1_search_for_matches. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    for i, line in enumerate(lines[:max_scan]):
        first_cell = line.split("\t", 1)[0].strip()
        if first_cell == "Pos":
            return i
    return None


def _empty_parse_stats() -> dict:
    return {
        "corrupted_rows":   0,
        "truncated_blocks": 0,
        "rows_per_allele":  {},
        "nan_counts":       {"Rank_EL": 0, "Rank_BA": 0, "Affinity_nM": 0},
        "peptide_rows":     [],  # Peptide handling.
    }


def parse_xls_file(fpath: str) -> tuple[pd.DataFrame, dict]:
    """Portable implementation of NetMHCIIpan_4_1_search_for_matches. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    fname = os.path.basename(fpath)

    if is_binary_excel(fpath):
        print(f"    ОШИБКА: {fname} — это настоящий бинарный Excel-файл, а не "
              f"текстовый tab-separated .xls. Такой формат этот скрипт пока не "
              f"поддерживает — пришли пример, чтобы добавить его разбор.")
        return pd.DataFrame(), _empty_parse_stats()

    with open(fpath, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    header_idx = find_header_row_index(lines)
    if header_idx is None:
        print(f"    ПРЕДУПРЕЖДЕНИЕ: строка заголовков ('Pos') не найдена в {fname}")
        return pd.DataFrame(), _empty_parse_stats()
    if header_idx == 0:
        print(f"    ПРЕДУПРЕЖДЕНИЕ: над строкой заголовков нет строки с именами "
              f"аллелей в {fname}")
        return pd.DataFrame(), _empty_parse_stats()

    allele_row = lines[header_idx - 1].rstrip("\n").split("\t")
    header_row = [c.strip() for c in lines[header_idx].rstrip("\n").split("\t")]

    # Allele handling.
    block_starts = [i for i, v in enumerate(header_row) if v == "Core"]
    if not block_starts:
        print(f"    ПРЕДУПРЕЖДЕНИЕ: не найдено ни одного блока аллеля "
              f"(колонка 'Core') в {fname}")
        return pd.DataFrame(), _empty_parse_stats()

    try:
        ave_idx = header_row.index("Ave")
    except ValueError:
        ave_idx = len(header_row)
    block_ends = block_starts[1:] + [ave_idx]

    blocks = []
    for start, end in zip(block_starts, block_ends):
        allele_raw = allele_row[start].strip() if start < len(allele_row) else ""
        if not allele_raw:
            print(f"    ПРЕДУПРЕЖДЕНИЕ: пустое имя аллеля для блока в колонке "
                  f"{start} ({fname}) — блок пропущен")
            continue
        sub_headers = header_row[start:end]
        try:
            rank_el_idx = start + sub_headers.index("Rank")
            rank_ba_idx = start + sub_headers.index("Rank_BA")
            afin_idx    = start + sub_headers.index("nM")
        except ValueError:
            print(f"    ПРЕДУПРЕЖДЕНИЕ: неожиданная структура блока для аллеля "
                  f"{allele_raw} ({sub_headers}) в {fname} — блок пропущен")
            continue
        blocks.append({
            "allele_raw":  allele_raw,
            "allele_std":  raw_allele_to_standard(allele_raw),
            "rank_el_idx": rank_el_idx,
            "rank_ba_idx": rank_ba_idx,
            "afin_idx":    afin_idx,
        })

    if not blocks:
        print(f"    ПРЕДУПРЕЖДЕНИЕ: ни один блок аллеля не удалось разобрать в {fname}")
        return pd.DataFrame(), _empty_parse_stats()

    def to_float(s: str) -> float:
        s = s.strip()
        if s in ("", "NA"):
            return float("nan")
        try:
            return float(s)
        except ValueError:
            return float("nan")

    stats = _empty_parse_stats()
    stats["rows_per_allele"] = {blk["allele_raw"]: 0 for blk in blocks}

    rows = []
    for line in lines[header_idx + 1:]:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        parts = line.split("\t")
        try:
            int(parts[0])
        except (ValueError, IndexError):
            stats["corrupted_rows"] += 1
            continue
        if len(parts) < 2:
            stats["corrupted_rows"] += 1
            continue
        peptide = parts[1].strip()
        stats["peptide_rows"].append(peptide)

        for blk in blocks:
            needed_max = max(blk["rank_el_idx"], blk["rank_ba_idx"], blk["afin_idx"])
            if needed_max >= len(parts):
                stats["truncated_blocks"] += 1
                continue
            rank_el = to_float(parts[blk["rank_el_idx"]])
            rank_ba = to_float(parts[blk["rank_ba_idx"]])
            afin    = to_float(parts[blk["afin_idx"]])
            if pd.isna(rank_el):
                stats["nan_counts"]["Rank_EL"] += 1
            if pd.isna(rank_ba):
                stats["nan_counts"]["Rank_BA"] += 1
            if pd.isna(afin):
                stats["nan_counts"]["Affinity_nM"] += 1
            stats["rows_per_allele"][blk["allele_raw"]] += 1
            rows.append({
                "Allele_raw":  blk["allele_raw"],
                "Allele":      blk["allele_std"],
                "Peptide":     peptide,
                "Rank_EL":     rank_el,
                "Rank_BA":     rank_ba,
                "Affinity_nM": afin,
                "source_file": fname,
            })

    return pd.DataFrame(rows), stats


def parse_xls_files(xls_paths: list[str]) -> tuple[pd.DataFrame, dict]:
    """Portable implementation of NetMHCIIpan_4_1_search_for_matches. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    dfs = []
    file_stats = {}
    for fpath in xls_paths:
        print(f"  Парсинг xls: {os.path.basename(fpath)}")
        df, stats = parse_xls_file(fpath)
        fname = os.path.basename(fpath)
        file_stats[fname] = stats
        if df.empty:
            print(f"    ПРЕДУПРЕЖДЕНИЕ: строки не извлечены из {fname}")
        else:
            print(f"    строк: {len(df):,}   аллелей в файле: {df['Allele_raw'].nunique()}")
        dfs.append(df)

    non_empty = [d for d in dfs if not d.empty]
    if not non_empty:
        print("  ПРЕДУПРЕЖДЕНИЕ: предсказания не найдены ни в одном xls-файле!")
        return pd.DataFrame(), file_stats

    df_all = pd.concat(non_empty, ignore_index=True)
    print(f"  Загружено строк предсказания (всего по антигену): {len(df_all):,}")
    return df_all, file_stats


def check_parse_integrity(antigen_name: str, file_stats: dict) -> bool:
    """Portable implementation of NetMHCIIpan_4_1_search_for_matches. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    print(f"\n  --- Проверка целостности разбора xls для {antigen_name} ---")
    ok = True
    peptide_lists = {}

    for fname, stats in file_stats.items():
        print(f"  [{fname}]")
        if stats["corrupted_rows"] > 0:
            ok = False
            print(f"    ✗ повреждённых строк (Pos/Peptide не распознаны): {stats['corrupted_rows']}")
        if stats["truncated_blocks"] > 0:
            ok = False
            print(f"    ✗ обрезанных блоков аллель×строка (не хватает колонок): {stats['truncated_blocks']}")
        nan_total = sum(stats["nan_counts"].values())
        if nan_total > 0:
            ok = False
            print(f"    ✗ NaN в метриках: {stats['nan_counts']}")

        rpa = stats["rows_per_allele"]
        counts = set(rpa.values())
        if len(counts) > 1:
            ok = False
            print(f"    ✗ число строк по аллелям внутри файла НЕ одинаково: {rpa}")
        elif rpa:
            n = next(iter(counts))
            print(f"    ✓ строк на аллель: {n} (одинаково для всех {len(rpa)} аллелей файла)")

        print(f"    строк данных прочитано: {len(stats['peptide_rows'])}")
        peptide_lists[fname] = stats["peptide_rows"]

    if len(peptide_lists) > 1:
        names = list(peptide_lists.keys())
        base_name, base_list = names[0], peptide_lists[names[0]]
        mismatch_found = False
        for other_name in names[1:]:
            other_list = peptide_lists[other_name]
            if other_list != base_list:
                ok = False
                mismatch_found = True
                print(f"    ✗ последовательности пептидных строк в {base_name} и "
                      f"{other_name} НЕ совпадают: {len(base_list)} против "
                      f"{len(other_list)} строк")

                first_diff = next(
                    (i for i, (a, b) in enumerate(zip(base_list, other_list)) if a != b),
                    None,
                )
                if first_diff is not None:
                    print(f"        первое расхождение на позиции {first_diff}: "
                          f"'{base_list[first_diff]}' ({base_name}) vs "
                          f"'{other_list[first_diff]}' ({other_name})")

                c_base, c_other = Counter(base_list), Counter(other_list)
                dup_diffs = {
                    p: (c_base[p], c_other[p])
                    for p in (set(c_base) | set(c_other))
                    if c_base[p] != c_other[p]
                }
                if dup_diffs:
                    sample = sorted(dup_diffs.items())[:5]
                    tail = f" (+{len(dup_diffs) - 5} ещё)" if len(dup_diffs) > 5 else ""
                    print(f"        пептиды с разным числом повторов "
                          f"(пептид: [{base_name}], [{other_name}]): {sample}{tail}")
        if not mismatch_found:
            print(f"    ✓ последовательности пептидных строк совпадают между "
                  f"{len(names)} файлами ({len(base_list)} строк, включая дубли)")

    print(f"  {'-'*56}")
    return ok


# ──────────────────────────────────────────────────────────────────────────────
# Validation.
# ──────────────────────────────────────────────────────────────────────────────



def validate_alleles(antigen_name: str, pred_all: pd.DataFrame) -> bool:
    expected=set(EXPECTED_ALLELES[antigen_name])
    by_file=pred_all.groupby("Allele_raw")["source_file"].unique()
    found=set(by_file.index)
    missing,extra=sorted(expected-found),sorted(found-expected)
    duplicated={a:sorted(files) for a,files in by_file.items() if a in expected and len(files)>1}
    print(f"\n  Prediction allele check [{antigen_name}]: expected={len(expected)}, found={len(found)}")
    if missing: print(f"    отсутствуют: {', '.join(missing)}")
    if extra: print(f"    лишние: {', '.join(extra)}")
    if duplicated:
        for a,files in duplicated.items(): print(f"    duplicate {a}: {', '.join(files)}")
    return not missing and not extra and not duplicated


# ──────────────────────────────────────────────────────────────────────────────
# Implementation detail; see the repository documentation.
# ──────────────────────────────────────────────────────────────────────────────



def read_validated_epitopes(excel_path: str) -> pd.DataFrame:
    df_raw=pd.read_excel(excel_path,sheet_name=0,header=0)
    allele_col,peptide_col="HLA allele.1","HTL epitopes"
    missing=[c for c in (allele_col,peptide_col) if c not in df_raw.columns]
    if missing: raise ValueError(f"В {excel_path} отсутствуют колонки: {missing}")
    rows=[]
    for _,row in df_raw[[allele_col,peptide_col]].dropna().iterrows():
        allele=normalize_allele(str(row[allele_col])); peptide=normalize_str(str(row[peptide_col])).upper()
        if re.match(r'^(DRB|DQB|DPB|HLA-D)',allele) and len(peptide)>=8 and re.match(r'^[A-Z]+$',peptide): rows.append({"Allele":allele,"Peptide":peptide})
    return pd.DataFrame(rows).drop_duplicates()


# ──────────────────────────────────────────────────────────────────────────────
# Threshold handling.
# ──────────────────────────────────────────────────────────────────────────────

def _single_mask(df: pd.DataFrame, thr_key: str) -> pd.Series:
    col, val, op = THRESHOLDS[thr_key]
    if op == "<=":
        return df[col] <= val
    return df[col] < val


def _or_mask(df: pd.DataFrame) -> pd.Series:
    mask = pd.Series(False, index=df.index)
    for thr_key in THRESHOLDS:
        mask = mask | _single_mask(df, thr_key)
    return mask


def get_mode_mask(mode: str, df: pd.DataFrame) -> pd.Series:
    if mode == "union_or":
        return _or_mask(df)
    return _single_mask(df, mode)


def _union_flags_vectorized(df: pd.DataFrame) -> pd.DataFrame:
    """Portable implementation of NetMHCIIpan_4_1_search_for_matches. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    out = pd.DataFrame(index=df.index)
    flag_cols = []
    for thr_key, colname in (("rank_el", "Pass_Rank_EL"),
                              ("rank_ba", "Pass_Rank_BA"),
                              ("affinity", "Pass_Affinity")):
        col, val, op = THRESHOLDS[thr_key]
        passed = (df[col] <= val) if op == "<=" else (df[col] < val)
        passed = passed.astype("boolean")
        passed[df[col].isna()] = pd.NA
        out[colname] = passed
        flag_cols.append(colname)
    out["Passed_Thresholds_Count"] = out[flag_cols].sum(axis=1).astype("Int64")
    return out




def build_mode_frames(mode: str, pred_all: pd.DataFrame, val_set: set) -> dict:
    passed=pred_all[get_mode_mask(mode,pred_all)].copy()
    agg={"Rank_EL":"min","Rank_BA":"min","Affinity_nM":"min","source_file":lambda s:"; ".join(dict.fromkeys(map(str,s)))}
    pair_df=passed.groupby(["Allele","Peptide"],as_index=False).agg(agg)
    pairs=set(zip(pair_df["Allele"],pair_df["Peptide"]))
    rm_pairs,up_pairs=val_set & pairs,pairs-val_set
    is_union=mode=="union_or"
    if is_union and not pair_df.empty:
        flags=_union_flags_vectorized(pair_df)
        for c in flags.columns: pair_df[c]=flags[c]
    pair_df["Source File"]=pair_df["source_file"] if "source_file" in pair_df else ""
    cols=["Allele","Peptide"]+MODE_METRIC_COLS[mode]
    if is_union: cols += ["Pass_Rank_EL","Pass_Rank_BA","Pass_Affinity","Passed_Thresholds_Count"]
    cols += ["Source File"]
    idx=pd.MultiIndex.from_frame(pair_df[["Allele","Peptide"]]) if not pair_df.empty else None
    def subset(ps):
        if not ps: return pd.DataFrame(columns=cols)
        wanted=pd.MultiIndex.from_tuples(sorted(ps),names=["Allele","Peptide"])
        return pair_df[idx.isin(wanted)][cols].reset_index(drop=True)
    return {"rm_df":subset(rm_pairs),"up_df":subset(up_pairs),"rm_pairs":rm_pairs,"up_pairs":up_pairs}


# ──────────────────────────────────────────────────────────────────────────────
# Output generation.
# ──────────────────────────────────────────────────────────────────────────────

def save_excel_multi(sheets: dict[str, pd.DataFrame], path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    row_counts = ", ".join(f"{name}={len(df)}" for name, df in sheets.items())
    print(f"  Сохранено: {path}\n    ({row_counts})")


# ──────────────────────────────────────────────────────────────────────────────
# Implementation detail; see the repository documentation.
# ──────────────────────────────────────────────────────────────────────────────



def process_antigen(antigen_name: str, config: dict) -> str:
    print(f"\n{'='*60}\nОбработка антигена: {antigen_name}\n{'='*60}")
    xls_files=config.get("xls_files") or []
    if not xls_files: return "остановлен — xls_files не заполнен"
    missing=[f for f in xls_files if not os.path.isfile(f)]
    if missing: return "остановлен — xls-файлы не найдены"
    if not os.path.isfile(config["excel"]): return "остановлен — reference Excel не найден"
    pred_all,file_stats=parse_xls_files(xls_files)
    validated=read_validated_epitopes(config["excel"])
    if pred_all.empty or validated.empty: return "остановлен — пустые входные данные"
    if not check_parse_integrity(antigen_name,file_stats): return "остановлен — целостность разбора xls нарушена"
    if not validate_alleles(antigen_name,pred_all): return "остановлен — несоответствие prediction alleles"
    val_set=set(zip(validated["Allele"],validated["Peptide"]))
    results={m:build_mode_frames(m,pred_all,val_set) for m in MODES}
    print(f"\n  {'Режим':<18}{'RM':>8}{'UP':>10}")
    for m in MODES: print(f"  {MODE_LABELS[m]:<18}{len(results[m]['rm_pairs']):>8}{len(results[m]['up_pairs']):>10}")
    out_dir=get_output_dir(antigen_name)
    for kind,key in (("RM","rm_df"),("UP","up_df")):
        save_excel_multi({MODE_SHEET_NAMES[m]:results[m][key] for m in MODES},os.path.join(out_dir,f"{kind}_{antigen_name}_{TOOL_NAME}.xlsx"))
    return "успешно"


def main():
    print("NetMHCIIpan 4.1 — поиск совпадений с валидированными HTL эпитопами")
    print("Пороги считаются отдельно (%Rank_EL <= 10 / %Rank_BA <= 10 / "
          "Affinity < 5000 nM) и как объединение (Union OR)\n")

    statuses = {}
    for antigen_name, config in ANTIGENS.items():
        statuses[antigen_name] = process_antigen(antigen_name, config)

    print(f"\n{'='*60}")
    print("Итог по антигенам:")
    any_error = False
    for antigen_name, status in statuses.items():
        ok = (status == "успешно")
        any_error = any_error or not ok
        print(f"  {'✓' if ok else '✗'} {antigen_name}: {status}")

    print()
    if any_error:
        print("Завершено с ошибками. Проверь сообщения выше.")
    else:
        print("Готово! Все антигены обработаны успешно.")


if __name__ == "__main__":
    main()
