#!/usr/bin/env python3
"""English documentation for this module or helper is provided in the repository README and in the surrounding code."""

# Implementation detail; see the repository documentation.

EVALUABLE_SET_MHC_I = {
    'p24': {
        'A*01:01', 'A*02:01', 'A*11:01', 'A*24:02', 'A*26:01', 'A*30:02', 'A*68:01', 'A*68:02',
        'B*07:02', 'B*08:01', 'B*15:01', 'B*35:01', 'B*40:01', 'B*44:02', 'B*44:03', 'B*51:01',
        'B*53:01', 'B*57:01', 'B*58:01'
    },
    'pp65': {
        'A*01:01', 'A*02:01', 'A*11:01', 'A*24:02', 'A*26:01', 'A*30:01', 'A*32:01', 'A*68:01',
        'B*07:02', 'B*15:01', 'B*35:01', 'B*40:01', 'B*44:02', 'B*44:03', 'B*51:01', 'B*53:01',
        'B*57:01', 'B*58:01'
    },
    'PtxS1': {
        'A*02:01', 'B*07:02', 'B*35:01', 'B*44:03', 'B*53:01', 'B*57:01'
    },
}

FULL_REFERENCE_SET_MHC_I = {
    'p24': {
        'A*01:01', 'A*02:01', 'A*02:07', 'A*11:01', 'A*11:03', 'A*24:02', 'A*24:07', 'A*25:01',
        'A*26:01', 'A*26:02', 'A*26:03', 'A*26:38', 'A*30:01', 'A*30:02', 'A*32:01', 'A*33:03',
        'A*68:01', 'A*68:02', 'A*74:01', 'B*07:02', 'B*08:01', 'B*13:02', 'B*14:01', 'B*14:02',
        'B*14:03', 'B*15:01', 'B*15:02', 'B*15:03', 'B*15:10', 'B*15:16', 'B*15:17', 'B*15:24',
        'B*15:40', 'B*27:05', 'B*35:01', 'B*35:02', 'B*35:03', 'B*35:05', 'B*35:08', 'B*35:11',
        'B*39:01', 'B*39:10', 'B*40:01', 'B*40:02', 'B*40:06', 'B*42:01', 'B*42:02', 'B*44:02',
        'B*44:03', 'B*44:15', 'B*45:01', 'B*48:01', 'B*50:01', 'B*51:01', 'B*52:01', 'B*53:01',
        'B*57:01', 'B*57:02', 'B*57:03', 'B*58:01', 'B*58:02', 'B*67:01', 'C*01:02', 'C*03:03',
        'C*03:04', 'C*04:01', 'C*06:02', 'C*07:01', 'C*08:02', 'C*18:01', 'E*01:01', 'E*01:03'
    },
    'pp65': {
        'A*01:01', 'A*02:01', 'A*02:07', 'A*11:01', 'A*11:03', 'A*24:02', 'A*24:07', 'A*25:01',
        'A*26:01', 'A*26:02', 'A*26:03', 'A*26:38', 'A*30:01', 'A*30:02', 'A*32:01', 'A*33:03',
        'A*68:01', 'A*68:02', 'A*74:01', 'B*07:02', 'B*13:02', 'B*15:01', 'B*15:02', 'B*15:03',
        'B*15:10', 'B*15:16', 'B*15:17', 'B*15:24', 'B*15:40', 'B*35:01', 'B*35:02', 'B*35:03',
        'B*35:05', 'B*35:08', 'B*35:11', 'B*40:01', 'B*40:02', 'B*40:06', 'B*42:01', 'B*42:02',
        'B*44:02', 'B*44:03', 'B*44:15', 'B*45:01', 'B*50:01', 'B*51:01', 'B*52:01', 'B*53:01',
        'B*57:01', 'B*57:02', 'B*57:03', 'B*58:01', 'B*58:02', 'B*67:01', 'C*01:02', 'C*04:01',
        'C*08:01', 'C*12:02', 'C*15:02'
    },
    'PtxS1': {
        'A*02:01', 'A*02:07', 'A*68:02', 'B*07:02', 'B*13:02', 'B*35:01', 'B*35:02', 'B*35:03',
        'B*35:05', 'B*35:08', 'B*35:11', 'B*40:01', 'B*40:02', 'B*40:06', 'B*42:01', 'B*42:02',
        'B*44:02', 'B*44:03', 'B*44:15', 'B*45:01', 'B*50:01', 'B*51:01', 'B*53:01', 'B*57:01',
        'B*57:02', 'B*57:03', 'B*58:01', 'B*58:02', 'B*67:01'
    },
}

# Implementation detail; see the repository documentation.

EVALUABLE_SET_MHC_II = {
    'p24': {
        'DRB1*01:01', 'DRB1*03:01', 'DRB1*04:01', 'DRB1*04:05', 'DRB1*07:01', 'DRB1*09:01',
        'DRB1*11:01', 'DRB1*13:02', 'DRB1*15:01', 'DRB3*01:01', 'DRB4*01:01', 'DRB5*01:01'
    },
    'pp65': {
        'DRB1*03:01', 'DRB1*04:01', 'DRB1*07:01', 'DRB1*11:01', 'DRB1*15:01', 'DRB3*01:01',
        'DRB3*02:02'
    },
    'PtxS1': {
        'DRB1*01:01', 'DRB1*11:01'
    },
}

FULL_REFERENCE_SET_MHC_II = {
    'p24': {
        'DRB1*01:01', 'DRB1*03:01', 'DRB1*03:02', 'DRB1*04:01', 'DRB1*04:04', 'DRB1*04:05',
        'DRB1*07:01', 'DRB1*08:01', 'DRB1*09:01', 'DRB1*10:01', 'DRB1*11:01', 'DRB1*13:01',
        'DRB1*13:02', 'DRB1*13:03', 'DRB1*13:04', 'DRB1*14:01', 'DRB1*15:01', 'DRB1*15:02',
        'DRB3*01:01', 'DRB3*03:01', 'DRB3*03:03', 'DRB4*01:01', 'DRB5*01:01'
    },
    'pp65': {
        'DRB1*03:01', 'DRB1*04:01', 'DRB1*04:02', 'DRB1*04:04', 'DRB1*07:01', 'DRB1*11:01',
        'DRB1*11:04', 'DRB1*15:01', 'DRB3*01:01', 'DRB3*02:02'
    },
    'PtxS1': {
        'DRB1*01:01', 'DRB1*11:01'
    },
}


def get_allele_sets(mhc_class: str, antigen: str) -> tuple[set[str], set[str]]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    normalized_class = str(mhc_class).strip().upper()
    by_class = {
        'I': (EVALUABLE_SET_MHC_I, FULL_REFERENCE_SET_MHC_I),
        'II': (EVALUABLE_SET_MHC_II, FULL_REFERENCE_SET_MHC_II),
    }
    if normalized_class not in by_class:
        raise ValueError(f"Processing details{mhc_class!r}Validation status")

    evaluable_by_antigen, full_by_antigen = by_class[normalized_class]
    if antigen not in evaluable_by_antigen or antigen not in full_by_antigen:
        available = sorted(set(evaluable_by_antigen) & set(full_by_antigen))
        raise KeyError(f"Allele status{antigen!r}Processing details{available}")

    evaluable = set(evaluable_by_antigen[antigen])
    full_reference = set(full_by_antigen[antigen])
    if not evaluable <= full_reference:
        extra = sorted(evaluable - full_reference)
        raise ValueError(
            f"Processing details"
            f"MHC {normalized_class}, {antigen}: {extra}"
        )
    return evaluable, full_reference


if __name__ == "__main__":
    # Validation.
    for cls_name, ev_dict, ref_dict in (
        ('MHC I', EVALUABLE_SET_MHC_I, FULL_REFERENCE_SET_MHC_I),
        ('MHC II', EVALUABLE_SET_MHC_II, FULL_REFERENCE_SET_MHC_II),
    ):
        print(f"\n{cls_name}")
        for ag in ev_dict:
            ev, ref = get_allele_sets(cls_name.replace('MHC ', ''), ag)
            assert len(ev) == len(set(ev)) and len(ref) == len(set(ref))
            print(f"  {ag:<6} evaluable={len(ev):>3}  full_reference={len(ref):>3}")
