# Data directory

The scripts expect one directory per antigen. The original workstation folders
named `Запуски тулов на ...` have been replaced by portable English paths:

| Original folder | Repository folder |
|---|---|
| `Запуски тулов на p24` | `data/p24` |
| `Запуски тулов на pp65` | `data/pp65` |
| `Запуски тулов на PtxS1` | `data/PtxS1` |

Preserve the file names produced by each prediction tool. A shortened example
is shown below; the same pattern is used for all three antigens.

```text
data/p24/
├── p24.xlsx
├── p24.txt
├── MHC I/
│   ├── IEDB I/
│   ├── NetCTL/
│   ├── NetMHC/
│   └── NetMHCpan/
├── MHC II/
│   ├── IEDB II/
│   ├── NetMHC II/
│   └── NetMHCIIpan/
├── Matches MHC I/
├── Matches MHC II/
├── Benchmarking MHC I/
├── Benchmarking MHC II/
├── FRANK MHC I/
├── FRANK MHC II/
├── McNemar test MHC I/
├── McNemar test MHC II/
└── Visualization/
```

Only source inputs need to be copied initially. Output directories are created
by the scripts when required. Consult `project_tree.txt` in `docs/` for the full
file inventory from the original analysis workspace.

Before publishing any input data, verify that redistribution is allowed and
remove personal, confidential, or credential-bearing files. For large files,
prefer an archival data repository and record its DOI in the main README.

