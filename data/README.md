# Data directory

The scripts expect one directory per antigen. Machine-specific source folders
have been replaced by these portable repository paths:

| Antigen | Repository folder |
|---|---|
| p24 | `data/p24` |
| pp65 | `data/pp65` |
| PtxS1 | `data/PtxS1` |

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

The repository version includes source inputs, selected intermediate outputs,
and final figures in this hierarchy. Output directories are created by the
scripts when required. `SHA256SUMS.txt` provides a complete integrity manifest
for the versioned data files; regenerate it whenever the published dataset is
intentionally changed.

Before publishing any input data, verify that redistribution is allowed and
remove personal, confidential, or credential-bearing files. For large files,
prefer an archival data repository and record its DOI in the main README.
