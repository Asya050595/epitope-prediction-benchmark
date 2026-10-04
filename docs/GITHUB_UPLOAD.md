# Uploading this package to GitHub

The target repository is:
<https://github.com/Asya050595/epitope-prediction-benchmark>

## Recommended command-line method

1. Download and extract `epitope-prediction-benchmark_github-ready.zip`.
2. Open a terminal in the extracted directory.
3. If the GitHub repository is still empty except for its license, run:

```bash
git init
git branch -M main
git remote add origin https://github.com/Asya050595/epitope-prediction-benchmark.git
git add .
git commit -m "Add reproducible benchmark analysis"
git pull --rebase origin main
git push -u origin main
```

If `origin` already exists, omit the `git remote add` command. If Git reports a
conflict in `README.md`, `.gitignore`, or `LICENSE`, stop and resolve it before
continuing; do not use a force push.

## Input files

The current package contains the scripts and a data-layout guide, but not the
raw input files. Copy redistributable input files into `data/` only after
checking their licenses and privacy status. Large input files should normally be
stored in Zenodo, Figshare, an institutional repository, or Git LFS rather than
ordinary Git history.

## Before creating a public release

- Replace the account-level author entry in `CITATION.cff` with the complete
  author list and ORCID identifiers.
- Confirm that the repository license matches the license selected on GitHub.
- Add the article citation and data DOI to `README.md` when available.
- Run `python scripts/check_repository.py`.
- Create a versioned GitHub release, for example `v1.0.0`.

