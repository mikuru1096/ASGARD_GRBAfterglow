# ASGARD ApJ manuscript

This directory contains an AASTeX manuscript draft for an ASGARD code and methods paper.

## Baseline

- Code baseline: `2cf7eaecda0b`
- Scope: ASGARD GRB afterglow public API, physical contracts, numerical kernels, validation artifacts, and explicitly unsupported backend boundaries.
- Figure backend: Python/matplotlib only.
- LaTeX backend: Windows-side TeX Live through `latexmk`.

## Rebuild

Generate the figures and source data from the repository root:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/generate_paper_figures.py'
```

Compile the manuscript from the repository root on Windows:

```powershell
rtk latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=paper/build paper/main.tex
```

## Manuscript evidence rules

- Main figures use generated `paper/source_data/*.csv` files.
- Existing tracked benchmark artifacts are used only when their source CSV is tracked.
- Untracked `output/asgard_doc/**` files are not manuscript evidence.
- Unsupported backend boundaries are stated as limitations, not hidden fallback behavior.
