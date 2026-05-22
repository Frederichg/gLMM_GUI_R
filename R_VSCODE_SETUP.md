# VS Code + R Setup (Windows)

This project now includes a first R-backed analyzer:
- Python GUI: glmm_analyzer_r.py
- R backend: run_analysis.R
- R package installer: install_r_packages.R

## 1) Install R

1. Install R from CRAN for Windows:
   https://cran.r-project.org/
2. During installation, keep default options.
3. Verify Rscript is available:
   - Typical path: C:/Program Files/R/R-<version>/bin/Rscript.exe

## 2) Install VS Code Extensions

Install these extensions in VS Code:
- R (REditorSupport.r)
- Python (ms-python.python)

Optional but useful:
- R Debugger (R-Debugger.r-debugger)
- Markdown All in One (yzhang.markdown-all-in-one)

## 3) Configure PATH (recommended)

Add R bin folder to PATH so the GUI can call Rscript automatically.
Example path:
C:/Program Files/R/R-4.4.2/bin

If you do not set PATH, use the GUI button "Set Rscript Path" once.

## 4) Install R Packages for this Analyzer

Run in terminal from this project folder:

```powershell
Rscript install_r_packages.R
```

If you see:

```text
'lib = "C:/Program Files/R/.../library" is not writable
```

the installer script now automatically falls back to your user library (for example under `Documents/R/win-library/<R-version>`). Re-run the same command.

If package installation fails due to build tools, install Rtools:
https://cran.r-project.org/bin/windows/Rtools/

## 5) Run the Analyzer

```powershell
python glmm_analyzer_r.py
```

Workflow:
1. Browse an Excel file.
2. Select sheet, Animal ID, IV(s) with type (within/between), and DV(s) with type.
3. Click Analyze.
4. Outputs are written under outputs/<timestamp>/.

## 6) Expected Outputs

Each run produces:
- results_summary.json
- report_*.md
- report_*.html (if rmarkdown/pandoc conversion works)
- report_*.pdf (if pandoc PDF pipeline is available)
- figures/ (DHARMa and EMMeans plots)

## 7) Notes for First Version

- This is a robust first version focused on reliable execution and diagnostics.
- If mixed-model convergence fails, the R backend falls back to GLM + cluster-robust SE (Plan B).
- EMMeans are currently generated using the first IV for maximum stability.
- PDF generation depends on your local pandoc/LaTeX availability; HTML and Markdown are fallback outputs.
