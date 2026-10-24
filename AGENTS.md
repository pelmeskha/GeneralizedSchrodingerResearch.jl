# AGENTS.md

Guidance for agents working in this repository.


### LaTeX documents toolchain

The user compiles with MiKTeX 22.10 (pdfTeX) on Windows, through TeXworks. There is no LaTeX install inside WSL, so compile with the MiKTeX binary via WSL interop:

```
/mnt/c/Users/vikto/AppData/Local/Programs/MiKTeX/miktex/bin/x64/pdflatex.exe
```

### Julia

Julia 1.12.6 (installed on Windows through juliaup). Run it from WSL with the juliaup shim:

```
/mnt/c/Users/vikto/AppData/Local/Microsoft/WindowsApps/julia.exe
```

Run scripts with the repo project using a Windows-style path:

```
"/mnt/c/Users/vikto/AppData/Local/Microsoft/WindowsApps/julia.exe" --project=C:/Users/vikto/Desktop/Projects/GeneralizedSchrodingerResearch.jl Coupled_GL/run/run_0_main.jl
```

Repository layout for Julia code:
- `src/` is the shared package `GeneralizedSchrodingerResearch`. Its entry `src/GeneralizedSchrodingerResearch.jl` loads the shared modules `BaseUtilities` (`src/base_utilities/`), `Utilities` (`src/utilities/`: `periodic_time`, `field_error`) and `Plotting` (`src/plotting/`), plus the two submodules below.
- `Coupled_GL/src/CoupledGL.jl` is submodule `GeneralizedSchrodingerResearch.CoupledGL`: analytical solutions, residuals, the coupled solver, `sn_reference()` and convergence helpers. `Coupled_GL/run/` holds its runnable scripts named `run_N_<name>.jl`. Their outputs go to `Coupled_GL/results/` (`RESULTS_DIR`).
- `Coupled_GL/LaTeX/` holds the Coupled GL documents, including the PHD git submodule.
- `GNSE_3_5_7/src/GNSE_3_5_7.jl` is submodule `GeneralizedSchrodingerResearch.GNSE_3_5_7`: NSE solutions, `solve.jl`, construct helpers, and stubs for `cuda_solve` and `live_plot_record`. `GNSE_3_5_7/run/` holds its runnable scripts; their outputs go to a flat `GNSE_3_5_7/results/` folder.
- `gpu/` is a separate package `GeneralizedSchrodingerResearchGPU` (`gpu/Project.toml`, `gpu/Manifest.toml`, `gpu/src/`). It depends on CUDA and GLMakie and implements the GPU solver `cuda_solve` and the live plot `live_plot_record`. Its manifest is the only one that records CUDA and GLMakie, because Julia 1.12 does not record weak dependencies. Scripts that call `cuda_solve` need `using GeneralizedSchrodingerResearchGPU` and must run with `--project=C:/Users/vikto/Desktop/Projects/GeneralizedSchrodingerResearch.jl/gpu`. The main project must not use CUDA.
- Scripts never `include` other code. They use `using GeneralizedSchrodingerResearch...` only. Run scripts from the repo root, using `Coupled_GL/run/...` as the script path.

Do not use `/mnt/c/Users/vikto/AppData/Local/Programs/Julia-1.9.3`, which is an older install.

Do not run scripts from a WSL path such as `/home/...` as `--project`. Windows Julia does not resolve it. Use `C:/...` paths for `--project`, and relative paths for `include`.

Test simulations in a temporary copy of the repo under `/mnt/c/Users/vikto/AppData/Local/Temp/`, not in `results/`.

### Python env sandbox
conda env: coupled_GL. Sympy installed.

### How to compile

1. Do not run `pdflatex.exe` from a WSL path (`/home/...`) or a `\\wsl.localhost\...` path. It fails with a fatal error. Build inside a Windows-side folder under `/mnt/c/`, such as `/mnt/c/Users/vikto/AppData/Local/Temp/<name>`.
2. Copy the `.tex` file and every figure it references into that folder.
3. Run pdflatex **twice** (the second run resolves the table of contents and references), with `-interaction=nonstopmode` and `-shell-escape`. The template loads `minted`, which needs shell-escape:
```
P="/mnt/c/Users/vikto/AppData/Local/Programs/MiKTeX/miktex/bin/x64/pdflatex.exe"
"$P" -interaction=nonstopmode -shell-escape disser.tex
"$P" -interaction=nonstopmode -shell-escape disser.tex
```
Run each command with the build folder as `cwd`.
4. Check for errors with `grep -a -n "^!" disser.log`. Empty output means no LaTeX errors. The MiKTeX log is at `C:\Users\vikto\AppData\Local\MiKTeX\miktex\log\pdflatex.log` (`/mnt/c/Users/vikto/AppData/Local/MiKTeX/miktex/log/pdflatex.log` from WSL).
5. Delete the temporary folder when finished.

### Rules

- Build in a temporary folder, not in the document folder. The user's `disser.pdf`, `disser.aux`, `disser.log` and `disser.toc` are the user's own build output. Do not overwrite them unless the user asks.
- Do not change the encoding of `disser.tex`. Verify edits with `iconv -f CP1251 -t UTF-8 disser.tex`.
- Do not remove or change imports in `.tex` files (`\usepackage` lines and similar) without asking the user first. Some packages are needed even when no visible command uses them, such as `float`.
- Place every figure with the `[H]` flag, as `\begin{figure}[H]`. This needs the `float` package.
- Do not repeat in the body text what a figure caption already states (parameters, grid, fit values). Refer to the figure instead.
- Verify a change by compiling it and checking the log for `^!` errors before reporting it done.
- Do not install TeX packages or run `sudo` without asking the user.
- No PDF renderer was found for this setup (Ghostscript was not located), so the PDF's visual layout cannot be inspected directly. Use the `.aux`/`.toc`/`.log` files for checks.
