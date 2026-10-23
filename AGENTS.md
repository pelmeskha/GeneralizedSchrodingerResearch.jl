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
"/mnt/c/Users/vikto/AppData/Local/Microsoft/WindowsApps/julia.exe" --project=C:/Users/vikto/Desktop/Projects/GeneralizedSchrodingerResearch.jl run/Coupled_GL/run_0.jl
```

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
