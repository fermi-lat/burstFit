---
package: burstFit
summary: gtburstfit — Bayesian Blocks segmentation of a binned burst light curve, then a χ² fit of a sum of Norris-type pulses plus constant background, with optional st_graph plots.
layer: analysis
depends_on: [optimizers, evtbin, st_app, st_stream, st_graph, tip]
test_depends_on: [facilities]
used_by: []
last_reviewed_commit: 318d1ca
---

# burstFit Architecture

**burstFit** is a small package (about 1,300 lines of C++) that contains:

1. A static **library** with two classes. `BurstModel` is a pulse-shape
   model built on `optimizers::Function`. `NegativeStat` is a sign-flipping
   `Statistic` wrapper.
2. The **`gtburstfit`** tool. It reads a binned light curve (for example
   `gtbin algorithm=LC` output), segments it with Bayesian Blocks
   (`evtbin::BayesianBinner`), fits pulses to it with Minuit, prints the
   results and optionally plots them through `st_graph`.

`gtburstfit` writes **no output file**. All results go to the terminal
(and the plot window).

The Doxygen mainpage (`src/mainpage.h`) and `doc/release.notes` were last
updated in 2015. The mainpage's pulse formula is missing its minus sign
(see *Key abstractions*), and the mainpage does not mention `fitguess`,
`plotres` or the manual parameters. Use `pfiles/gtburstfit.par` and
`src/gtburstfit/gtburstfit.cxx` as the source of truth.

## Purpose and scope

**Does:**
- Read a light curve table: `TIME` and/or a cell-size column (`TIMEDEL`
  or `CELLSIZE`) or `TIMEDEL` keyword, and counts from `COUNTS` or `SUM`.
- Segment the light curve into Bayesian Blocks. The algorithm is in
  `evtbin` (`BayesianBinner`), not here.
- Guess an initial pulse model from peaks and valleys in the block
  profile, or take a single pulse from user parameters.
- Fit the model to the **original bins** (not the blocks) by minimizing
  Pearson χ² with `optimizers::Minuit`.
- Plot data, blocks, fit, markers and residuals, and let the user place
  markers by right-clicking in GUI mode to define a new starting model.

**Does not:**
- Bin events. Use `gtbin` (evtbin) to make the light curve.
- Fit unbinned (time-tagged) data. Event lists are accepted and blocked,
  but the fit is skipped with a warning.
- Write FITS or text output, or report parameter uncertainties.
- Provide a Python interface. The separate `BayesianBlocks` package has
  its own, unrelated Bayesian Blocks implementation with SWIG bindings
  (see `src/BayesianBlocks/doc/architecture/`).

## Position in the system

| Direction | Package | Why |
|-----------|---------|-----|
| Depends on | `optimizers` | `BurstModel` derives from `Function`; `NegativeStat` derives from `Statistic`; the tool uses `ChiSq`, `Minuit`, `Parameter`, `dArg`. Linked `PUBLIC` because the installed headers include `optimizers/Function.h` and `Statistic.h`. |
| Depends on | `evtbin` | `BayesianBinner` (blocks) and `Hist1D` (block-averaged profile, input to `BurstModel(const Hist1D*)`). |
| Depends on | `st_app` (and `hoops`, `ape` through it) | Tool framework, `AppParGroup`, `StAppGui` for the marker GUI. |
| Depends on | `st_graph` | All plotting: `Engine`, `IFrame`, `IPlot`, `Sequence`, markers. burstFit is the only real plotting client of `st_graph` (see `src/st_graph/doc/architecture/`). |
| Depends on | `st_stream` | `StreamFormatter` output; `BurstModel::write`. |
| Depends on | `tip` | Reads the input table. |
| Test only | `facilities` | `commonUtilities::getDataPath("burstFit")` in `test_burstFit`. Reached transitively (not listed in CMake). |
| Link only | `st_facilities` | Linked `PRIVATE` into the library, `gtburstfit` and the test. Nothing calls it: the test includes `st_facilities/Env.h` and has `using namespace st_facilities;` but uses no symbol. |

**CMake vs. use.** The library's sources (`BurstModel.cxx`,
`NegativeStat.cxx`) use only `optimizers`, `evtbin` and `st_stream`. The
`PRIVATE` links of the library to `st_app`, `st_graph`, `tip` and
`st_facilities` are needed only by the executables, which link them again
themselves.

**Consumers.** None. No package in `src/` includes `burstFit/` headers or
calls `gtburstfit`. `sane`'s `gt_apps` does not wrap it, and
`pyBurstAnalysisGUI` (`gtburst`) does its own analysis.

## Key abstractions

```
optimizers::Function                optimizers::Statistic
        ^                                   ^
        |                                   |
  burstFit::BurstModel   <--model--  optimizers::ChiSq(x, y, model)
                                            ^ wraps (not owned)
                                  burstFit::NegativeStat  --> optimizers::Minuit
```

| Class | Header | Role |
|-------|--------|------|
| `BurstModel` | `burstFit/BurstModel.h` | Sum of *N* pulses plus a constant. Parameters are named `Amp_i`, `Time0_i`, `Tau1_i`, `Tau2_i` (i = 0…N−1) and `Bckgnd`, stored 4 per pulse then background last (`Parameter_e`). `getNumPeaks()`, `getCoefficient(i, "Amp"/"Time0"/"Tau1"/"Tau2")`. Two constructors: from an explicit parameter vector (4N+1 entries, then `setBounds`), or from a `Hist1D` (automatic guess; see gotchas). |
| `NegativeStat` | `burstFit/NegativeStat.h` | Wraps a `Statistic*` and returns `−value()` and negated derivatives. Needed because `optimizers` optimizers **maximize** `Statistic::value()` while `ChiSq::value()` is +χ² (see `src/optimizers/doc/architecture/`). Does not own the wrapped statistic; `clone()` shares the pointer. |

### The pulse shape

For each pulse, with `t = x − Time0`:

```
f(t) = Amp · exp(2·sqrt(Tau1/Tau2) − Tau1/t − t/Tau2)    for t > 0, else 0
model(x) = Σ f_i + Bckgnd
```

This is the Norris et al. (2005) pulse, normalized so that `Amp` is the
peak height above background. The peak is at `t = sqrt(Tau1·Tau2)`. The
code (`BurstModel::value`) follows the IDL routine `dofit_prodexpback.pro`.
Both the class comment in `BurstModel.h` and `mainpage.h` give the
exponent with the wrong sign or form. Analytic derivatives
(`derivByParamImp`) are implemented for all parameters. They were checked
against the formula above when this doc was written.

### Automatic guess (`fitguess=AUTO`/`CALC`)

`BurstModel::findPeaks` is based on `findpvs3.pro`. It works on the
block-averaged `Hist1D`:
- The background is the mean of the first and last block. The code throws
  if either of those is the maximum, or if there are fewer than 3 blocks
  or no local maximum.
- Candidate peaks are local maxima.
- Valleys are minima between candidates that lie below
  `0.9·(max − background)`, plus virtual valleys at the first rise and
  last decay (5% of the mean peak height).
- The final peaks are the maxima between valleys.

`guessInitialParameters` then sets, for each peak: `Amp` = height above
background; `Time0` = start of the first block after the valley that is
above the first block's value; `Tau1 = Tau2 = peak midpoint − Time0`. It
also sets bounds (`Amp` ∈ [0, 2·Amp], `Tau` ∈ [0, 10·Tau], `Bckgnd` ∈
[0, 3·Bckgnd]). **At present this guess never produces any pulses**; see
*Invariants and gotchas*.

## Data flow and control flow

```
BurstFitApp::run   (src/gtburstfit/gtburstfit.cxx; everything is in this one file)
  prompt()  evfile, fitguess, [amp,time0,tau1,tau2,bckgnd if MAN], evtable, fit, plot, plotres, ncpprior
  tip readTable(evfile, evtable)
  detect columns: TIME? ; TIMEDEL | CELLSIZE column? ; else TIMEDEL keyword ; COUNTS | SUM (else 1 per row)
  build domain[i] (cell start), intervals[i], cell_pop[i]
    (drop the last row if there is no cell size)
  if fitguess != MAN:
    BayesianBinner(intervals, cell_pop, ncpprior)           // evtbin
    Hist1D(bb) filled with cell_pop·cellwidth/blockwidth    // mean counts per original cell, per block
  model = MAN ? BurstModel(5 user pars) : BurstModel(&hist)
  ChiSq(domain, cell_pop, model) → NegativeStat → Minuit::find_min()
    (only if fit=yes AND a cell-size column exists; exceptions printed, not rethrown)
  print χ²/dof and parameters before and after; print block table
  if plot: st_graph Engine → data (black), blocks (red), fit (green),
           3 markers per pulse (start, peak, decay), optional residuals frame; engine.run()
```

**GUI mode** (`gui=yes`, `BurstFitGui : StAppGui`). Each right-click on
the plot adds a marker, cycling through *Pulse start* → *Peak* →
*Decay*. When the user runs again, `runApp` converts every three markers
into one pulse:
- `Time0` = start x
- `Amp` = |peak y − Bckgnd|
- `Tau2` = |decay x − Time0|
- `Tau1` = (peak x − Time0)² / Tau2

It replaces `m_model` with the result and calls `run()` again. Because
`run()` creates a new model only if none exists, `fitguess=AUTO` keeps
this marker-built model, while `CALC` discards it and re-guesses from the
blocks. On the command line, `AUTO` and `CALC` behave the same.

## Public interfaces

### Applications

| Tool | Par file | Notes |
|------|----------|-------|
| `gtburstfit` | `pfiles/gtburstfit.par` | Asked: `evfile`, `fitguess` ∈ {`AUTO`, `MAN`, `CALC`}, and `amp`, `time0`, `tau1`, `tau2`, `bckgnd` (prompted only for `MAN`). Hidden: `evtable=RATE`, `fit=yes`, `plotres=no`, `ncpprior=9`, `plot=yes`, `title` (**unused**: the code builds its own title), and the standard `chatter`/`clobber`/`debug`/`gui`/`mode` (see `src/st_app/doc/architecture/`). |
| `test_burstFit` | `pfiles/test_burstFit.par` | Test program; installed to `bin/` like the tool. |

User-facing help text lives in `fermitools-fhelp`.

**Exit status.** The `st_app` contract applies (see
`src/st_app/doc/architecture/README.md`). A missing input file or a
failed automatic guess (for example "Background is as large or larger
than any peaks") exits 1. A failed **fit** or **plot** only prints a
message, and the tool exits 0.

### C++ library

`burstFit/BurstModel.h`, `burstFit/NegativeStat.h`. No other package uses
them. The library is `STATIC`.

### Python

None.

### Environment variables

No `getenv` calls. Par files are found through `PFILES`, and the test's
data through `FERMI_INST_DIR` (see "Runtime environment" in ScienceTools'
`docs/architecture/README.md`). Plotting needs a working Tk display for
`st_graph`.

### File formats

- **Input:** a FITS table, by default the `RATE` extension of a `gtbin`
  `LC` file (`TIME`, `TIMEDEL`, `COUNTS`). It also accepts old-style
  `CELLSIZE`/`SUM` tables (the `data/binData*.fits` test files, whose
  extension is unnamed, so they need `evtable=1`), and FT1/TTE event lists
  (`evtable=EVENTS`; each row counts 1 and no fit is done).
- **Output:** terminal text only.

## Build and test

- **CMake targets:** `burstFit` (STATIC library), `gtburstfit`,
  `test_burstFit`.
- **Install:**
  - headers go to `include/burstFit`
  - `pfiles/` goes to `syspfiles`
  - `data/` goes to `${FERMI_INSTALL_REFDATADIR}/burstFit`
    (`share/fermitools/refdata/fermi/burstFit`), which is where
    `getDataPath("burstFit")` looks
  - there is no `data/outref`
- **`test_burstFit`:**
  - It is an `st_app` application, so `PFILES` must be set. It is **not
    registered with CTest**. ScienceTools' `recipe/tests/ST-unit-test`
    runs it in both test lists. It has no `check_test_burstFit_output`,
    because the test writes no files.
  - It reads 18 files from `getDataPath("burstFit")`: `ft1tiny.fits`,
    `n0_p05.lc` and 16 `binData*.fits` (BATSE-style `CELLSIZE`/`SUM`
    tables). For each one, it repeats `gtburstfit`'s logic (copied, not
    shared) and fits.
  - It contains **no assertions**. It fails only if something throws,
    and fit exceptions are caught and logged. A regression in the fit
    results would go unnoticed. `plot=no` by default; `debug=yes` in its
    par file.
  - It writes nothing to disk.
- **Legacy:** `SConscript`, `burstFitLib.py` (SCons) and
  `azure-pipelines.yml` are no longer used.

## Invariants and gotchas

**Automatic guess is broken**
- **`fitguess=AUTO`/`CALC` produce a background-only model (confirmed at
  runtime with Fermitools 2.5.1 on `data/n0_p05.lc`).**
  `BurstModel(const Hist1D*)` calls `guessInitialParameters` *before*
  `findPeaks`. At that point `m_peak_index` is empty, so the guess sizes
  the parameter vector to 4·0 + 1 and creates only `Bckgnd`. The fit then
  fits a constant, which hits its `3·Bckgnd` upper bound (on `n0_p05.lc`:
  11.75 → 35.25, reduced χ² 274 → 69). No pulse markers are drawn. Commit
  `a84b4cf` (2015-03-03, "propagate changes from optimizers::Function
  refactoring") reversed the order; before that `findPeaks` ran first.
  `test_burstFit` cannot detect this because it has no assertions. As a
  workaround, use `fitguess=MAN` or GUI markers. A manual single-pulse run
  on the same file reaches reduced χ² ≈ 2.1.

**Input interpretation**
- When a `TIME` column exists, each cell's interval is **[TIME[i],
  TIME[i+1])**. `TIMEDEL` is used only for the last row. Because `gtbin`
  writes the bin **midpoint** to `TIME`, cells (and plotted bars, which use
  `TIME` as the lower edge) are shifted by half a bin. A gap in the light
  curve is absorbed into the cell before it.
- The model is evaluated and fitted at `TIME[i]` (the bin midpoint for
  `gtbin` output), against **counts per bin**, not rate. Fitted `Amp` and
  `Bckgnd` are therefore in counts per original bin.
- Fitting needs a cell-size **column** (`TIMEDEL`/`CELLSIZE`). If only a
  `TIMEDEL` *keyword* is present, the cells are built from it, but the fit
  is skipped with the misleading warning "not supported for time-tagged
  event data".
- With no cell size at all (event lists), the last row is dropped, each
  row counts 1, and each event's cell runs to the next event's time.

**Fit statistic**
- `ChiSq` is Pearson χ² with the **model** in the denominator, so
  `BurstModel::value` clamps 0 to `epsilon` to avoid dividing by 0. The
  Minuit `ErrorDef` is 0.5 (log-likelihood convention). That is wrong for
  χ², but `gtburstfit` reports no errors, so it affects only Minuit's
  convergence tolerance (see `src/optimizers/doc/architecture/`).
- Parameter bounds are relative to the starting values (`Amp` ≤ 2×,
  `Tau` ≤ 10×, `Bckgnd` ≤ 3×, `Time0` within ±3·`Tau2` for manual input).
  A poor starting guess therefore caps the fit. A starting value of 0
  gives bounds [0, 0], except for `Amp`, which gets a 2·epsilon window.
- `NegativeStat` mirrors the wrapped statistic's parameters at
  construction. Build it *after* the model's parameters are final.

**Blocks**
- `ncpprior` is not the whole prior. `BayesianBinner::computeBlocks`
  subtracts `ncp_prior + 6.7` per block (`src/evtbin/src/BayesianBinner.cxx`),
  so the default 9 means an effective penalty of 15.7.
- The algorithm is O(N²) in the number of input cells. That is fine for
  light curves, but slow for long event lists.
- The "Average Counts" printed per block is the mean counts per
  **original cell** in that block, weighted by cell width.

**Error handling**
- Fit and plot failures are printed and swallowed, so the tool exits 0.
  Scripts cannot detect a failed fit from the exit status.

## Extension points

- **Different pulse shape:** subclass `optimizers::Function` as
  `BurstModel` does: name the parameters `<coeff>_<i>`, and implement
  `value`, `derivByParamImp` and `clone`. Then construct it in
  `BurstFitApp::run` in place of `BurstModel`. The marker GUI and the
  marker drawing assume the 4-parameter pulse and the `Parameter_e` layout.
- **Different statistic:** any `optimizers::Statistic` can replace
  `ChiSq`. Wrap it in `NegativeStat` if its `value()` is something to
  minimize.
- **Writing results:** there is no output file. Adding one means adding a
  par-file parameter and a `tip` writer in `gtburstfit.cxx`.
- If you touch `BurstModel`, also add real assertions to `test_burstFit`,
  for example the number of peaks found in `n0_p05.lc` and a χ² threshold.

## Further reading

- `src/mainpage.h`: the original parameter descriptions. It has the wrong
  sign in the pulse formula and is missing newer parameters.
- `src/evtbin/doc/architecture/`: `BayesianBinner` and `Hist1D`.
- `src/optimizers/doc/architecture/`: the `Statistic` sign convention,
  `ChiSq`, `Minuit`.
- `src/st_graph/doc/architecture/`: the plotting engine and markers.
- J. D. Scargle, *Studies in Astronomical Time Series Analysis V/VI*
  (Bayesian Blocks).
- J. P. Norris et al. 1996, ApJ 459, 393, and 2005, ApJ 627, 324 (pulse
  model).
