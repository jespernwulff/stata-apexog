# apexog

**Sensitivity analysis for the Average Partial Effect (APE) and Average Interaction Effect (AIE) in a control-function probit, under plausibly-exogenous instruments.**

`apexog` extends Conley, Hansen & Rossi's (2012) "Plausibly Exogenous" framework — originally implemented for `ivregress 2sls` in Damian Clarke's [`plausexog`](https://github.com/damiancclarke/plausexog) — to nonlinear control-function settings (Stata's `cfprobit`).  It answers the question:

> *How much exclusion-restriction violation γ ≠ 0 would it take to bring the APE (or AIE) down to zero, or to lose significance?*

Inference uses a closed-form **Murphy-Topel two-step sandwich** computed in Mata; at γ = 0 the resulting APE SE matches `cfprobit` + `margins, dydx()` to machine precision.

---

## Installation

```stata
net install apexog, from("https://raw.githubusercontent.com/jespernwulff/stata-apexog/main/")
```

`apexog` requires **Stata 18** or later. **StataNow is not required**: apexog implements the control-function probit directly (`regress`/`probit`/`fracreg`/`poisson` plus Mata) and never calls `cfprobit` at runtime — the `cfprobit` comparisons live only in the validation suite.

To uninstall:

```stata
ado uninstall apexog
```

---

## Quick example

```stata
* simulate a control-function probit DGP
clear
set seed 12345
set obs 5000
gen double X     = rnormal()
gen double Zinst = rnormal()
gen double u     = rnormal()
gen double e     = 0.5*u + sqrt(1-0.5^2)*rnormal()
gen double d     = 0.8*Zinst + 0.4*X + u             // d is endog (corr(u,e)=0.5)
gen byte    y     = (-0.2 + 0.4*X + 0.5*d + e) > 0

* APE of d, with the exclusion restriction allowed to flex over gamma in [-0.3, 0.3]
apexog uci y X (c.d = Zinst), gmin(-0.3) gmax(0.3) grid(11) target(c.d)
```

You'll get a `Baseline (gamma = 0)` row (the cfprobit-style APE), followed by the **Union-of-CIs bounds** for the APE across the assumed γ range.

For the **Average Interaction Effect** of, say, a continuous moderator × a binary endogenous treatment:

```stata
apexog uci y X (i.d = Zinst), gmin(-0.3) gmax(0.3) grid(11) target(c.X#i.d)
```

The `i.` / `c.` prefixes drive both the first-stage default (`(i.d=...)` → `firststage(probit)`) and the AIE formula choice:

- `target(c.X#c.Z)` — cts × cts (cross-partial)
- `target(c.X#i.d)` — cts × binary (firstdiff in the binary)
- `target(i.X#i.d)` — binary × binary (2nd discrete difference)

When the **interaction itself is endogenous** — the canonical case of an endogenous treatment `d` interacted with a moderator `X`, so `d#X` is endogenous too — declare the interaction as a *second endogenous regressor* and instrument it (typically with the base instruments and their moderator interactions), exactly as in `cfprobit` / `plausexog`. The control-function product (`cfinteract`) is added automatically, and the AIE is reported on the probability scale:

```stata
apexog uci y X (c.d c.d#c.X = Zinst c.Zinst#c.X), target(c.d#c.X) gmin(-0.3) gmax(0.3) grid(11)
```

Factor-variable instruments (e.g. `dummy##c.eldisadv`) are supported: they are expanded and cleaned to full rank (base levels, the moderator main effect, and bare-vs-factor duplicates are dropped) so the Murphy-Topel variance is well-defined. A single `gmin`/`gmax`/`mu`/`omega` value is broadcast to all instruments.

---

## Features

- **Methods**: `uci` (Union of Confidence Intervals over a γ grid) and `ltz` (Local-to-Zero Monte Carlo). The LTZ interval convolves the estimator's sampling noise with the γ prior per draw (`APE*_i = APE_i + SE_i·ε_i` — the simulation analog of `plausexog`'s analytic `V + AΩA′`, CHR 2012 eq. 7); the prior-only spread of the point APEs is also reported, as a sensitivity band (v0.14.0).
- **First stages**: `linear`, `probit`, `fprobit`, `poisson` — one per endogenous regressor, or one shared.
- **Multiple endogenous regressors** with per-endog instrument sets (multi-clause syntax `(d1 = Z1 Z2)(d2 = Z2 Z3)`).
- **Target selection** via `target()`:
  - `target(d)`, `target(c.d)` — APE for that endog
  - `target(d1 d2)` — APE for both
  - `target(c.X#i.d)` — AIE for the interaction (prefixes drive the formula)
- **Endogenous-interaction AIE**: declare an endogenous treatment × moderator interaction as a second instrumented endogenous regressor — `(c.d c.d#c.X = ...)` with `target(c.d#c.X)` — with the control-function product (`cfinteract`) added automatically. Matches `cfprobit` + `ginteff` to machine precision.
- **Factor-variable instruments** are expanded and cleaned to full rank (base/collinear terms dropped); a single `gmin`/`gmax`/`mu`/`omega` value broadcasts to all instruments.
- **Auto-detection** of binary endogenous regressors → discrete-change APE.
- **LTZ priors**: Gaussian via `mu()` + `omega()`, or non-Gaussian via `distribution(normal | uniform | chi2 | t | gamma | poisson | special)`.
- **Gamma units**: γ enters as an offset on the probit-index scale; `gunits(latent|pr|structural)` says which of the three elicitation scales your numbers live on. `gunits(pr)` takes **probability points** (the instrument's direct effect on Pr(y=1), the scale of `plausexog` run on a binary outcome; divided internally by m̂ = mean φ(η̂) from the baseline fit). `gunits(structural)` takes **structural-latent units** (total error variance 1, the scale simulation DGPs and textbook latent-variable models are written on; divided by σ̃ = (1 + b_cf′ V̂(ĉ) b_cf)^(−1/2), estimated from the control-function coefficients). Feeding a structural-units prior unrescaled costs coverage (2Φ(zσ̃)−1 instead of 1−α). Both factors are printed in the header and returned in `e(mphi)` / `e(stilde)` (v0.15.0, v0.17.0).
- **Inference**: Murphy-Topel two-step sandwich; `vce(robust)` (default), `vce(cluster <id>)` with a numeric or string id (`e(N_clust)`, `e(clustvar)`).
- **UCI grid guard**: `uci` stops early with an actionable message when grid^(#instruments) would exceed Stata's matrix limit — use `ltz`, fewer instruments, or a coarser grid (v0.16.1).
- **Graphing**: line plot of APE vs γ (`uci`), histogram of the convolved LTZ draws (`ltz`), or — with `graphomega()` — the CHR "bounds vs prior width" delta-sweep figure (`ltz`).
- **Estimation command**: posts `e(b)`/`e(V)` (baseline target estimates with full cross-target Murphy-Topel covariance) and `e(sample)`, so `test`, `lincom`, `estimates store`, `... if e(sample)`, and replay all work (v0.16.0).
- **Input validation** (v0.18.0): explicit errors instead of silent misbehaviour — a depvar not coded 0/1, overlapping variable lists, `gmin() > gmax()`, a first-stage model that does not suit the endogenous regressor's support, a first stage that drops observations; `level(95)` and `level(.95)` both accepted; `vce(cl id)` abbreviations; abbreviated variable names; the full command in `e(cmdline)`; internal objects cleaned up at exit. Validated by a 57-file sanity suite plus coverage Monte Carlos for cluster-robust inference, non-Gaussian priors, and multiple endogenous regressors.

---

## Documentation

After installing, see:

```stata
help apexog
```

The help file documents every option, both methods, the Mata internals (Murphy-Topel derivation), and the three AIE formulas.

---

## Citation

If you use `apexog` in published work, please cite the methodology paper and (optionally) the package:

> Conley, T. G., Hansen, C. B., & Rossi, P. E. (2012). Plausibly Exogenous. *Review of Economics and Statistics*, 94(1), 260–272. [DOI:10.1162/REST_a_00139](https://doi.org/10.1162/REST_a_00139)

```bibtex
@misc{apexog,
  author = {Wulff, Jesper},
  title  = {apexog: APE/AIE bounds in cfprobit under plausibly-exogenous IVs},
  year   = {2026},
  url    = {https://github.com/jespernwulff/stata-apexog}
}
```

The companion to the linear-2SLS version of the same framework is Damian Clarke's [`plausexog`](https://github.com/damiancclarke/plausexog).

For the AIE quantity computed by `target(c.X#i.d)` etc. when γ = 0, the reference is:

> Radean, M. (2023). ginteff: A new command for estimating and graphing average semielasticities and average marginal and interaction effects after a nonlinear regression model. *Stata Journal*, 23(4), 1015–1043.

---

## License

[MIT](LICENSE) © Jesper Wulff.
