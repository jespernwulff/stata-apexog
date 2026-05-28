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

`apexog` requires **StataNow** (Stata 18.5 or later) because it builds on the `cfprobit` family.

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

---

## Features

- **Methods**: `uci` (Union of Confidence Intervals over a γ grid) and `ltz` (Local-to-Zero Monte Carlo).
- **First stages**: `linear`, `probit`, `fprobit`, `poisson` — one per endogenous regressor, or one shared.
- **Multiple endogenous regressors** with per-endog instrument sets (multi-clause syntax `(d1 = Z1 Z2)(d2 = Z2 Z3)`).
- **Target selection** via `target()`:
  - `target(d)`, `target(c.d)` — APE for that endog
  - `target(d1 d2)` — APE for both
  - `target(c.X#i.d)` — AIE for the interaction (prefixes drive the formula)
- **Auto-detection** of binary endogenous regressors → discrete-change APE.
- **LTZ priors**: Gaussian via `mu()` + `omega()`, or non-Gaussian via `distribution(normal | uniform | chi2 | t | gamma | poisson | special)`.
- **Inference**: Murphy-Topel two-step sandwich; `vce(robust)` (default), `vce(cluster <id>)`.
- **Graphing**: line plot of APE vs γ (`uci`) or histogram of simulated APEs (`ltz`).

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
