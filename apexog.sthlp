{smcl}
{* *! version 0.13.2  2026-06-07}{...}
{vieweralsosee "[R] cfprobit" "help cfprobit"}{...}
{vieweralsosee "[R] cfregress" "help cfregress"}{...}
{viewerjumpto "Syntax" "apexog##syntax"}{...}
{viewerjumpto "Description" "apexog##description"}{...}
{viewerjumpto "Options" "apexog##options"}{...}
{viewerjumpto "Examples" "apexog##examples"}{...}
{viewerjumpto "Stored results" "apexog##results"}{...}
{viewerjumpto "Methods and formulas" "apexog##methods"}{...}
{viewerjumpto "References" "apexog##references"}{...}
{title:Title}

{phang}
{bf:apexog} {hline 2} APE bounds in a control-function probit under
plausibly-exogenous instruments (Conley, Hansen, & Rossi 2012)


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:apexog} {it:method} {depvar} [{indepvars}] {cmd:(}{it:endog}{cmd: =} {it:iv}{cmd:)}
{ifin}
[{cmd:,} {it:options}]

{phang}
where {it:method} is one of:

{synoptset 20 tabbed}{...}
{synopthdr:method}
{synoptline}
{synopt:{opt uci}}Union of Confidence Intervals over a grid of assumed gamma values{p_end}
{synopt:{opt ltz}}Local-to-Zero simulation from a Gaussian prior on gamma{p_end}
{synoptline}

{synoptset 32 tabbed}{...}
{synopthdr:options}
{synoptline}
{syntab:Model}
{synopt:{opt est:imator(cfprobit)}}main-equation estimator; only {cmd:cfprobit} is supported (default){p_end}
{synopt:{opt fir:ststage(}{it:type}{cmd:)}}first-stage model: {opt linear} (default), {opt probit}, {opt fprobit}, or {opt poisson}.  Accepts one type (applied to all endogs) or one per endog.{p_end}
{synopt:{opt aie(}{it:X Z}{cmd:)}}switch the sensitivity target from the APE to the Average Interaction Effect of {it:X} and {it:Z} (cts x cts, cts x binary, or binary x binary; type auto-detected).  See {it:ginteff} (Radean 2023).{p_end}
{synopt:{opt cfint:eract}}add the products of the control functions to the main probit (replicates {help cfprobit}'s {cmd:cfinteract}).  Auto-on when an interaction of an endogenous treatment with a moderator is declared as an endogenous regressor (see below).{p_end}

{syntab:UCI (method = uci)}
{synopt:{opt gmin(numlist)}}lower bound of assumed gamma range; one value per instrument (a single value is broadcast to all){p_end}
{synopt:{opt gmax(numlist)}}upper bound of assumed gamma range; one value per instrument (a single value is broadcast to all){p_end}
{synopt:{opt grid(#)}}number of grid points per instrument; default {cmd:grid(11)}{p_end}

{syntab:LTZ (method = ltz)}
{synopt:{opt mu(numlist)}}Gaussian prior mean for each instrument's gamma (single value broadcast to all){p_end}
{synopt:{opt omega(numlist)}}Gaussian prior variance for each instrument's gamma (single value broadcast to all){p_end}
{synopt:{opt iter:ations(#)}}number of Monte-Carlo draws; default {cmd:iterations(1000)}{p_end}
{synopt:{opt seed(#)}}random-number seed for reproducibility{p_end}
{synopt:{opt graphomega(numlist)}}delta sweep: rerun the LTZ simulation for each prior variance in the list (one value per step, broadcast to all instruments; Gaussian priors only) and store/plot bounds vs delta{p_end}
{synopt:{opt graphmu(numlist)}}prior means for the sweep steps (default: zeros){p_end}
{synopt:{opt graphdelta(numlist)}}x-axis labels for the sweep steps (default: the omega values){p_end}

{syntab:Inference}
{synopt:{opt gu:nits(latent|pr)}}units of the user-supplied gamma values: {opt latent} (default, latent-index units) or {opt pr} (probability points; divided by m_phi from the baseline fit){p_end}
{synopt:{opt level(#)}}confidence level; default {cmd:level(.95)}{p_end}
{synopt:{opt vce(}{it:vcetype}{cmd:)}}{it:vcetype} may be {opt robust} (default) or {opt cluster} {it:clustvar}{p_end}

{syntab:Target selection (v0.11)}
{synopt:{opt tar:get(spec list)}}selects the sensitivity target.  Three forms:{p_end}
{synopt:{cmd:    target(d1 d2)}}APE for the listed endogenous regressors (subset of varlist2){p_end}
{synopt:{cmd:    target(c.d1)}}same as {cmd:target(d1)}; the {cmd:i.}/{cmd:c.} prefix on an APE spec is stripped{p_end}
{synopt:{cmd:    target(c.X#i.d)}}AIE for the interaction; prefixes drive the formula:{p_end}
{synopt:{cmd:        c.X#c.Z}}cts × cts → cross-partial{p_end}
{synopt:{cmd:        c.X#i.d}}cts × binary → firstdiff in the binary{p_end}
{synopt:{cmd:        i.X#i.d}}binary × binary → 2nd discrete difference{p_end}
{synopt:}target() is mutually exclusive with {opt aie()}.  Default (target omitted): APE for every endog in the model.{p_end}

{syntab:Display}
{synopt:{opt noheader}}suppress the formatted output header{p_end}
{synopt:{opt gra:ph}}produce a graph of UCI bounds vs gamma (uci) or histogram of LTZ draws (ltz){p_end}
{synopt:{opt graphopt:ions(}{it:string}{cmd:)}}extra options passed to the {cmd:twoway} call (overrides defaults){p_end}
{synopt:{opt debug}}print the gmm/probit fits and Mata intermediates per iteration{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:apexog} answers the question: {it:in a control-function probit
model, how large would the exclusion-restriction violation gamma have to
be before the estimated Average Partial Effect of the endogenous regressor
crosses zero or loses significance?}  It adapts Conley, Hansen, and Rossi's
(2012) "Plausibly Exogenous" framework -- originally formulated for 2SLS
coefficients -- to target the APE (probability scale) in
{help cfprobit:cfprobit}.

{pstd}
With {opt aie(X Z)} or a {opt target()} interaction spec, the same UCI/LTZ
machinery instead targets the Average Interaction Effect (the AIE quantity
reported by ginteff; Radean 2023).  All three interaction types are
supported -- continuous x continuous (cross-partial), continuous x binary
(first difference), and binary x binary (second discrete difference) -- with
the type auto-detected from the data, or forced via the {cmd:c.}/{cmd:i.}
prefixes in {opt target()}.

{pstd}
apexog does NOT require StataNow: it implements the control-function probit
directly (regress/probit/fracreg/poisson plus Mata) and never calls
{cmd:cfprobit} at runtime; the cfprobit comparisons live in the validation
suite only.  Typing {cmd:apexog} without arguments after estimation replays
the posted baseline estimates.

{pstd}
Two methods are available:

{phang}
{cmd:uci} -- Union of Confidence Intervals.  For each gamma on a grid over
{cmd:[gmin, gmax]}, the main-equation latent index is shifted by Z*gamma
(via an {cmd:offset} term) and the APE + its delta-method standard error
are recomputed.  Reports the union of per-gamma CIs.  This is the
conservative interval: it provides correct coverage as long as the true
gamma lies somewhere in {cmd:[gmin, gmax]}.

{phang}
{cmd:ltz} -- Local-to-Zero.  Draws gamma_i from the prior (Gaussian
N({cmd:mu}, {cmd:omega}) or {opt distribution()}), recomputes the APE and its
Murphy-Topel SE for each draw, and convolves each point APE with its sampling
noise: APE*_i = APE_i + SE_i * eps_i, eps_i ~ N(0,1).  The reported interval
is the percentile interval of the CONVOLVED draws -- the simulation
implementation of CHR (2012) eq. (7), where the estimator's sampling
distribution and the exogeneity-error term A*gamma are independent and must
be added (in the linear model this is the familiar V + A*Omega*A' interval
reported by {help plausexog}).  As omega -> 0 the LTZ interval converges to
the baseline gamma = 0 confidence interval, never below it.  The prior-only
spread of the point APEs (no sampling noise) is also displayed, as a
sensitivity band, and returned in {cmd:e(lo_prior_vec)}/{cmd:e(hi_prior_vec)}.

{pstd}
Standard errors use the Murphy-Topel two-step sandwich (Murphy & Topel
1985), computed in Mata.  At gamma = 0 the resulting APE SE matches
{help cfprobit:cfprobit} + {help margins:margins, dydx()} to machine
precision (validated in {it:./_sanity/06_mt_vs_gmm.do}).


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt estimator(cfprobit)} fixes the main-equation estimator to cfprobit
(probit with control function).  Currently the only supported value.
Reserved for forward compatibility with future alternative estimators
(e.g. ivprobit).

{phang}
{opt firststage(type)} controls the first-stage regression for the
endogenous variable.  Supported types:

{phang2}{opt linear} {hline 2} OLS regression (default).  Endogenous d
is treated as continuous.{p_end}

{phang2}{opt probit} {hline 2} Probit regression.  Endogenous d must be
binary {0, 1}.  The control function is the generalized residual
{cmd:(d - Phi(Wpi))*phi(Wpi) / [Phi*(1-Phi)]}.{p_end}

{phang2}{opt fprobit} {hline 2} Fractional probit QMLE (Stata's
{help fracreg:fracreg probit}).  Endogenous d may take any value in
[0, 1].{p_end}

{phang2}{opt poisson} {hline 2} Poisson regression for count-valued
endogenous d.  The control function is the raw residual {cmd:d - exp(Wpi)}.{p_end}

{phang}
{ul:Phase-1 (v0.10) shortcut:} the type can also be declared in the
endogenous-variable clause itself via the standard Stata factor-variable
prefix:

{phang2}{cmd:(c.d = z)} -- treats d as continuous -> default {opt firststage(linear)}{p_end}
{phang2}{cmd:(i.d = z)} -- treats d as binary -> default {opt firststage(probit)}{p_end}

The {opt firststage()} option still overrides the prefix.  The prefix is
stripped before passing d to the main probit (binary-ness is then
auto-detected from the values, so the discrete-change APE is used).

{phang}
{opt aie(X Z)} switches the sensitivity target from the APE to the
Average Interaction Effect (AIE) for the interaction of {it:X} and
{it:Z}, in the style of ginteff (Radean 2023).  {cmd:apexog}
adds the interaction term {it:X} * {it:Z} to the regressor list
automatically -- do NOT also write {cmd:c.X##c.Z} or {cmd:i.X##i.Z}
in the varlist.  Both {it:X} and {it:Z} must appear among the
exogenous OR endogenous regressors.  EITHER OR BOTH may be endogenous
-- this is the canonical cfprobit use case (interaction of an
endogenous treatment with an exogenous or endogenous moderator).

Each interactor's type (continuous vs binary) is auto-detected and
the appropriate AIE formula is used:

{phang2}{ul:continuous-by-continuous} -- cross-partial:{p_end}
{p 12 12 2}AIE = mean(phi(eta) * (b_XZ - eta * a * b)){p_end}
{p 12 12 2}with a_i = b_X + b_XZ * Z_i, b_i = b_Z + b_XZ * X_i.{p_end}

{phang2}{ul:continuous-by-binary} (say {it:X} cts, {it:Z} binary)
-- firstdiff in {it:Z} of the slope in {it:X}:{p_end}
{p 12 12 2}AIE = mean( phi(eta|Z=1)*(b_X + b_XZ) - phi(eta|Z=0)*b_X ){p_end}

{phang2}{ul:binary-by-binary} -- 2nd discrete difference of Phi:{p_end}
{p 12 12 2}AIE = mean( Phi(eta|1,1) - Phi(eta|1,0) - Phi(eta|0,1) + Phi(eta|0,0) ){p_end}

When an interactor is endogenous, the {it:X} * {it:d} tempvar is
automatically EXCLUDED from that endog's first-stage regression
(otherwise the first stage would be circular).  The Murphy-Topel
correction propagates first-stage uncertainty into the AIE delta-method
SE in all three cases.  Multi-level factor (3+ unique non-{0,1} values)
interactors are not supported in v0.9.

{dlgtab:UCI}

{phang}
{opt gmin(numlist)} and {opt gmax(numlist)} {hline 2} define the assumed
range for gamma in each instrument.  Both must be specified for {cmd:uci};
each takes a {it:numlist} with one value per instrument (in the order the
instruments appear in the parenthesized clause).

{phang}
{opt grid(#)} {hline 2} number of evenly-spaced grid points per instrument.
Total fits per UCI run is {cmd:grid^k} where {it:k} is the number of
instruments.

{dlgtab:LTZ}

{phang}
{opt mu(numlist)} and {opt omega(numlist)} {hline 2} Gaussian prior on
gamma: gamma_k ~ N(mu_k, omega_k).  Both must be specified for {cmd:ltz};
{cmd:omega} is the variance (not the standard deviation), and must be
positive.

{phang}
{opt iterations(#)} {hline 2} number of simulation draws.  Default 1000.

{phang}
{opt seed(#)} {hline 2} sets the random-number seed via {cmd:set seed}
before drawing the gamma samples.  Use this for reproducibility.

{phang}
{opt graphomega(numlist)} [{opt graphmu(numlist)} {opt graphdelta(numlist)}]
{hline 2} the DELTA SWEEP: after the main LTZ results, rerun the simulation
for each prior variance in {opt graphomega()} (a single value per step,
broadcast to all instruments; Gaussian priors only) and report the
convolved CI and prior-only band per step -- the CHR (2012) Figures 1-2
"bounds as a function of the prior width" display, for the FIRST target.
Results are stored in {cmd:e(sweep)} (columns: delta mu lb ub lb_prior
ub_prior); with {opt graph}, the sweep plot replaces the LTZ histogram.
{opt graphmu()} sets per-step prior means (default zeros); {opt graphdelta()}
relabels the x-axis (default: the omega values).  {opt gunits(pr)} applies
to sweep values too.  Each step costs a full set of {opt iterations()}
draws.

{dlgtab:Inference}

{phang}
{opt gunits(latent|pr)} {hline 2} units in which the user-supplied gamma
quantities ({opt gmin()}/{opt gmax()}, {opt mu()}/{opt omega()}, and
{opt distribution()} draws) are expressed.

{phang2}{cmd:gunits(latent)} (default) {hline 2} gamma is a coefficient in the
latent index (probit standard-deviation units), i.e. the direct effect of the
instrument on y*.  This is the structural scale on which the offset enters
the model.{p_end}

{phang2}{cmd:gunits(pr)} {hline 2} gamma values are PROBABILITY POINTS: the
instrument's direct effect on Pr(y=1) per unit of Z.  Internally apexog
divides every gamma quantity by m_phi = mean(normalden(eta_hat)) from the
baseline (gamma = 0) fit -- the average dPr/dindex slope, the same factor
that maps the latent coefficient into the APE ({opt omega()}, a variance, is
divided by m_phi^2; for non-Gaussian {opt distribution()} priors the draw
itself is rescaled, which handles chi2/t/gamma/poisson/special whose
parameters are not scale parameters).  This puts apexog's gamma on the same
scale as {help plausexog}'s when the outcome is binary and plausexog is run
as a linear probability model, so the two commands can be compared for the
same elicited violation.  The conversion is first-order (exact as
gamma -> 0; m_phi is held at its baseline value).  All gamma values
displayed and stored ({cmd:e(grid)}, graphs) are in LATENT units; a header
note records the conversion, and m_phi is returned in {cmd:e(mphi)}.{p_end}

{phang}
{opt level(#)} {hline 2} confidence level (both methods; for {cmd:ltz} it is
the percentile level of the convolved draws).  Must be in (0, 1).
Default 0.95.

{phang}
{opt vce(vcetype)} {hline 2} variance-covariance estimator:

{phang2}{cmd:vce(robust)} {hline 2} heteroskedasticity-consistent (default).
The Murphy-Topel correction is applied to account for the generated
regressor.{p_end}

{phang2}{cmd:vce(cluster} {it:clustvar}{cmd:)} {hline 2} cluster-robust on
{it:clustvar}.  The cluster correction is propagated through both the
first-stage VCE and the Murphy-Topel meat (cluster-summed scores).{p_end}

{phang}
For bootstrap inference, wrap {cmd:apexog} in Stata's
{help bootstrap:bootstrap} prefix:

{phang2}
{cmd:. bootstrap, reps(500): apexog uci y X (d=Z), gmin(-0.3) gmax(0.3)}
{p_end}


{marker examples}{...}
{title:Examples}

{phang}{cmd:. use myproject_data}{p_end}

{pstd}{ul:UCI with one IV, gamma in [-0.3, 0.3]:}{p_end}
{phang}{cmd:. apexog uci grad income age (hsgpa = hscomp), gmin(-0.3) gmax(0.3)}{p_end}

{pstd}{ul:UCI with a binary endogenous regressor and probit first stage:}{p_end}
{phang}{cmd:. apexog uci employed age educ (vetstatus = draftnumber),}{p_end}
{phang}{cmd:    firststage(probit) gmin(-0.1) gmax(0.1) grid(11)}{p_end}

{pstd}{ul:LTZ with a tight Gaussian prior centered at 0:}{p_end}
{phang}{cmd:. apexog ltz grad income age (hsgpa = hscomp),}{p_end}
{phang}{cmd:    mu(0) omega(0.0025) iterations(2000) seed(12345)}{p_end}

{pstd}{ul:Cluster-robust SEs:}{p_end}
{phang}{cmd:. apexog uci y X (d = Z), gmin(-0.5) gmax(0.5)}{p_end}
{phang}{cmd:    vce(cluster firm_id)}{p_end}

{pstd}{ul:Visualize the UCI:}{p_end}
{phang}{cmd:. apexog uci y X (d = Z), gmin(-0.5) gmax(0.5) grid(21) graph}{p_end}

{pstd}{ul:Phase-1 syntax with i. prefix and target() narrowing (v0.10):}{p_end}
{phang}{cmd:. apexog uci y X (i.d = Z), gmin(-0.1) gmax(0.1)} {it:// = firststage(probit)}{p_end}
{phang}{cmd:. apexog uci y X (d1 d2 = Z1 Z2), gmin(0 0) gmax(0 0) target(d1)} {it:// only d1's APE is reported}{p_end}

{pstd}{ul:Phase-2 syntax: target(c.X#i.d) for AIE with prefix-driven type (v0.11):}{p_end}
{phang}{cmd:. apexog uci y X Z (d = Z2), gmin(-0.3) gmax(0.3) target(c.X#c.Z)} {it:// cts × cts AIE}{p_end}
{phang}{cmd:. apexog uci y X (i.d = Z), gmin(-0.3) gmax(0.3) target(c.X#i.d)} {it:// cts × bin AIE (firstdiff)}{p_end}
{phang}{cmd:. apexog uci y X (i.d = Z), gmin(-0.3) gmax(0.3) target(i.X#i.d)} {it:// bin × bin AIE}{p_end}

{pstd}{ul:Sensitivity of the Average Interaction Effect for two exogenous
moderators under plausible exclusion-restriction violations:}{p_end}
{phang}{cmd:. apexog uci y X Z (d = Z2), gmin(-0.3) gmax(0.3) grid(11) aie(X Z)}{p_end}
{phang}{cmd:. apexog ltz y X Z (d = Z2), mu(0) omega(0.0025) iterations(2000) aie(X Z)}{p_end}

{pstd}{ul:Canonical case (v0.13): interaction of an ENDOGENOUS treatment d with
a moderator X.  Because d*X is itself endogenous, declare it as a SECOND
endogenous regressor and instrument it (typically with the base instrument(s)
and their moderator interactions), exactly as in {help cfprobit} / {help plausexog}.
{cmd:cfinteract} is auto-on and the AIE is reported on the probability scale:}{p_end}
{phang}{cmd:. apexog uci y X (c.d c.d#c.X = Zinst c.Zinst#c.X),}{p_end}
{phang}{cmd:      target(c.d#c.X) gmin(-0.3 -0.3) gmax(0.3 0.3) grid(11)}{p_end}

{pstd}{it:Note:} the older {cmd:aie(d X)} form (with d in a single endog clause
{cmd:(d = Zinst)}) treats {cmd:d*X} as an exogenous control -- a mis-specification
when d is endogenous.  Prefer the explicit two-endog form above.{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:apexog} is an estimation command: it posts {cmd:e(b)} (the baseline,
gamma = 0, target estimates -- one APE per endogenous regressor, or the
AIE), {cmd:e(V)} (their full Murphy-Topel delta-method covariance,
cross-target terms included), and {cmd:e(sample)}.  {cmd:test},
{cmd:lincom}, {cmd:estimates store}, {cmd:... if e(sample)}, and replay
all work.  In addition it stores (common to both methods):

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}sample size{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}
{synopt:{cmd:e(baseline_ape)}}APE (or AIE if {opt aie()}) at gamma = 0{p_end}
{synopt:{cmd:e(baseline_se)}}MT-corrected SE at gamma = 0{p_end}
{synopt:{cmd:e(mphi)}}mean(normalden(eta_hat)) at the baseline fit: average dPr/dindex slope; the latent-to-probability conversion factor used by {opt gunits(pr)}{p_end}
{synopt:{cmd:e(lb_}{it:endog}{cmd:)}}lower bound of the APE interval (per endog){p_end}
{synopt:{cmd:e(ub_}{it:endog}{cmd:)}}upper bound of the APE interval (per endog){p_end}
{synopt:{cmd:e(lb_aie)} / {cmd:e(ub_aie)}}lower/upper bound of the AIE interval (with {opt aie()}){p_end}

{p2col 5 20 24 2: Macros (additional in {opt aie()} mode)}{p_end}
{synopt:{cmd:e(target)}}{cmd:ape} (default) or {cmd:aie}{p_end}
{synopt:{cmd:e(aie_vars)}}the two variables in {opt aie(X Z)}, space-separated{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:apexog}{p_end}
{synopt:{cmd:e(method)}}{cmd:uci} or {cmd:ltz}{p_end}
{synopt:{cmd:e(gunits)}}{cmd:latent} or {cmd:pr} (units of the user-supplied gamma values){p_end}
{synopt:{cmd:e(estimator)}}{cmd:cfprobit}{p_end}
{synopt:{cmd:e(firststage)}}{cmd:linear}, {cmd:probit}, {cmd:fprobit}, or {cmd:poisson}{p_end}
{synopt:{cmd:e(vce)}}{cmd:Murphy-Topel}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(endog)}}name of endogenous regressor{p_end}
{synopt:{cmd:e(exog)}}names of exogenous regressors{p_end}
{synopt:{cmd:e(instruments)}}names of instruments{p_end}

{pstd}
For {cmd:method = uci}, also:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(baseline_lo)} / {cmd:e(baseline_hi)}}baseline (gamma=0) CI bounds{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(grid)}}{cmd:(grid^k) x (k+3)} matrix with columns (g1, ..., gk, ape, lo, hi){p_end}

{pstd}
For {cmd:method = ltz}, also:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(iterations)}}number of simulation draws{p_end}

{p2col 5 20 24 2: Matrices (J = number of targets)}{p_end}
{synopt:{cmd:e(ape_sim)}}{cmd:iterations x J} per-draw point APEs/AIEs (prior only){p_end}
{synopt:{cmd:e(se_sim)}}{cmd:iterations x J} per-draw Murphy-Topel SEs{p_end}
{synopt:{cmd:e(ape_sim_cnv)}}{cmd:iterations x J} convolved draws APE_i + SE_i*eps_i (basis of the headline CI){p_end}
{synopt:{cmd:e(lo_pct_vec)} / {cmd:e(hi_pct_vec)}}{cmd:J x 1} convolved percentile bounds (= {cmd:e(lb_*)}/{cmd:e(ub_*)}){p_end}
{synopt:{cmd:e(lo_prior_vec)} / {cmd:e(hi_prior_vec)}}{cmd:J x 1} prior-only percentile bounds (sensitivity band){p_end}
{synopt:{cmd:e(sweep)}}(with {opt graphomega()}) {cmd:S x 6} delta-sweep matrix: delta, mu, lb, ub, lb_prior, ub_prior (first target){p_end}


{marker methods}{...}
{title:Methods and formulas}

{pstd}
{ul:Factor-variable instruments.}  Factor instruments (e.g.
{cmd:dummy##c.eldisadv}) are expanded and then cleaned to full rank: base
levels ({cmd:0b.}), the moderator main effect, and bare-vs-factor duplicates
that are collinear with the controls are dropped (data-based, via {cmd:_rmcoll}),
and each surviving instrument is materialised as a plain variable.  This keeps
the first-stage design matrix full rank so the Murphy-Topel VCE is well-defined.
The reported instrument count (and hence the required number of
{cmd:gmin()}/{cmd:gmax()}/{cmd:mu()}/{cmd:omega()} values, unless a single value
is broadcast) is the cleaned, full-rank count; {cmd:e(instruments)} lists the
retained terms.

{pstd}
Following Conley, Hansen, and Rossi (2012), the main equation is
augmented with a possible violation gamma of the exclusion restriction:

{p 8 8 2}
Pr(y_i = 1 | y_e, x, z) = Phi(x_i' beta_1 + y_e,i * beta_d + chat_i * beta_chat + z_i' gamma)

{pstd}
where chat is the control function from the first-stage regression of
y_e on (x, z) (an OLS residual for linear first stage, a generalized
residual for probit/fprobit/poisson).

{pstd}
For each value of gamma, {cmd:apexog} re-fits the probit with
{cmd:offset(z * gamma)} and computes:

{p 8 8 2}APE = mean(phi(eta_hat)) * beta_d_hat (continuous d){p_end}

{p 8 8 2}APE = mean(Phi(eta_hat | d=1) - Phi(eta_hat | d=0)) (binary d){p_end}

{pstd}
Standard errors use the Murphy-Topel two-step sandwich:

{p 8 8 2}V_beta_MT = H_beta^(-1) * sum_i s_total_i s_total_i' * H_beta^(-1){p_end}

{p 8 8 2}s_total_i = s_main_i + J_2pi * psi_pi_i{p_end}

{pstd}
where s_main_i is the main-equation probit score, J_2pi is the
cross-Hessian d(s_main)/dpi, and psi_pi_i is the first-stage influence
function.  The cross-Hessian propagates the uncertainty in the
first-stage parameters into beta's variance.

{pstd}
For cluster-robust SEs, the meat is replaced by the sum of OUTER products
of cluster-summed scores:

{p 8 8 2}meat_cluster = sum_c (sum_{i in c} s_total_i) (sum_{i in c} s_total_i)'{p_end}

{pstd}
{ul:LTZ interval construction.}  CHR (2012) eq. (7) approximates the
estimator's distribution as the INDEPENDENT sum of its sampling distribution
and the exogeneity-error term: betahat ~ N(beta, V) + A*gamma, gamma ~ F.
{cmd:apexog} implements the simulation form of this on the APE/AIE scale:
each draw refits the probit at gamma_i (giving the bias-adjusted point
estimate APE_i, the nonlinear analog of betahat - A*gamma_i, with
Murphy-Topel SE_i), and the sampling-noise term is then convolved in:

{p 8 8 2}APE*_i = APE_i + SE_i * eps_i,   eps_i ~ N(0,1) independent of gamma_i{p_end}

{pstd}
The reported {opt level(#)} interval is the percentile interval of
{c -(}APE*_i{c )-}.  Because the noise term is symmetric, this equals CHR's
pivot construction [APE_0 - c_(1-a/2), APE_0 - c_(a/2)] with c_p the
percentiles of the simulated deviation distribution, and in the linear model
it reproduces {help plausexog}'s analytic V + A*Omega*A' interval (validated
in {it:./_sanity/37_ltz_convolution.do}).  The prior-only percentile spread of
the point APEs {c -(}APE_i{c )-} is reported separately as a sensitivity band;
it answers "how does the bias-adjusted estimate move as gamma varies over the
prior", not "what is a confidence interval for the APE".  For multiple
targets the eps draws are independent across targets, which is exact for the
per-target marginal intervals reported.

{pstd}
{ul:Gamma scale.}  The violation enters the LATENT index:
Pr(y=1) = Phi(... + z'gamma), so gamma is in probit standard-deviation
units, whereas {help plausexog} applied to a binary outcome (an LPM) takes
gamma in probability units.  The two scales are linked by the average link
slope m_phi = mean(normalden(eta_hat)): the average direct effect of Z on
Pr(y=1) implied by a latent gamma is m_phi * gamma.  {opt gunits(pr)}
applies the conversion gamma_latent = gamma_pr / m_phi (variances divided by
m_phi^2) using the baseline fit, which makes apexog and plausexog bounds
comparable on a common probability-scale violation axis.  The conversion is
local: m_phi is evaluated at gamma = 0 and the mapping is exact in the
limit gamma -> 0.  m_phi is reported in the output header and stored in
{cmd:e(mphi)}.

{pstd}
{ul:Average Interaction Effect (with {opt aie(X Z)})}.  Three cases
based on the auto-detected types of {it:X} and {it:Z}:

{phang2}{ul:cts x cts} -- cross-partial averaged over the sample:{p_end}
{p 12 12 2}AIE = mean( phi(eta) * (b_XZ - eta * a * b) ),{p_end}
{p 12 12 2}a_i = b_X + b_XZ*Z_i, b_i = b_Z + b_XZ*X_i.{p_end}

{phang2}{ul:cts x bin} (e.g. {it:Z} binary) -- firstdiff in {it:Z} of the
{it:X}-derivative:{p_end}
{p 12 12 2}AIE = mean( phi(eta|Z=1)*(b_X+b_XZ) - phi(eta|Z=0)*b_X ){p_end}

{phang2}{ul:bin x bin} -- 2nd discrete difference of Phi:{p_end}
{p 12 12 2}AIE = mean( Phi(eta|1,1) - Phi(eta|1,0) - Phi(eta|0,1) + Phi(eta|0,0) ){p_end}

The Jacobian wrt the full beta-vector (with special corrections for
{cmd:b_X}, {cmd:b_Z}, and {cmd:b_XZ}, derived per case) is formed in
closed form and passed through the same V_beta_MT.  The delta-method
SE then incorporates the generated-regressor uncertainty from the first
stage.  Verified against hand-coded analytic references to machine
precision in all three cases (sanity 24/27/28).


{marker references}{...}
{title:References}

{phang}
Conley, T. G., C. B. Hansen, and P. E. Rossi.  2012.
Plausibly Exogenous.  {it:Review of Economics and Statistics} 94: 260-272.

{phang}
Murphy, K. M., and R. H. Topel.  1985.
Estimation and Inference in Two-Step Econometric Models.  {it:Journal of
Business & Economic Statistics} 3: 370-379.

{phang}
Wooldridge, J. M.  2015.
Control function methods in applied econometrics.  {it:Journal of Human
Resources} 50: 420-445.

{phang}
Clarke, D.  2020.
{cmd:plausexog}: Stata module for plausibly exogenous 2SLS.  Statistical
Software Components S458310, Boston College Department of Economics.
{it:[This package targets the structural coefficient; apexog adapts
the same framework to target the APE in a cfprobit model.]}

{phang}
Radean, M.  2023.
ginteff: A new command for estimating and graphing average
semielasticities and average marginal and interaction effects after a
nonlinear regression model.  {it:Stata Journal} 23: 1015-1043.
{it:[apexog's aie() option targets the same AIE quantity as
ginteff, but adds Conley-Hansen-Rossi sensitivity bounds and a
Murphy-Topel correction for the generated regressor.]}


{title:Author}

{pstd}
Jesper N. Wulff, Department of Economics and Business Economics,
Aarhus University.  {browse "mailto:jwulff@econ.au.dk":jwulff@econ.au.dk}


{title:License}

{pstd}
{cmd:apexog} is distributed under the MIT License.  Copyright (c) 2026
Jesper N. Wulff.  See
{browse "https://github.com/jespernwulff/stata-apexog/blob/main/LICENSE":LICENSE}
for the full text.


{title:Also see}

{psee}
Manual:  {bf:[R] cfprobit}, {bf:[R] cfregress}, {bf:[R] ivprobit},
{bf:[R] margins}

{psee}
Related: {help plausexog} (Clarke's 2SLS-targeting version)
