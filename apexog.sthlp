{smcl}
{* *! version 0.21.0  2026-11-08}{...}
{vieweralsosee "[R] cfprobit" "help cfprobit"}{...}
{vieweralsosee "[R] cfregress" "help cfregress"}{...}
{vieweralsosee "[R] biprobit" "help biprobit"}{...}
{viewerjumpto "Syntax" "apexog##syntax"}{...}
{viewerjumpto "Description" "apexog##description"}{...}
{viewerjumpto "Options" "apexog##options"}{...}
{viewerjumpto "Examples" "apexog##examples"}{...}
{viewerjumpto "Input requirements" "apexog##limits"}{...}
{viewerjumpto "Stored results" "apexog##results"}{...}
{viewerjumpto "Methods and formulas" "apexog##methods"}{...}
{viewerjumpto "References" "apexog##references"}{...}
{title:Title}

{phang}
{bf:apexog} {hline 2} APE bounds in a control-function probit (or a recursive
bivariate probit) under plausibly-exogenous instruments (Conley, Hansen, & Rossi 2012)


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:apexog} {it:method} {depvar} [{indepvars}] {cmd:(}{it:endog}{cmd: =} {it:iv}{cmd:)}
[{cmd:(}{it:endog2}{cmd: =} {it:iv2}{cmd:)} ...]
{ifin}
[{cmd:,} {it:options}]

{phang}
{it:depvar} must be coded 0/1.  Each {it:endog} may carry a {cmd:c.} or {cmd:i.} prefix (which sets the default first stage) or be a two-way interaction such as {cmd:c.d#c.X} of an endogenous treatment with a moderator (see the canonical case under Examples).  Weights and the {cmd:svy} prefix are not supported.  Typing {cmd:apexog} alone after estimation replays the posted results; {cmd:apexog, level(90)} replays at another level.

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
{synopt:{opt est:imator(cfprobit)}}main-equation estimator; only {cmd:cfprobit} is supported (default).  {opt firststage(biprobit)} switches the whole estimator to the recursive bivariate probit.{p_end}
{synopt:{opt fir:ststage(}{it:type}{cmd:)}}first-stage model: {opt linear} (default), {opt probit}, {opt fprobit}, {opt poisson}, or {opt biprobit} (recursive bivariate probit for one 0/1 endogenous regressor: joint MLE instead of a control function; v0.19).  Accepts one type (applied to all endogs) or one per endog.{p_end}
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
{synopt:{opt dist:ribution(type params)}}non-Gaussian prior on gamma instead of {opt mu()}/{opt omega()}: {opt normal} {it:mean sd}, {opt uniform} {it:min max}, {opt chi2} {it:df}, {opt t} {it:df}, {opt gamma} {it:shape scale}, {opt poisson} {it:mean}, or {opt special} {it:varname} (empirical draws); one parameter set per instrument{p_end}
{synopt:{opt graphomega(numlist)}}delta sweep: rerun the LTZ simulation for each prior variance in the list (one value per step, broadcast to all instruments; Gaussian priors only) and store/plot bounds vs delta{p_end}
{synopt:{opt graphmu(numlist)}}prior means for the sweep steps (default: zeros){p_end}
{synopt:{opt graphdelta(numlist)}}x-axis labels for the sweep steps (default: the omega values){p_end}

{syntab:Inference}
{synopt:{opt gu:nits(latent|pr|structural)}}units of the user-supplied gamma values: {opt latent} (default, probit-index units; synonym {opt index}), {opt pr} (probability points; divided by m_phi from the baseline fit), or {opt structural} (structural-latent units, total error variance 1; divided by sigma_tilde from the baseline fit -- a conversion derived for a LINEAR first stage; not for use with a probit, fprobit, or poisson first stage, see {help apexog##gunits:gunits()}){p_end}
{synopt:{opt level(#)}}confidence level as a fraction ({cmd:level(.95)}, the default) or as a percent ({cmd:level(95)}){p_end}
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
{synopt:{opt noheader}}suppress the results display (notes are still printed; everything is stored in {cmd:e()}){p_end}
{synopt:{opt gra:ph}}produce a graph of UCI bounds vs gamma (uci) or histogram of LTZ draws (ltz){p_end}
{synopt:{opt graphopt:ions(}{it:string}{cmd:)}}extra options passed to the {cmd:twoway} call (overrides defaults){p_end}
{synopt:{opt debug}}print the per-iteration probit fits and Mata intermediates{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:apexog} answers the question:
{it:in a control-function probit model, how large would the}
{it:exclusion-restriction violation gamma have to be before the estimated}
{it:Average Partial Effect of the endogenous regressor crosses zero or loses}
{it:significance?}  It adapts Conley, Hansen, and Rossi's
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
directly (regress/probit/fracreg/poisson/biprobit plus Mata) and never calls
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
precision (checked in the package's validation suite).  With
{opt firststage(biprobit)} (v0.19) the main equation is instead the recursive
bivariate probit fit jointly with the 0/1 treatment, and the SE is the MLE
delta method; see the Options and Methods sections for when to prefer it.


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

{phang2}{opt biprobit} {hline 2} (v0.19) Recursive bivariate probit instead of
a control function.  Endogenous d must be binary {c -(}0, 1{c )-} and the only
endogenous regressor.  At every gamma {cmd:apexog} fits
{cmd:biprobit (depvar = indepvars d, offset(z'gamma)) (d = iv indepvars)}
by maximum likelihood and reports the average-structural-function discrete
change mean(Phi(eta|d=1) - Phi(eta|d=0)) from the outcome equation, with a
delta-method SE from the joint-MLE {cmd:e(V)} (robust or cluster, as
requested).  There is no generated regressor and hence no Murphy-Topel step.
This is the exact estimator when the structural errors are jointly normal,
where the generalized-residual control function of {opt probit} is only an
approximation; see {help apexog##gr:Methods and formulas}.  Not available
with an AIE target.  {opt gunits(structural)} is accepted and identical to the
default for this engine (sigma_tilde = 1).{p_end}

{pstd}
{ul:Which one for a binary treatment?}  {opt probit} reproduces
{help cfprobit} exactly and relies on Wooldridge's sufficiency assumption
D(eps | z, d) = D(eps | GR); {opt biprobit} relies on joint normality of the
outcome and treatment errors.  Neither set of assumptions nests the other
(Wooldridge 2014).  Under joint normality the two-step APE is biased, and
increasingly so in the error correlation (the large-N table under Methods),
so when that is the maintained model prefer {opt biprobit}.  The results
header prints a note whenever a generalized-residual first stage is in use.

{phang}
{ul:Phase-1 (v0.10) shortcut:} the type can also be declared in the
endogenous-variable clause itself via the standard Stata factor-variable
prefix:

{phang2}{cmd:(c.d = z)} -- treats d as continuous -> default {opt firststage(linear)}{p_end}
{phang2}{cmd:(i.d = z)} -- treats d as binary -> default {opt firststage(probit)}{p_end}

The {opt firststage()} option still overrides the prefix.  The prefix is
stripped before passing d to the main probit (binary-ness is then
auto-detected from the values, so the discrete-change APE is used).  The
{cmd:i.} prefix keeps its probit (control-function) default; the recursive
bivariate probit must be requested explicitly with {opt firststage(biprobit)}.

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
instruments (the CHR curse of dimensionality).  If {cmd:grid^k} exceeds
{cmd:c(max_matsize)}, {cmd:uci} stops with an actionable error: coarsen
{opt grid()}, use fewer instruments, or switch to {cmd:ltz} (which draws
from the prior and needs no grid).

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
{opt distribution(type params)} {hline 2} non-Gaussian prior on gamma for
{cmd:ltz}, mutually exclusive with {opt mu()}/{opt omega()}.  {it:type} is one of
{cmd:normal}, {cmd:uniform}, {cmd:chi2}, {cmd:t}, {cmd:gamma}, {cmd:poisson}, or
{cmd:special}; {it:params} lists the parameters for EACH instrument in turn:
two per instrument for {cmd:normal} (mean, sd), {cmd:uniform} (min, max), and
{cmd:gamma} (shape, scale); one per instrument for {cmd:chi2} (df), {cmd:t}
(df), and {cmd:poisson} (mean); and, for {cmd:special}, the name of a variable
per instrument whose observed values are sampled with replacement as the
empirical prior.  The conventions mirror {help plausexog}.  Examples with one
instrument: {cmd:distribution(uniform -0.1 0.1)}, {cmd:distribution(t 5)},
{cmd:distribution(special mygammas)}.  {opt gunits()} rescales the draws
themselves, so every type works with {opt gunits(pr)} and
{opt gunits(structural)}.  The delta sweep ({opt graphomega()}) requires the
Gaussian prior.

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
relabels the x-axis (default: the omega values).  {opt gunits(pr)} and
{opt gunits(structural)} apply to sweep values too.  Each step costs a full
set of {opt iterations()} draws.

{dlgtab:Inference}

{marker gunits}{...}
{phang}
{opt gunits(latent|pr|structural)} {hline 2} units in which the user-supplied
gamma quantities ({opt gmin()}/{opt gmax()}, {opt mu()}/{opt omega()},
{opt distribution()} draws, and sweep values) are expressed.  Internally the
violation always enters as an offset on the PROBIT-INDEX scale; this option
tells apexog which of the three elicitation scales your numbers live on, and
it converts them after the baseline fit.  Getting the scale right matters:
Simulations show that feeding a structural-units prior to
the index unrescaled is a scale miscalibration whose wide-prior LTZ coverage
is 2*Phi(z*sigma_tilde)-1 (about 0.88 instead of 0.95 at an endogeneity
correlation of 0.6).

{phang2}{cmd:gunits(latent)} (default; synonym {cmd:index}) {hline 2} gamma
is a coefficient in the PROBIT-INDEX latent equation: the scale on which the
probit normalizes the error sd CONDITIONAL on the control functions to 1.
This is the scale the offset uses internally, so no conversion is applied.
Note this is NOT the structural-latent scale (see {cmd:structural} below);
the two differ by the factor sigma_tilde.{p_end}

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
same elicited violation.  Probability-scale elicitation is also IMMUNE to
the index-vs-structural units gap described under {cmd:structural}.{p_end}

{phang2}{cmd:gunits(structural)} {hline 2} gamma values are on the
STRUCTURAL latent scale: the researcher reasons about
y* = xb + gamma*z + eps with the TOTAL error variance normalized to 1 (the
scale on which simulation DGPs and textbook latent-variable models are
usually written).  The probit index instead normalizes the error sd
CONDITIONAL on the control functions, so index units = structural units /
sigma_tilde with sigma_tilde = sd(eps | CF) < 1 whenever there is real
endogeneity.  apexog estimates sigma_tilde from the baseline fit's
control-function coefficient(s) as sigma_tilde =
(1 + b_cf' Var(chat) b_cf)^(-1/2) -- for a single control function with a
unit-variance first-stage residual this is (1 + theta_cf^2)^(-1/2) -- and
divides every gamma quantity by it ({opt omega()} by sigma_tilde^2;
non-Gaussian draws are rescaled draw-by-draw).  cfinteract
control-function product terms are ignored in sigma_tilde.{p_end}

{phang2}{ul:Linear first stage only.}  The sigma_tilde conversion (and the
always-on header line "1 structural-latent unit of gamma = 1/sigma_tilde
latent-index units") is derived for, and exact with, a LINEAR first stage,
where the control function estimates the structural first-stage error.  With
{opt firststage(probit)}, {opt fprobit}, or {opt poisson} it does NOT apply:
in Monte Carlo simulations
(n = 5000, R = 1000 per cell, gamma_true
entered in structural units) the structural gamma enters the probit index
essentially UNSCALED -- refitting at the raw structural gamma gives Wald
coverage of the true APE 0.936 against 0.589 for the refit at
gamma/sqrt(1 - rho^2) (binary d, probit first stage, gamma_true = 0.2) and
0.950 against 0.881 (fractional d, fprobit first stage), and the
interpolated zero-bias offset lies within 5% of the structural gamma in all
three generalized-residual designs (binary, fractional, and count d).  So do
not use {opt gunits(structural)} with those first stages: elicit on the index
scale (the default) or on the probability scale ({opt gunits(pr)}), which is
immune to the two latent scales.  {cmd:apexog} still accepts the combination
for backward compatibility but prints a warning line in the results header;
nothing else changes.  {cmd:e(stilde)} is then the linear-first-stage formula
evaluated on a generalized residual (0.92, 0.92, 0.98 at gamma_true = 0 in
those designs), NOT a conversion factor, and should be read as such.  For a
0/1 treatment {opt firststage(biprobit)} has sigma_tilde = 1 by construction,
consistent with this evidence.{p_end}

{phang2}For all three settings the conversion is first-order (exact as
gamma -> 0; m_phi and sigma_tilde are held at their baseline values).  All
gamma values displayed and stored ({cmd:e(grid)}, graphs) are in
PROBIT-INDEX (latent) units; a header note records any conversion, and both
conversion factors are always reported in the header and returned in
{cmd:e(mphi)} and {cmd:e(stilde)}.  Both are computed from the gamma = 0
baseline fit, which is misspecified whenever the true violation is nonzero,
so they describe that fit rather than being structural constants: in the
Monte Carlo simulations the mean {cmd:e(stilde)} moves with the true violation
(linear first stage 0.80, 0.91, 1.00 at gamma_true = 0, 0.2, 0.5 against a
true sqrt(1 - rho^2) of 0.80; Poisson first stage 0.98, 0.93, 0.61): the
omitted z'gamma term is partly absorbed by the control-function coefficient
through the first stage.{p_end}

{phang2}{ul:Exact probability scale of the refits (v0.20).}  Whatever the
input scale, {cmd:apexog} also reports, at every refit, the EXACT
probability-scale value of each instrument's violation, computed from the
model refitted at that violation over the estimation sample, holding the
treatment, the control functions and the other instruments' offsets at their
observed values: for a 0/1 instrument the average change in Pr(y=1) when it
switches from 0 to 1, for any other instrument the average derivative
mean(phi(eta))*gamma.  The first-order value m_phi*gamma used by
{opt gunits(pr)} approximates it; the two agree for small violations and can
differ by several percent for a large violation of a 0/1 instrument on an
outcome near 0 or 1 (about 10% in an application with a 0/1 instrument and an
outcome near 1).  With a linear first stage, as many instruments as
endogenous regressors and a continuous instrument they coincide exactly,
because the offset then lies in the span of the regressors and control
functions and the refit leaves the fitted index unchanged.  The UCI output
prints the exact range beside each gamma range and stores the values in
{cmd:e(grid_pr)}; with {opt gunits(pr)} and {opt graph} the x-axis is the
exact scale.  The LTZ output prints the exact mean and sd of the draws and
stores {cmd:e(gamma_sim)} (latent draws) and {cmd:e(gamma_pr_sim)}.  The
input conversion itself is unchanged (first order).{p_end}

{phang}
{opt level(#)} {hline 2} confidence level (both methods; for {cmd:ltz} it is
the percentile level of the convolved draws), given as a fraction in (0, 1)
or, as in other Stata commands, as a percent in [1, 100): {cmd:level(.9)} and
{cmd:level(90)} are equivalent.  {cmd:e(level)} always stores the fraction.
Default 0.95.

{phang}
{opt vce(vcetype)} {hline 2} variance-covariance estimator:

{phang2}{cmd:vce(robust)} {hline 2} heteroskedasticity-consistent (default).
The Murphy-Topel correction is applied to account for the generated
regressor, and the effect's direct dependence on the first-stage estimate is
included (see {help apexog##methods:Methods and formulas}).{p_end}

{phang2}{cmd:vce(cluster} {it:clustvar}{cmd:)} {hline 2} cluster-robust on
{it:clustvar}, which may be numeric or string; the usual abbreviations such as
{cmd:vce(cl id)} are accepted.  The cluster correction is propagated through
both the first-stage VCE and the Murphy-Topel meat (cluster-summed scores).
The number of clusters is stored in {cmd:e(N_clust)} and the variable name in
{cmd:e(clustvar)}.{p_end}

{phang}
For bootstrap inference, wrap {cmd:apexog} in Stata's
{help bootstrap:bootstrap} prefix:

{phang2}
{cmd:. bootstrap, reps(500): apexog uci promoted tenure manager (training = subsidy), gmin(-0.3) gmax(0.3)}
{p_end}


{marker examples}{...}
{title:Examples}

{pstd}
The examples below run on this small simulated dataset (a binary promotion
outcome, two continuous and one binary endogenous treatment, two instruments,
and a few exogenous controls/moderators).  Paste it in to follow along:

{phang}{cmd:. clear}{p_end}
{phang}{cmd:. set seed 8675309}{p_end}
{phang}{cmd:. set obs 4000}{p_end}
{phang}{cmd:. gen double tenure     = rnormal()} {it:// exogenous control (cts)}{p_end}
{phang}{cmd:. gen double experience = rnormal()} {it:// exogenous moderator (cts)}{p_end}
{phang}{cmd:. gen byte   manager    = runiform() < 0.5} {it:// exogenous moderator (binary)}{p_end}
{phang}{cmd:. gen double subsidy    = rnormal()} {it:// instrument}{p_end}
{phang}{cmd:. gen double grant      = rnormal()} {it:// second instrument}{p_end}
{phang}{cmd:. gen double a1         = rnormal()} {it:// UNOBSERVED (drives training/certified)}{p_end}
{phang}{cmd:. gen double a2         = rnormal()} {it:// UNOBSERVED (drives coaching)}{p_end}
{phang}{cmd:. gen double e          = 0.4*a1 + 0.4*a2 + sqrt(0.68)*rnormal()}{p_end}
{phang}{cmd:. gen double training   = 0.7*subsidy + 0.3*tenure + a1} {it:// endogenous (cts)}{p_end}
{phang}{cmd:. gen double coaching   = 0.6*grant + 0.3*experience + a2} {it:// endogenous (cts)}{p_end}
{phang}{cmd:. gen byte   certified  = 0.9*subsidy + 0.3*tenure + a1 > 0.4} {it:// endogenous (binary)}{p_end}
{phang}{cmd:. gen double ystar = -0.3 + 0.5*training + 0.3*coaching + 0.3*tenure ///}{p_end}
{phang}{cmd:        + 0.2*experience + 0.3*manager + 0.4*training*manager ///}{p_end}
{phang}{cmd:        + 0.25*tenure*experience + e}{p_end}
{phang}{cmd:. gen byte   promoted = ystar > 0} {it:// binary outcome}{p_end}
{phang}{cmd:. gen int    team = ceil(_n/20)} {it:// 200 clusters of 20}{p_end}

{pstd}{ul:UCI for the APE of a continuous treatment, gamma in [-0.3, 0.3]:}{p_end}
{phang}{cmd:. apexog uci promoted tenure manager (training = subsidy), gmin(-0.3) gmax(0.3)}{p_end}

{pstd}{ul:UCI with a binary endogenous treatment and probit first stage:}{p_end}
{phang}{cmd:. apexog uci promoted tenure manager (certified = subsidy),}{p_end}
{phang}{cmd:    firststage(probit) gmin(-0.1) gmax(0.1) grid(11)}{p_end}

{pstd}{ul:Same treatment with the recursive bivariate probit (exact under joint}
{ul:normality, where the probit first stage above is only an approximation):}{p_end}
{phang}{cmd:. apexog uci promoted tenure manager (certified = subsidy),}{p_end}
{phang}{cmd:    firststage(biprobit) gmin(-0.1) gmax(0.1) grid(11)}{p_end}
{phang}{cmd:. apexog ltz promoted tenure manager (certified = subsidy),}{p_end}
{phang}{cmd:    firststage(biprobit) mu(0) omega(0.0025) iterations(1000) seed(12345)}{p_end}

{pstd}{ul:LTZ with a tight Gaussian prior centered at 0:}{p_end}
{phang}{cmd:. apexog ltz promoted tenure manager (training = subsidy),}{p_end}
{phang}{cmd:    mu(0) omega(0.0025) iterations(2000) seed(12345)}{p_end}

{pstd}{ul:Same prior elicited on the other two gamma scales (see gunits()):}{p_end}
{phang}{cmd:. apexog ltz promoted tenure manager (training = subsidy),}{p_end}
{phang}{cmd:    mu(0) omega(0.0025) iterations(2000) seed(12345) gunits(structural)} {it:// structural-latent units, /sigma_tilde}{p_end}
{phang}{cmd:. apexog ltz promoted tenure manager (training = subsidy),}{p_end}
{phang}{cmd:    mu(0) omega(1.5e-06) iterations(2000) seed(12345) gunits(pr)} {it:// probability points, /m_phi}{p_end}

{pstd}{ul:Cluster-robust SEs:}{p_end}
{phang}{cmd:. apexog uci promoted tenure manager (training = subsidy), gmin(-0.5) gmax(0.5)}{p_end}
{phang}{cmd:    vce(cluster team)}{p_end}

{pstd}{ul:Visualize the UCI:}{p_end}
{phang}{cmd:. apexog uci promoted tenure manager (training = subsidy), gmin(-0.5) gmax(0.5) grid(21) graph}{p_end}

{pstd}{ul:Phase-1 syntax with i. prefix and target() narrowing (v0.10):}{p_end}
{phang}{cmd:. apexog uci promoted tenure manager (i.certified = subsidy), gmin(-0.1) gmax(0.1)} {it:// = firststage(probit)}{p_end}
{phang}{cmd:. apexog uci promoted tenure (training coaching = subsidy grant), gmin(0 0) gmax(0 0) target(training)} {it:// only training's APE}{p_end}

{pstd}{ul:Phase-2 syntax: target(c.X#i.d) for AIE with prefix-driven type (v0.11):}{p_end}
{phang}{cmd:. apexog uci promoted tenure experience (training = subsidy), gmin(-0.3) gmax(0.3) target(c.tenure#c.experience)} {it:// cts-by-cts AIE}{p_end}
{phang}{cmd:. apexog uci promoted tenure (i.certified = subsidy), gmin(-0.3) gmax(0.3) target(c.tenure#i.certified)} {it:// cts-by-binary (firstdiff)}{p_end}
{phang}{cmd:. apexog uci promoted manager (i.certified = subsidy), gmin(-0.3) gmax(0.3) target(i.manager#i.certified)} {it:// binary-by-binary}{p_end}

{pstd}{ul:Sensitivity of the Average Interaction Effect for two exogenous}
{ul:moderators under plausible exclusion-restriction violations:}{p_end}
{phang}{cmd:. apexog uci promoted tenure experience (training = subsidy), gmin(-0.3) gmax(0.3) grid(11) aie(tenure experience)}{p_end}
{phang}{cmd:. apexog ltz promoted tenure experience (training = subsidy), mu(0) omega(0.0025) iterations(2000) aie(tenure experience)}{p_end}

{pstd}{ul:Canonical case (v0.13): an ENDOGENOUS treatment interacted with a}
{ul:moderator.  training*manager is itself endogenous, so declare it as a}
{ul:SECOND endogenous regressor and instrument it (here with subsidy and}
{ul:subsidy*manager), exactly as in cfprobit / plausexog.  cfinteract is}
{ul:auto-on; the AIE is reported on the probability scale:}{p_end}
{phang}{cmd:. apexog uci promoted tenure manager (c.training c.training#c.manager = subsidy c.subsidy#i.manager),}{p_end}
{phang}{cmd:      target(c.training#i.manager) gmin(-0.3 -0.3) gmax(0.3 0.3) grid(11)}{p_end}

{pstd}{it:Note:} the older {cmd:aie(d X)} form (with d in a single endog clause
{cmd:(d = Zinst)}) treats {cmd:d*X} as an exogenous control -- a mis-specification
when d is endogenous.  Prefer the explicit two-endog form above.{p_end}

{pstd}{it:Moderator type.}  For an endogenous-treatment interaction, the moderator
type is set by the prefix in {opt target()}: {cmd:target(c.d#c.W)} gives the
cts-by-cts cross-partial, {cmd:target(c.d#i.W)} the cts-by-binary first
difference (validated to machine precision against the manual decomposition for a
binary exogenous moderator).  A binary {it:endogenous treatment} (probit first
stage with a generalized-residual control function for the interaction) is not yet
supported on this path.{p_end}

{pstd}{it:cfinteract scope.}  {cmd:cfinteract} adds the {bf:pairwise} products of
the control functions.  For two endogenous regressors (the canonical
endogenous-interaction AIE) this reproduces {help cfprobit}'s {cmd:cfinteract}
exactly.  With three or more control functions, {cmd:cfprobit}'s {cmd:##}
expansion also includes the three-way and higher products, which apexog does not
add; for exact {cmd:cfprobit} equivalence use at most two endogenous
regressors.{p_end}


{marker limits}{...}
{title:Input requirements and limitations}

{pstd}
{cmd:apexog} checks its inputs up front and stops with an explicit message
(return code 198) rather than failing deep inside a fit or, worse, returning a
silently wrong answer:

{phang2}- {it:depvar} must be coded 0/1 on the estimation sample (Stata's {cmd:probit} would treat any nonzero value as 1, but the two-step variance needs the 0/1 coding); recode first, e.g. {cmd:gen byte y01 = (y != 0)}.{p_end}
{phang2}- {it:depvar}, the controls, the endogenous regressors, and the instruments must be disjoint.  An instrument that merely duplicates a control (e.g. the main effect of a {cmd:##} term) is dropped by the full-rank cleaning with a note.{p_end}
{phang2}- {opt gmin()} may not exceed {opt gmax()}; {opt grid()} must be at least 2; {opt omega()} must be positive; {opt iterations()} must be at least 10.{p_end}
{phang2}- {opt firststage(probit)} needs a 0/1 endogenous regressor, {opt firststage(fprobit)} one in [0,1], and {opt firststage(poisson)} a nonnegative one.{p_end}
{phang2}- {opt firststage(biprobit)} needs a 0/1 endogenous regressor that is the ONLY endogenous regressor and no AIE target; several endogenous regressors, {opt aie()}/{opt target(c.X#i.d)}, or {opt estimator(biprobit)} with it stop with an explicit message.{p_end}
{phang2}- {ul:Nonlinear first stages are approximations.}  With {opt firststage(probit)}, {opt fprobit}, or {opt poisson} the control function is Wooldridge's (2014, 2015) generalized residual.  It is consistent for the APE under his sufficiency assumption D(eps | z, d) = D(eps | GR), which does NOT hold under joint normality of the structural errors; there the two-step APE is biased, increasingly so in the endogeneity correlation (-2.7%, -11.6%, -50.6% at error correlations 0.3, 0.5, 0.8 in the large-N check under {help apexog##gr:Methods and formulas}).  {cmd:apexog} reproduces {cmd:cfprobit} exactly, so this is a property of the estimator, not of the implementation.  For a 0/1 treatment under joint normality use {opt firststage(biprobit)}; the results header carries the same note whenever one of these first stages is used.{p_end}
{phang2}- {ul:gunits(structural) is a linear-first-stage conversion.}  The sigma_tilde rescaling behind {opt gunits(structural)} (and the "1 structural-latent unit of gamma = 1/sigma_tilde latent-index units" header line) is derived for a linear first stage.  With {opt firststage(probit)}, {opt fprobit}, or {opt poisson} Monte Carlo simulations put the correct offset at the UNSCALED structural gamma (Wald coverage at the true violation 0.936 vs 0.589 for binary d and 0.950 vs 0.881 for fractional d, raw gamma vs gamma/sqrt(1 - rho^2); zero-bias offset within 5% of the structural gamma), so the option should not be used with those first stages.  {cmd:apexog} accepts the combination but prints a warning line in the results header (results unchanged), and {cmd:e(stilde)} is then the linear formula evaluated on a generalized residual, not a conversion factor.  Elicit on the index scale or with {opt gunits(pr)} instead; see {help apexog##gunits:gunits()}.{p_end}
{phang2}- A first stage that drops observations (perfect prediction in a probit first stage) or leaves the control function undefined stops the command.  Perfect prediction or a collinear control function in the MAIN probit returns missing bounds with a warning instead.{p_end}
{phang2}- Weights, the {cmd:svy} prefix, and {cmd:vce(bootstrap)} are not supported; wrap {cmd:apexog} in the {cmd:bootstrap} prefix for bootstrap inference.{p_end}
{phang2}- Variable names are unabbreviated before use, so abbreviations work; names longer than 27 characters are truncated in {cmd:e(lb_}{it:endog}{cmd:)}/{cmd:e(ub_}{it:endog}{cmd:)} and in matrix column labels (Stata's 32-character name limit).{p_end}
{phang2}- Options that belong to the other method ({opt gmin()}/{opt gmax()} with {cmd:ltz}; {opt mu()}/{opt omega()}/{opt distribution()}/{opt seed()} with {cmd:uci}) are ignored with a note.{p_end}


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
{synopt:{cmd:e(baseline_se)}}MT-corrected SE at gamma = 0 (MLE delta-method SE with {opt firststage(biprobit)}){p_end}
{synopt:{cmd:e(mphi)}}mean(normalden(eta_hat)) at the baseline fit: average dPr/dindex slope; the latent-to-probability conversion factor used by {opt gunits(pr)}{p_end}
{synopt:{cmd:e(stilde)}}sigma_tilde = (1 + b_cf' Var(chat) b_cf)^(-1/2) at the gamma = 0 baseline fit.  With a linear first stage this is the probit-index-to-structural-latent scale ratio, the conversion factor used by {opt gunits(structural)}; with a probit, fprobit, or poisson first stage it is the same formula evaluated on a generalized residual and NOT a conversion factor (see {help apexog##gunits:gunits()}); = 1 with {opt firststage(biprobit)}.  A baseline-fit quantity: it moves with the true violation{p_end}
{synopt:{cmd:e(rho)}}(with {opt firststage(biprobit)}) estimated correlation of the outcome and treatment errors at the baseline fit{p_end}
{synopt:{cmd:e(lb_}{it:endog}{cmd:)}}lower bound of the APE interval (per endog; a name longer than 27 characters is truncated){p_end}
{synopt:{cmd:e(ub_}{it:endog}{cmd:)}}upper bound of the APE interval (per endog){p_end}
{synopt:{cmd:e(lb_aie)} / {cmd:e(ub_aie)}}lower/upper bound of the AIE interval (with {opt aie()}){p_end}
{synopt:{cmd:e(N_clust)}}number of clusters (with {opt vce(cluster)}){p_end}

{p2col 5 20 24 2: Macros (additional in {opt aie()} mode)}{p_end}
{synopt:{cmd:e(target)}}{cmd:ape} (default) or {cmd:aie}{p_end}
{synopt:{cmd:e(aie_vars)}}the two variables in {opt aie(X Z)}, space-separated{p_end}
{synopt:{cmd:e(aie_spec)}}the {opt target()} interaction spec, e.g. {cmd:c.X#i.d} (AIE via {opt target()}){p_end}
{synopt:{cmd:e(target_list)}}the endogenous regressors whose APE is reported (APE mode){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:apexog}{p_end}
{synopt:{cmd:e(cmdline)}}the full command as typed{p_end}
{synopt:{cmd:e(method)}}{cmd:uci} or {cmd:ltz}{p_end}
{synopt:{cmd:e(gunits)}}{cmd:latent}, {cmd:pr}, or {cmd:structural} (units of the user-supplied gamma values){p_end}
{synopt:{cmd:e(estimator)}}{cmd:cfprobit}, or {cmd:biprobit} with {opt firststage(biprobit)}{p_end}
{synopt:{cmd:e(firststage)}}{cmd:linear}, {cmd:probit}, {cmd:fprobit}, {cmd:poisson}, or {cmd:biprobit}{p_end}
{synopt:{cmd:e(vce)}}{cmd:Murphy-Topel}, or {cmd:MLE} with {opt firststage(biprobit)}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(endog)}}names of the endogenous regressors{p_end}
{synopt:{cmd:e(exog)}}names of exogenous regressors{p_end}
{synopt:{cmd:e(instruments)}}names of instruments (after full-rank cleaning){p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (with {opt vce(cluster)}){p_end}

{p2col 5 20 24 2: Matrices (both methods; J = number of targets)}{p_end}
{synopt:{cmd:e(b)} / {cmd:e(V)}}baseline (gamma = 0) target estimates and their Murphy-Topel (biprobit: MLE delta-method) covariance{p_end}
{synopt:{cmd:e(baseline_ape_vec)} / {cmd:e(baseline_se_vec)}}{cmd:J x 1} baseline estimates and SEs (rows labelled by target){p_end}

{pstd}
For {cmd:method = uci}, also:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(baseline_lo)} / {cmd:e(baseline_hi)}}baseline (gamma=0) CI bounds{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(grid)}}{cmd:(grid^k) x (k+3J)} matrix with columns (g1, ..., gk) followed by (ape, lo, hi) for each target{p_end}
{synopt:{cmd:e(grid_pr)}}{cmd:(grid^k) x k} exact probability-scale value of each instrument's violation at each row of {cmd:e(grid)} (columns gpr1, ..., gprk; see {help apexog##gunits:gunits()}){p_end}
{synopt:{cmd:e(baseline_lo_vec)} / {cmd:e(baseline_hi_vec)}}{cmd:J x 1} baseline CI bounds{p_end}

{pstd}
For {cmd:method = ltz}, also:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(iterations)}}number of simulation draws{p_end}

{p2col 5 20 24 2: Matrices (J = number of targets)}{p_end}
{synopt:{cmd:e(ape_sim)}}{cmd:iterations x J} per-draw point APEs/AIEs (prior only){p_end}
{synopt:{cmd:e(se_sim)}}{cmd:iterations x J} per-draw Murphy-Topel SEs{p_end}
{synopt:{cmd:e(ape_sim_cnv)}}{cmd:iterations x J} convolved draws APE_i + SE_i*eps_i (basis of the headline CI){p_end}
{synopt:{cmd:e(gamma_sim)}}{cmd:iterations x k} the gamma draws, one column per instrument (latent-index units, after any {opt gunits()} conversion){p_end}
{synopt:{cmd:e(gamma_pr_sim)}}{cmd:iterations x k} exact probability-scale value of each draw at its refit (see {help apexog##gunits:gunits()}){p_end}
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
{ul:The variance of the effect (version 0.21.0).}  The APE or AIE depends on
the first-stage estimate twice: through the probit coefficients (the
Murphy-Topel channel above) and directly, because the control function
chat_i(pi_hat) sits in the index at which the effect is evaluated.  Its
influence function therefore has two terms,

{p 8 8 2}l_i = J_beta' psi_beta_i + J_pi' psi_pi_i,   psi_beta_i = H_beta^(-1) s_total_i{p_end}

{p 8 8 2}J_pi = mean_i [ dq_i/d(eta_i) * d(eta_i)/d(chat_i) * d(chat_i)/d(pi) ]{p_end}

{pstd}
where J_beta is the delta-method gradient of the effect in the probit
coefficients and J_pi its direct gradient in the first-stage parameters
(d(eta)/d(chat) = the control-function coefficient, plus the cfinteract
product terms).  SE = sqrt(sum_i l_i^2), with l_i summed within clusters
under {opt vce(cluster)}; {cmd:e(V)} is the same outer product for all
targets.  Versions through 0.20.1 used only the first term, SE =
sqrt(J_beta' V_beta_MT J_beta), as does {cmd:cfprobit} followed by
{cmd:margins} (which treats the stored control function as a fixed
regressor); so at gamma = 0 apexog's SE now differs from that benchmark by
the second term, about 1% of the SE in the validation designs (up to 3% under
strong endogeneity), in either direction.  Point estimates are unchanged.
The two-term formula was checked against a formula-free oracle: with one
instrument and a linear first stage the fitted probit index does not depend
on pi_hat, so the APE equals mean(phi) * (b_d + b_z/pi_z) from a plain
probit of y on the controls, d and z, whose stacked delta-method variance
the two-term formula reproduces to machine precision.
{opt firststage(biprobit)} is a joint MLE and is not affected.

{pstd}
For cluster-robust SEs, the meat is replaced by the sum of OUTER products
of cluster-summed scores:

{p 8 8 2}meat_cluster = sum_c (sum_{i in c} s_total_i) (sum_{i in c} s_total_i)'{p_end}

{marker gr}{...}
{pstd}
{ul:Nonlinear first stages: the generalized-residual control function is an approximation.}
For {opt firststage(probit)}, {opt fprobit}, and {opt poisson} the control
function is Wooldridge's (2014, 2015) generalized residual GR (for a 0/1
treatment GR = (d - Phi(Wpi)) phi(Wpi) / [Phi(Wpi) (1 - Phi(Wpi))]).  Adding GR
to the probit for {it:depvar} identifies the APE under Wooldridge's
sufficiency assumption

{p 8 8 2}D(eps | z, d) = D(eps | GR){p_end}

{pstd}
(Wooldridge 2015, eq. 53 and the text around it on p. 442; Wooldridge 2014,
the discussion following his eq. 66).  That assumption is NOT implied by joint
normality of the structural errors (eps, v): conditioning on the binary d
rather than on v leaves Var(eps | d, z) heteroskedastic and eps | d, z
non-normal, so the probit in (x, d, GR) is misspecified there.  Wooldridge
(2014) is explicit that with a binary endogenous regressor the control-function
approach uses "different parametric assumptions than those used in the
bivariate probit model -- and neither set subsumes the other", and that the
model should be viewed "only as an approximation"; Wooldridge (2015, p. 442)
adds that under the bivariate-probit assumptions "there is no known CF method
that consistently estimates the parameters"; [R] cfprobit (Example 2) says that
with a probit first stage one has "moved away from a structural model" and
defers to Wooldridge (2015) for when the APEs are "valid, or approximately
valid".  Because {cmd:apexog} reproduces {cmd:cfprobit} exactly, the property
is inherited.  Its size is not academic.  In the bivariate-normal design of the
validation suite (y* = -0.2 + 0.6 d + 0.5 x + z'gamma + eps, d = 1[0.8 z + 0.4 x
+ v > 0], x, z, v standard normal, corr(eps, v) = rho, gamma = 0, N = 400,000;
true discrete-change
APE 0.21065):

{p 8 8 2}rho = 0.3:   GR-CF two-step APE 0.20488 (bias -2.7%)     recursive bivariate probit 0.21071 (+0.0%){p_end}
{p 8 8 2}rho = 0.5:   GR-CF two-step APE 0.18632 (bias -11.6%)    recursive bivariate probit 0.20982 (-0.4%){p_end}
{p 8 8 2}rho = 0.8:   GR-CF two-step APE 0.10411 (bias -50.6%)    recursive bivariate probit 0.21083 (+0.1%){p_end}

{pstd}
The Murphy-Topel SE (about 0.003 here) is correctly sized for the two-step
estimator, so each bias is many standard errors; the validation suite
reproduces the rho = 0.8 row with both engines.  The results header therefore prints a note
whenever a generalized-residual first stage is used.  For a 0/1 treatment whose
error is plausibly jointly normal with the outcome's, use
{opt firststage(biprobit)}; for fractional and count treatments no exact analog
is shipped and the same caveat applies.  Relatedly, the sigma_tilde conversion
behind {opt gunits(structural)} is derived for a linear first stage and should
not be applied with these first stages (simulation evidence under
{ul:Gamma scale} below and {help apexog##gunits:gunits()}).

{pstd}
{ul:Recursive bivariate probit engine ({opt firststage(biprobit)}).}  With one
0/1 endogenous regressor d and (eps, v) bivariate normal with correlation rho,
the structural model with a possible exclusion-restriction violation is the
recursive bivariate probit

{p 8 8 2}y* = x'beta_1 + d beta_d + z'gamma + eps,   d* = z'pi_z + x'pi_x + v,   y = 1[y* > 0],   d = 1[d* > 0].{p_end}

{pstd}
At each gamma {cmd:apexog} fits {cmd:biprobit (y = x d, offset(z'gamma)) (d = z x)}
by maximum likelihood (robust or cluster-robust {cmd:e(V)}, as requested) and
takes

{p 8 8 2}APE = mean( Phi(x'beta_1 + beta_d + z'gamma) - Phi(x'beta_1 + z'gamma) ){p_end}

{pstd}
from the outcome equation: because eps is standard normal, Phi(.) is the
average structural function itself and nothing has to be averaged out.  The
delta-method SE uses the closed-form Jacobian d/dbeta_1 = mean((phi_1 - phi_0)
x), d/dbeta_d = mean(phi_1), d/dbeta_0 = mean(phi_1 - phi_0) (zero for the
treatment equation and athrho) against the joint-MLE {cmd:e(V)}; there is no
generated regressor and hence no Murphy-Topel step.  At gamma = 0 the point
estimate and SE equal {cmd:margins, dydx(d) predict(pmarg1)} after
{cmd:biprobit} with {cmd:i.d}.  The UCI and LTZ machinery is
unchanged: {opt gmin()}/{opt gmax()}, {opt mu()}/{opt omega()},
{opt distribution()}, {opt gunits(pr)}, the delta sweep, and the graphs work
exactly as for the control-function engine.  Fits after the baseline are
warm-started from the baseline estimates.  {cmd:e(rho)} returns the estimated
error correlation at the baseline fit.

{pstd}
{it:Gamma scale with biprobit.}  The outcome-equation error of the recursive
bivariate probit has unit variance, so the probit-index and structural-latent
scales coincide: sigma_tilde = 1, {cmd:e(stilde)} = 1, and
{opt gunits(structural)} is accepted as identical to the default.  The
probability-point scale ({opt gunits(pr)}) uses m_phi from the outcome equation
as usual.

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
it reproduces {help plausexog}'s analytic V + A*Omega*A' interval (checked
in the package's validation suite).  The prior-only percentile spread of
the point APEs {c -(}APE_i{c )-} is reported separately as a sensitivity band;
it answers "how does the bias-adjusted estimate move as gamma varies over the
prior", not "what is a confidence interval for the APE".  For multiple
targets the eps draws are independent across targets, which is exact for the
per-target marginal intervals reported.

{pstd}
{ul:Gamma scale.}  The violation enters as an offset in the PROBIT-INDEX
latent equation: Pr(y=1) = Phi(... + z'gamma), where the probit normalizes
the error sd CONDITIONAL on the control functions to 1.  Three elicitation
scales are supported, linked by two baseline-fit conversion factors that
are always reported in the header:

{phang2}(i) {ul:probit-index units} (the default, {opt gunits(latent)}):
the scale of the fitted probit coefficients.  No conversion.{p_end}

{phang2}(ii) {ul:probability points} ({opt gunits(pr)}): the instrument's
direct effect on Pr(y=1).  Linked to index units by the average link slope
m_phi = mean(normalden(eta_hat)): gamma_index = gamma_pr / m_phi (variances
by m_phi^2).  This makes apexog and {help plausexog} (run as an LPM on the
binary outcome) comparable on a common violation axis.{p_end}

{phang2}(iii) {ul:structural-latent units} ({opt gunits(structural)}): the
scale of y* = xb + gamma*z + eps with TOTAL error variance 1 -- the scale
simulation DGPs and textbook latent-variable models are usually written on.
Linked to index units by sigma_tilde = sd(eps | CF): gamma_index =
gamma_structural / sigma_tilde (variances by sigma_tilde^2), with
sigma_tilde estimated from the control-function coefficient(s) as
sigma_tilde = (1 + b_cf' Var(chat) b_cf)^(-1/2).  This link is derived for
a LINEAR first stage, whose residual estimates the structural first-stage
error; it does not carry over to the generalized-residual first stages (see
the scope paragraph below).{p_end}

{pstd}
Confusing scales (ii) and (iii) with (i) is not innocuous.  Because
sigma_tilde < 1 under endogeneity, a structural-units prior fed to the
index unrescaled is uniformly too narrow by the factor sigma_tilde; for a
wide Gaussian prior the LTZ coverage of a nominal 1-alpha interval is then
2*Phi(z*sigma_tilde)-1 rather than 1-alpha (about 0.88 instead of 0.95 when
the endogeneity correlation is 0.6).  Probability-scale elicitation
({opt gunits(pr)}) is immune to this gap because the pr -> index conversion
already carries the correct scale.  All conversions are local: m_phi and
sigma_tilde are evaluated at gamma = 0 and the mappings are exact in the
limit gamma -> 0.  Both factors are stored in {cmd:e(mphi)} and
{cmd:e(stilde)}.  The exact probability-scale value of each refit's
violation (the average discrete change for a 0/1 instrument, the average
derivative otherwise) is reported alongside and stored in {cmd:e(grid_pr)}
and {cmd:e(gamma_pr_sim)}; see {help apexog##gunits:gunits()}.

{pstd}
{ul:Scope of the sigma_tilde conversion.}  The derivation of (iii) uses a
linear first stage, whose residual estimates the structural first-stage
error v, so that sd(eps | CF) = sqrt(1 - rho^2).  Monte Carlo simulations
(n =
5000, R = 1000 per cell, gamma_true = 0, 0.2, 0.5 entered in the structural
latent index, UCI box [-0.3, 0.3] on the index scale) confirm it for the
linear first stage -- the Wald interval from the refit at the index-scale
offset gamma_true/sqrt(1 - rho^2) covers the conditional-on-CF APE 0.946 and
0.932 at gamma_true = 0.2 and 0.5, the refit at the raw structural gamma only
0.526 and 0.011 -- and reject it for the generalized-residual first stages:
with a binary d and a probit first stage (bivariate-normal design, rho = 0.5)
the refit at the raw gamma covers 0.936 against 0.589 for the refit at
gamma/sqrt(1 - rho^2) at gamma_true = 0.2; with a fractional d and an fprobit
first stage (rho = 0.4) 0.950 against 0.881; and interpolating the bias
between the two offsets puts the zero-bias offset at 0.195 and 0.518 (binary
d), 0.200 and 0.524 (fractional d), and 0.206 and 0.516 (count d, Poisson
first stage, against the structural slope APE) for gamma_true = 0.2 and 0.5,
within 5% of the structural gamma in every case.  So for those engines the
structural gamma enters the probit index essentially unscaled: neither
sqrt(1 - rho^2) nor {cmd:e(stilde)} (0.92, 0.92, 0.98 at gamma_true = 0 in
those designs) is the model-implied factor, and for binary and fractional d
dividing by it moves the offset the wrong way (the coverage losses 0.936 ->
0.589 and 0.950 -> 0.881 above).  {opt gunits(structural)}
should therefore be used with a linear first stage only; with
{opt firststage(probit)}, {opt fprobit}, or {opt poisson} {cmd:apexog} accepts
it unchanged but prints a warning line in the header.  Consistently, the
{opt firststage(biprobit)} engine has sigma_tilde = 1.  Note also that m_phi
and sigma_tilde are properties of the gamma = 0 baseline fit, which is
misspecified when the true violation is nonzero: in the same simulations the
mean {cmd:e(stilde)} is 0.80, 0.91, 1.00 with the linear first stage (true
sqrt(1 - rho^2) = 0.80) and 0.98, 0.93, 0.61 with the Poisson first stage at
gamma_true = 0, 0.2, 0.5: the omitted z'gamma term is partly absorbed by the
control-function coefficient through the first stage.

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
precision in all three cases.


{marker references}{...}
{title:References}

{phang}
Conley, T. G., C. B. Hansen, and P. E. Rossi.  2012.
Plausibly Exogenous.  {it:Review of Economics and Statistics} 94: 260-272.

{phang}
Murphy, K. M., and R. H. Topel.  1985.
Estimation and Inference in Two-Step Econometric Models.
{it:Journal of Business & Economic Statistics} 3: 370-379.

{phang}
Wooldridge, J. M.  2014.
Quasi-maximum likelihood estimation and testing for nonlinear models with
endogenous explanatory variables.  {it:Journal of Econometrics} 182: 226-234.
{it:[Section 6, discussion following eq. (66): with a binary endogenous regressor}
{it:the control-function and bivariate-probit assumptions differ and "neither set}
{it:subsumes the other"; the CF model is "only an approximation".]}

{phang}
Wooldridge, J. M.  2015.
Control function methods in applied econometrics.
{it:Journal of Human Resources} 50: 420-445.
{it:[p. 442: the sufficiency assumption D(u | z, y2) = D(u | r2) (his eq. 53) under}
{it:which the generalized-residual probit identifies the APE; "under the standard}
{it:bivariate probit assumptions, there is no known CF method that consistently}
{it:estimates the parameters".]}

{phang}
Clarke, D.  2020.
{cmd:plausexog}: Stata module for plausibly exogenous 2SLS.  Statistical
Software Components S458310, Boston College Department of Economics.
{it:[This package targets the structural coefficient; apexog adapts the same}
{it:framework to target the APE in a cfprobit model.]}

{phang}
Radean, M.  2023.
ginteff: A generalized command for computing interaction effects.
{it:Stata Journal} 23: 301-335.
{it:[apexog's aie() option targets the same AIE quantity as ginteff, but adds}
{it:Conley-Hansen-Rossi sensitivity bounds and a Murphy-Topel correction for}
{it:the generated regressor.]}


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
Manual:  {bf:[R] cfprobit}, {bf:[R] cfregress}, {bf:[R] biprobit}, {bf:[R] ivprobit},
{bf:[R] margins}

{psee}
Related: {help plausexog} (Clarke's 2SLS-targeting version)
