{smcl}
{* *! version 0.12.0  2026-05-27}{...}
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
{synopt:{opt aie(}{it:X Z}{cmd:)}}switch the sensitivity target from the APE to the Average Interaction Effect of the continuous-by-continuous interaction {cmd:c.}{it:X}{cmd:##c.}{it:Z}.  See {it:ginteff} (Radean 2023).{p_end}

{syntab:UCI (method = uci)}
{synopt:{opt gmin(numlist)}}lower bound of assumed gamma range; one value per instrument{p_end}
{synopt:{opt gmax(numlist)}}upper bound of assumed gamma range; one value per instrument{p_end}
{synopt:{opt grid(#)}}number of grid points per instrument; default {cmd:grid(11)}{p_end}

{syntab:LTZ (method = ltz)}
{synopt:{opt mu(numlist)}}Gaussian prior mean for each instrument's gamma{p_end}
{synopt:{opt omega(numlist)}}Gaussian prior variance for each instrument's gamma{p_end}
{synopt:{opt iter:ations(#)}}number of Monte-Carlo draws; default {cmd:iterations(1000)}{p_end}
{synopt:{opt seed(#)}}random-number seed for reproducibility{p_end}

{syntab:Inference}
{synopt:{opt level(#)}}confidence/credible level; default {cmd:level(.95)}{p_end}
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
With {opt aie(X Z)}, the same UCI/LTZ machinery instead targets the
Average Interaction Effect for the continuous-by-continuous interaction
{cmd:c.X##c.Z} (the AIE quantity reported by ginteff; Radean 2023).

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
{cmd:ltz} -- Local-to-Zero.  Draws gamma_i from a Gaussian prior
N({cmd:mu}, {cmd:omega}) and recomputes the APE for each draw.  Reports
the percentile interval of the simulated APE distribution.

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

{dlgtab:Inference}

{phang}
{opt level(#)} {hline 2} confidence level (uci) or credible level (ltz).
Must be in (0, 1).  Default 0.95.

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

{pstd}{ul:Canonical case: interaction of an endogenous treatment d with an
exogenous moderator X (the X*d term is endogenous; the wrapper handles
this automatically):}{p_end}
{phang}{cmd:. apexog uci y X (d = Zinst), gmin(-0.3) gmax(0.3) grid(11) aie(d X)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:apexog} stores the following in {cmd:e()} (common to both methods):

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}sample size{p_end}
{synopt:{cmd:e(level)}}confidence/credible level{p_end}
{synopt:{cmd:e(baseline_ape)}}APE (or AIE if {opt aie()}) at gamma = 0{p_end}
{synopt:{cmd:e(baseline_se)}}MT-corrected SE at gamma = 0{p_end}
{synopt:{cmd:e(lb_}{it:endog}{cmd:)}}lower bound of the APE interval (per endog){p_end}
{synopt:{cmd:e(ub_}{it:endog}{cmd:)}}upper bound of the APE interval (per endog){p_end}
{synopt:{cmd:e(lb_aie)} / {cmd:e(ub_aie)}}lower/upper bound of the AIE interval (with {opt aie()}){p_end}

{p2col 5 20 24 2: Macros (additional in {opt aie()} mode)}{p_end}
{synopt:{cmd:e(target)}}{cmd:ape} (default) or {cmd:aie}{p_end}
{synopt:{cmd:e(aie_vars)}}the two variables in {opt aie(X Z)}, space-separated{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:apexog}{p_end}
{synopt:{cmd:e(method)}}{cmd:uci} or {cmd:ltz}{p_end}
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

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(ape_sim)}}{cmd:iterations x 1} vector of simulated APE draws{p_end}


{marker methods}{...}
{title:Methods and formulas}

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
