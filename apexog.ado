*! apexog: APE/AIE bounds for cfprobit under plausibly exogenous IV (CHR 2012)
*! Version 0.21.2  2026-11-14
*! Copyright (c) 2026 Jesper N. Wulff, Aarhus University.  Licensed under
*! the MIT License -- see https://github.com/jespernwulff/stata-apexog
*! SPDX-License-Identifier: MIT
*! Targets the Average Partial Effect (APE) -- or, with target(c.X#i.d) or the
*! back-compat aie() shorthand, the Average Interaction Effect (AIE) in the
*! style of ginteff (Radean 2023) -- in a control-function probit model,
*! using the Conley/Hansen/Rossi (2012) UCI and LTZ frameworks.
*! Documentation: help apexog, https://jespernwulff.github.io/apexog/
*
* Version history (user-visible changes; help apexog has the details)
*   0.21.2  Display and documentation only: no estimate, standard error or
*           interval changes.  The header's standard-error label reads
*           "two-step" (the Murphy-Topel parameter covariance plus the direct
*           first-stage term, as computed since 0.21.0).  Header note when a
*           0/1 endogenous regressor has a linear first stage (a descriptive
*           fitted contrast).  With gunits(structural): a line that the
*           conversion uses the baseline sigma_tilde and, for one continuous
*           endogenous regressor with one instrument and a linear first stage,
*           the limit
*           |pi_z|/sigma_v below which a structural value has a unique
*           index-scale counterpart, with a warning when the box or the
*           prior draws reach it.  New e(gs_max).
*   0.21.1  firststage(fprobit): the first-stage weight (minus the derivative
*           of the generalized residual in the first-stage index, also the
*           first-stage Hessian weight) now uses the formula for a fractional
*           regressor, chat*(Wpi+chat) + d(1-d)*[phi/(Phi(1-Phi))]^2.
*           Versions up to 0.21.0 used the binary-probit identity
*           chat*(Wpi+chat), which holds only for a 0/1 regressor.  Point
*           estimates unchanged; fprobit standard errors move by about 1-3%
*           on average in the validation design.  All other first stages,
*           and firststage(fprobit) with a 0/1 regressor, are unchanged.
*   0.21.0  The standard error of every two-step target (APE and AIE, all
*           first stages, cfinteract, vce(cluster)) now adds the direct
*           dependence of the effect on the estimated first-stage parameters
*           through the control function (influence function
*           J_theta' psi_2 + J_pi' psi_1).  Point estimates unchanged; SEs
*           move by about 1% in the validation designs, in either direction.
*           firststage(biprobit) is unaffected.
*   0.20.1  APE standard errors with cfinteract now include the terms of the
*           control-function-product coefficients (point estimates, AIE
*           targets and fits without cfinteract are unchanged).
*   0.20.0  Exact probability-scale value of each instrument's violation at
*           every refit: e(grid_pr) for uci; e(gamma_sim) and e(gamma_pr_sim)
*           for ltz.  The gunits(pr) input conversion stays first order.
*   0.19.1  Header warning when gunits(structural) is combined with a probit,
*           fprobit or poisson first stage (the sigma_tilde conversion is
*           derived for a linear first stage).  No numeric change.
*   0.19.0  firststage(biprobit): recursive bivariate probit for one 0/1
*           endogenous regressor; e(rho).  Header note that the probit,
*           fprobit and poisson first stages use a generalized-residual
*           control function (Wooldridge 2014, 2015), an approximation under
*           joint normality.
*   0.18.0  Input validation (0/1 outcome, disjoint variable lists,
*           gmin <= gmax, first-stage support and sample checks); level(95),
*           vce(cl id), string cluster ids, e(N_clust), e(clustvar), full
*           e(cmdline).  Fixes: the debug option, the uci graph with an AIE
*           target, long variable names, if/in in the instrument cleaning.
*   0.17.0  gunits(structural); sigma_tilde in the header and e(stilde).
*   0.16.1  Early error when the uci grid exceeds c(max_matsize).
*   0.16.0  Posts e(b), e(V) and e(sample) (test, lincom, estimates store,
*           replay); graphomega() prior-width sweep, e(sweep); guards for
*           perfect prediction and collinear control functions.
*   0.15.0  gunits(pr): gamma in probability points; e(mphi).
*   0.14.0  The ltz interval convolves sampling noise with the gamma prior
*           (CHR 2012, eq. 7); the prior-only band is reported alongside.
*   0.13.x  Endogenous interaction AIE, (c.d c.d#c.X = z c.z#c.X) with
*           target(c.d#c.X) and cfinteract; factor-variable instruments
*           reduced to full rank.
*   0.12.0  Renamed from plausexog_ape to apexog.
*   0.11.0  target() accepts interaction specs for the AIE: c.X#c.Z,
*           c.X#i.d, i.X#i.d (the prefixes set the AIE formula).
*   0.10.0  target(); i./c. prefixes in the endogenous clause set the default
*           first stage.
*   0.9.0   aie(X Z): average interaction effect target (continuous or
*           binary interactors, endogenous or exogenous).
*   0.8.0   One instrument set per endogenous regressor, (d1 = z1) (d2 = z2);
*           APE bounds for every endogenous regressor.
*   0.7.0   distribution() priors for ltz; a firststage() list, one type per
*           endogenous regressor.
*   0.6.0   Multiple endogenous regressors.
*   0.5.0   vce(cluster); graph option; help file.
*   0.4.0   firststage(fprobit) and firststage(poisson).
*   0.3.0   firststage(probit); discrete-change APE for a 0/1 treatment.
*   0.2.0   Murphy-Topel two-step standard errors in Mata.
*
* Known limitation: first-stage options inside the parentheses, as in
* cfprobit's (d1 = z1, probit), are not supported (Stata's syntax command
* reads the inner comma as an option separator); use firststage() instead.

cap program drop apexog
cap program drop _apexogfit
cap program drop _apexogerr
cap program drop _apexogclean

program apexog, eclass
    version 18
    * Replay: `apexog' (no arguments) after estimation redisplays e(b)/e(V).
    if replay() {
        if "`e(cmd)'" != "apexog" error 301
        syntax [, level(real -1)]
        local _rlvl ""
        if `level' > 0 {
            if `level' < 1 local level = 100*`level'
            local _rlvl "level(`level')"
        }
        di as text "apexog `e(method)' (replay): baseline (gamma=0) target estimates"
        ereturn display, `_rlvl'
        exit
    }
    * Full command line (before `syntax' consumes it), for e(cmdline).
    local _cmdline_full `"`0'"'
    #delimit ;
    syntax anything(name=0 id="variable list")
    [if] [in]
    [,
        ESTimator(string)
        FIRststage(string)
        grid(integer 11)
        gmin(numlist)
        gmax(numlist)
        GUnits(string)
        level(real 0.95)
        mu(numlist min=1)
        omega(numlist min=1)
        DISTribution(string)
        seed(numlist min=1 max=1)
        iterations(integer 1000)
        vce(string)
        noheader
        debug
        graphomega(numlist)
        graphmu(numlist)
        graphdelta(numlist)
        GRAph
        GRAPHOPTions(string asis)
        AIE(string)
        TARget(string)
        CFINTeract
    ]
    ;
    #delimit cr

    *==========================================================================
    * (1) Parse the positional varlist.
    * Accepts both:
    *   apexog METHOD yvar [exog] (endog = iv) [, opts]
    *   apexog METHOD yvar [exog] (endog1 = iv1 [, cfopts]) (endog2 = iv2 [, cfopts]) ... [, opts]
    * (cfprobit-style multi-clause syntax; each clause has its own endog list,
    *  IV list, and optional cfopts.  Currently cfopts may include a
    *  first-stage type: linear | probit | fprobit | poisson.)
    *==========================================================================
    local 0_raw `"`0'"'
    local 0: subinstr local 0 "(" " ( ", all
    local 0: subinstr local 0 ")" " ) ", all
    local 0: subinstr local 0 "=" " = ", all
    local 0: subinstr local 0 "," " , ", all
    tokenize `0'

    local method `1'
    macro shift

    if "`method'" != "uci" & "`method'" != "ltz" {
        _apexogerr "Method must be either 'uci' or 'ltz'.  You typed: '`method''."
    }
    * level(): accept apexog's native fraction form level(.95) AND Stata's
    * percent form level(95).  e(level) always stores the fraction.
    if `level' >= 1 & `level' < 100 local level = `level'/100
    if `level' <= 0 | `level' >= 1 {
        _apexogerr "Confidence level must be in (0,1) (or 1-99.99 as a percent); you gave level(`level')."
    }
    if "`gunits'" == "probability" local gunits "pr"
    if "`gunits'" == "index"       local gunits "latent"
    if "`gunits'" == ""            local gunits "latent"
    if !inlist("`gunits'", "latent", "pr", "structural") {
        _apexogerr "gunits() must be latent (default; synonym index), pr, or structural.  You typed: gunits(`gunits')."
    }
    if "`graphomega'`graphmu'`graphdelta'" != "" & "`method'" != "ltz" {
        _apexogerr "graphomega()/graphmu()/graphdelta() (the delta sweep) require method ltz."
    }
    if "`method'" == "ltz" & "`gmin'`gmax'" != "" {
        di as text "  [note] gmin()/gmax() are UCI options; ignored with method ltz."
    }
    if "`method'" == "uci" & "`mu'`omega'`distribution'`seed'" != "" {
        di as text "  [note] mu()/omega()/distribution()/seed() are LTZ options; ignored with method uci."
    }

    local yvar `1'
    macro shift
    if "`yvar'" == "" | "`yvar'" == "(" {
        _apexogerr "No dependent variable found.  Use: METHOD depvar [exog] (endog = iv) [, opts]"
    }
    unab yvar : `yvar'
    confirm numeric variable `yvar'

    * Exog vars: everything from after yvar up to the first "(".
    local varlist1
    while `"`1'"' != "(" & `"`1'"' != "" {
        local varlist1 `varlist1' `1'
        macro shift
    }

    * Iterate over clauses.  Each clause = "( endogs = ivs [, cfopts] )".
    local n_clauses = 0
    local all_endog ""
    local all_iv    ""
    local cfopts_inline ""        // accumulate any inline cfopts for diagnostics
    * v0.13: metadata for an interaction term declared as an endogenous
    * regressor (endogenous treatment x moderator).  Empty when none present.
    local aieendog_tempvar ""
    local aieendog_X    ""
    local aieendog_Z    ""
    local aieendog_Xpfx ""
    local aieendog_Zpfx ""
    local aieendog_spec ""
    while `"`1'"' == "(" {
        local ++n_clauses
        macro shift               // skip (
        local endog_`n_clauses' ""
        while `"`1'"' != "=" & `"`1'"' != "" {
            local endog_`n_clauses' `endog_`n_clauses'' `1'
            macro shift
        }
        if `"`1'"' != "=" {
            _apexogerr "Bad clause: expected `=' inside parenthesized clause `n_clauses'."
        }
        macro shift               // skip =
        local iv_`n_clauses' ""
        while `"`1'"' != "," & `"`1'"' != ")" & `"`1'"' != "" {
            local iv_`n_clauses' `iv_`n_clauses'' `1'
            macro shift
        }
        if `"`1'"' == "," {
            macro shift           // skip ,
            local cfopts_`n_clauses' ""
            while `"`1'"' != ")" & `"`1'"' != "" {
                local cfopts_`n_clauses' `cfopts_`n_clauses'' `1'
                macro shift
            }
        }
        if `"`1'"' != ")" {
            _apexogerr "Bad clause: expected `)' to close parenthesized clause `n_clauses'."
        }
        macro shift               // skip )

        * Validate
        if "`endog_`n_clauses''" == "" _apexogerr "Clause `n_clauses' has no endogenous variable(s)."
        if "`iv_`n_clauses''" == ""    _apexogerr "Clause `n_clauses' has no instrument(s)."

        * Phase 1 / v0.13: parse endog tokens BEFORE fvexpand.
        *   - i./c. prefix on a plain var declares a TYPE (first-stage default).
        *   - a token containing `#' is a two-way INTERACTION declared as an
        *     endogenous regressor (the canonical endogenous-treatment x
        *     moderator case).  We materialise it as a product tempvar so it can
        *     serve as a regression LHS, and record its operands for AIE wiring.
        local _stripped_endog ""
        local _pos = 0
        foreach _tok of local endog_`n_clauses' {
            local ++_pos
            if strpos("`_tok'", "#") {
                if strpos("`_tok'", "##") {
                    _apexogerr "Use a single # for an interaction endogenous regressor (e.g. c.d#c.X) and list the main effects separately; got `_tok'."
                }
                local _opnd : subinstr local _tok "#" " ", all
                local _nop : word count `_opnd'
                if `_nop' != 2 {
                    _apexogerr "Only two-way interaction endogs are supported; got `_tok'."
                }
                if "`aieendog_tempvar'" != "" {
                    _apexogerr "Only one interaction endogenous regressor is supported; got a second: `_tok'."
                }
                local _o1 : word 1 of `_opnd'
                local _o2 : word 2 of `_opnd'
                local _b1 "`_o1'"
                local _f1 ""
                if substr("`_o1'",1,2)=="i." {
                    local _f1 "i"
                    local _b1 = substr("`_o1'",3,.)
                }
                else if substr("`_o1'",1,2)=="c." {
                    local _f1 "c"
                    local _b1 = substr("`_o1'",3,.)
                }
                local _b2 "`_o2'"
                local _f2 ""
                if substr("`_o2'",1,2)=="i." {
                    local _f2 "i"
                    local _b2 = substr("`_o2'",3,.)
                }
                else if substr("`_o2'",1,2)=="c." {
                    local _f2 "c"
                    local _b2 = substr("`_o2'",3,.)
                }
                confirm numeric variable `_b1'
                confirm numeric variable `_b2'
                tempvar _ixz
                qui gen double `_ixz' = `_b1' * `_b2'
                local aieendog_tempvar "`_ixz'"
                local aieendog_X    "`_b1'"
                local aieendog_Z    "`_b2'"
                local aieendog_Xpfx "`_f1'"
                local aieendog_Zpfx "`_f2'"
                local aieendog_spec "`_tok'"
                * a generated product is a continuous regressor: linear 1st stage
                local endog_prefix_`n_clauses'_`_pos' "c"
                local _stripped_endog "`_stripped_endog' `_ixz'"
            }
            else {
                local _pfx ""
                local _bare "`_tok'"
                if substr("`_tok'", 1, 2) == "i." {
                    local _pfx  "i"
                    local _bare = substr("`_tok'", 3, .)
                }
                else if substr("`_tok'", 1, 2) == "c." {
                    local _pfx  "c"
                    local _bare = substr("`_tok'", 3, .)
                }
                local endog_prefix_`n_clauses'_`_pos' "`_pfx'"
                local _stripped_endog "`_stripped_endog' `_bare'"
            }
        }
        local endog_`n_clauses' : list clean _stripped_endog

        * Expand factor-variable lists (now on the bare names).  fvunab first:
        * abbreviated names resolve and nonexistent names error up front.
        fvunab _endog_ua : `endog_`n_clauses''
        fvexpand `_endog_ua'
        local endog_`n_clauses' `r(varlist)'
        fvunab _iv_ua : `iv_`n_clauses''
        fvexpand `_iv_ua'
        local iv_`n_clauses' `r(varlist)'

        local all_endog `all_endog' `endog_`n_clauses''
        local all_iv    `all_iv'    `iv_`n_clauses''
    }
    if `n_clauses' < 1 {
        _apexogerr "No clauses found.  Use: METHOD yvar [exog] (endog = iv) [, opts]"
    }

    * Expand exog (fvunab first: abbreviations resolve, bad names error early)
    if "`varlist1'" != "" {
        fvunab varlist1 : `varlist1'
    }
    fvexpand `varlist1'
    local varlist1 `r(varlist)'

    * Build varlist2 (all endogs in order) and varlist_iv (UNIQUE IVs in order)
    local varlist2 : list uniq all_endog
    local varlist_iv : list uniq all_iv

    local count1   : word count `varlist1'
    local count2   : word count `varlist2'
    local count_iv : word count `varlist_iv'

    * Guard: depvar, endogenous regressors, controls, and instruments must be
    * disjoint -- an endog among the controls, or a regressor doubling as its
    * own instrument, would otherwise surface only as a cryptic collinearity
    * failure deep inside the fit.  (An instrument that merely duplicates a
    * control, e.g. a ## main effect, is handled by the full-rank cleaning
    * below with a note, as before.)
    local _yin_x : list yvar in varlist1
    local _yin_d : list yvar in varlist2
    local _yin_z : list yvar in varlist_iv
    if `_yin_x' | `_yin_d' | `_yin_z' {
        _apexogerr "The dependent variable `yvar' also appears among the regressors or instruments."
    }
    local _ovl : list varlist2 & varlist1
    if "`_ovl'" != "" {
        _apexogerr "Variable(s) listed as both endogenous and exogenous regressor: `_ovl'."
    }
    local _ovl : list varlist2 & varlist_iv
    if "`_ovl'" != "" {
        _apexogerr "Variable(s) listed as both endogenous regressor and instrument: `_ovl'."
    }

    * Build per-endog: which clause does this endog belong to, hence its IV list,
    * cfopts (first-stage type), and Phase-1 i./c. type prefix.  Endogs in
    * `varlist2' are listed in order of FIRST appearance across clauses; for
    * each, we look up the clause and the within-clause position.
    local j = 0
    foreach endogvar of local varlist2 {
        local ++j
        local fv_prefix_for_endog_`j' ""
        forvalues c = 1/`n_clauses' {
            local _in : list endogvar in endog_`c'
            if `_in' {
                local iv_for_endog_`j' "`iv_`c''"
                local cfopts_for_endog_`j' "`cfopts_`c''"
                * Find this endogvar's within-clause position to look up its prefix.
                local _pos_in_clause = 0
                local _i = 0
                foreach _v of local endog_`c' {
                    local ++_i
                    if "`_v'" == "`endogvar'" {
                        local _pos_in_clause = `_i'
                        continue, break
                    }
                }
                if `_pos_in_clause' > 0 {
                    local fv_prefix_for_endog_`j' "`endog_prefix_`c'_`_pos_in_clause''"
                }
                continue, break
            }
        }
        if "`iv_for_endog_`j''" == "" {
            _apexogerr "Internal: couldn't map endog `endogvar' to a clause."
        }
    }

    *==========================================================================
    * (1b) v0.13.1: clean instruments to FULL RANK for the Murphy-Topel W.
    *      fvexpand of factor instruments (e.g. dummy##c.eldisadv) leaves base
    *      levels (0b.) -- structural zeros -- and bare-vs-factor / moderator
    *      duplicates that are collinear with the exog controls.  Those degenerate
    *      columns make the Mata first-stage W rank-deficient (singular M_pi ->
    *      a non-PSD VCE / missing SE).  We drop them (data-based, _rmcoll against
    *      the controls) and MATERIALISE each surviving instrument as a plain
    *      tempvar so no later factor re-expansion can reintroduce a base level.
    *      `ivlabels' keeps the original term names for display/ereturn.  This is
    *      a no-op for already-clean continuous instruments.
    *==========================================================================
    *--- (1a) Parse vce() -- BEFORE the estimation sample is marked, so that a
    *    missing cluster id is excluded from the sample used by the instrument
    *    cleaning and by every fit alike.  Accept vce(robust) [default] and
    *    vce(cluster <varname>), with the usual abbreviations (vce(r),
    *    vce(cl id)).  A STRING cluster id is grouped into a numeric tempvar
    *    (the Mata views need a numeric id).
    if "`vce'" == "" local vce "robust"
    local vcefirst  : word 1 of `vce'
    local vcesecond : word 2 of `vce'
    local vcethird  : word 3 of `vce'
    local clustvar_name ""
    if "`vcethird'" != "" {
        _apexogerr "vce() must be robust or 'cluster <varname>'.  You gave: `vce'."
    }
    if "`vcefirst'" == substr("robust", 1, max(1, length("`vcefirst'"))) {
        if "`vcesecond'" != "" _apexogerr "vce(robust) takes no argument.  You gave: `vce'."
        local vcetype  = "robust"
        local clustvar = ""
    }
    else if "`vcefirst'" == substr("cluster", 1, max(2, length("`vcefirst'"))) {
        if "`vcesecond'" == "" {
            _apexogerr "vce(cluster ...) requires a cluster variable name."
        }
        unab vcesecond : `vcesecond'
        local clustvar_name "`vcesecond'"
        capture confirm numeric variable `vcesecond'
        if _rc {
            tempvar _clnum
            qui egen long `_clnum' = group(`vcesecond')
            local clustvar "`_clnum'"
        }
        else local clustvar "`vcesecond'"
        local vcetype  = "cluster"
    }
    else if "`vcefirst'" == "bootstrap" {
        _apexogerr "vce(bootstrap) is not supported directly.  Wrap apexog in" /*
            */ " Stata's `bootstrap' prefix instead:" /*
            */ "   bootstrap, reps(500): apexog uci y X (d=Z), gmin(-0.3) gmax(0.3)" /*
            */ " (MT SEs are asymptotically equivalent to bootstrap; this is mainly" /*
            */ " useful for finite-sample inference.)"
    }
    else {
        _apexogerr "vce() must be robust or 'cluster <varname>'.  You gave: `vce'."
    }

    local _raw_iv_n : word count `varlist_iv'
    * The estimation sample (if/in + missing values in every variable used,
    * incl. the cluster id) is fixed HERE, once, so the instrument cleaning
    * below and every fit afterwards see exactly the same observations.
    marksample touse
    markout `touse' `yvar' `varlist1' `varlist2' `varlist_iv' `aieendog_tempvar' `clustvar'
    qui count if `touse'
    if r(N) == 0 _apexogerr "No observations remain after if/in and missing-data filtering."
    * exog controls as plain tempvars (collinearity baseline)
    local _ctrl_tv ""
    foreach _c of local varlist1 {
        fvrevar `_c'
        local _ctrl_tv `_ctrl_tv' `r(varlist)'
    }
    * shared dedup'd pool of clean instruments (tempvar + original-name label)
    local _ivpool_tv    ""
    local _ivpool_names ""
    forvalues j = 1/`count2' {
        fvexpand `iv_for_endog_`j''
        local _allexp `r(varlist)'
        local _cand_names ""
        local _cand_tv    ""
        foreach _t of local _allexp {
            fvrevar `_t'
            tempvar _cc
            qui gen double `_cc' = `r(varlist)'
            local _cand_names `_cand_names' `_t'
            local _cand_tv    `_cand_tv' `_cc'
        }
        qui _rmcoll `_ctrl_tv' `_cand_tv' if `touse'
        local _surv `r(varlist)'
        local _ivj ""
        local _ii = 0
        foreach _tv of local _cand_tv {
            local ++_ii
            local _nm : word `_ii' of `_cand_names'
            local _ok : list _tv in _surv
            if `_ok' {
                local _pos : list posof "`_nm'" in _ivpool_names
                if `_pos' == 0 {
                    local _ivpool_names `_ivpool_names' `_nm'
                    local _ivpool_tv    `_ivpool_tv' `_tv'
                    local _pos : word count `_ivpool_names'
                }
                local _ptv : word `_pos' of `_ivpool_tv'
                local _ivj `_ivj' `_ptv'
            }
        }
        local iv_for_endog_`j' "`_ivj'"
    }
    local varlist_iv "`_ivpool_tv'"
    local ivlabels   "`_ivpool_names'"
    local count_iv : word count `varlist_iv'
    if `count_iv' < 1 {
        _apexogerr "No usable instruments remain after dropping degenerate/collinear terms."
    }
    if `_raw_iv_n' > `count_iv' {
        local _ndrop = `_raw_iv_n' - `count_iv'
        di as text "  [note] `_ndrop' of `_raw_iv_n' instrument term(s) dropped (base-level, e.g. 0b., or collinear with the controls / other instruments)."
        di as text "         `count_iv' full-rank instrument(s) retained: `ivlabels'"
        di as text "         gmin()/gmax()/mu()/omega() take `count_iv' value(s) -- or a single value, broadcast to all."
    }

    local countmin : word count `gmin'
    local countmax : word count `gmax'
    local count_mu : word count `mu'
    local count_om : word count `omega'

    * Convenience (v0.13.1): a SINGLE gmin/gmax/mu/omega value is broadcast to all
    * instruments.  Handy now that the instrument count is the cleaned, full-rank
    * count (which the user need not pre-compute): `gmin(-0.01) gmax(0.01)' applies
    * the same assumed range to every instrument.
    if `count_iv' > 1 {
        if `countmin' == 1 {
            local _v `gmin'
            local gmin ""
            forvalues k = 1/`count_iv' {
                local gmin "`gmin' `_v'"
            }
            local countmin = `count_iv'
        }
        if `countmax' == 1 {
            local _v `gmax'
            local gmax ""
            forvalues k = 1/`count_iv' {
                local gmax "`gmax' `_v'"
            }
            local countmax = `count_iv'
        }
        if `count_mu' == 1 {
            local _v `mu'
            local mu ""
            forvalues k = 1/`count_iv' {
                local mu "`mu' `_v'"
            }
            local count_mu = `count_iv'
        }
        if `count_om' == 1 {
            local _v `omega'
            local omega ""
            forvalues k = 1/`count_iv' {
                local omega "`omega' `_v'"
            }
            local count_om = `count_iv'
        }
    }

    *==========================================================================
    * (2) Validate options and v0.2 scope restrictions
    *==========================================================================
    if "`estimator'" == "" local estimator "cfprobit"
    if "`estimator'" == "biprobit" {
        _apexogerr "The recursive bivariate probit is requested with firststage(biprobit); estimator() stays at its default cfprobit."
    }
    if "`estimator'" != "cfprobit" {
        _apexogerr "v0.2 supports only estimator(cfprobit); you asked for estimator(`estimator')."
    }

    *--- Resolve per-endog first-stage type ------------------------------------
    * Precedence (highest first):
    *   1. Per-clause cfopts:  (d_k = z_k, <type>)  --  any word matching a type
    *   2. Global firststage() option:
    *        firststage(<type>)            -- applies to all endogs
    *        firststage(<t1> <t2> ... <tJ>) -- one per endog
    *   3. (Phase 1) i./c. prefix on the endog in the clause:
    *        (i.d = ...) -> probit ;  (c.d = ...) -> linear
    *   4. Default: linear
    if "`firststage'" != "" {
        local n_fs : word count `firststage'
        if `n_fs' != 1 & `n_fs' != `count2' {
            _apexogerr "firststage() must have either 1 type (applied to all endogs) or `count2' types (one per endog).  You gave: `firststage'."
        }
        forvalues k = 1/`n_fs' {
            local fst_k : word `k' of `firststage'
            if !inlist("`fst_k'", "linear", "probit", "fprobit", "poisson", "biprobit") {
                _apexogerr "Unknown firststage type '`fst_k''.  Valid: linear, probit, fprobit, poisson, biprobit."
            }
        }
    }
    else {
        local n_fs = 0
    }
    * Resolve per-endog type
    forvalues k = 1/`count2' {
        * (1) Per-clause cfopts override
        local _type_k ""
        if "`cfopts_for_endog_`k''" != "" {
            foreach _w of local cfopts_for_endog_`k' {
                if inlist("`_w'", "linear", "probit", "fprobit", "poisson", "biprobit") {
                    local _type_k "`_w'"
                }
            }
        }
        * (2) Global firststage() option
        if "`_type_k'" == "" & `n_fs' > 0 {
            if `n_fs' == 1 local _type_k "`firststage'"
            else           local _type_k : word `k' of `firststage'
        }
        * (3) Phase-1 i./c. prefix on the endog clause var
        if "`_type_k'" == "" {
            if "`fv_prefix_for_endog_`k''" == "i" local _type_k "probit"
            if "`fv_prefix_for_endog_`k''" == "c" local _type_k "linear"
        }
        * (4) Default
        if "`_type_k'" == "" local _type_k "linear"
        local firststage_`k' "`_type_k'"
    }
    * Display the resolved per-endog types (for the e() / header)
    local firststage_resolved ""
    forvalues k = 1/`count2' {
        local firststage_resolved "`firststage_resolved' `firststage_`k''"
    }
    local firststage_resolved : list clean firststage_resolved
    local firststage "`firststage_resolved'"
    * v0.19.1: flag any generalized-residual first stage (probit/fprobit/
    * poisson).  The sigma_tilde conversion behind gunits(structural) is a
    * linear-first-stage result, so the headers print a warning line when
    * the two are combined.  Display only: no computation uses this flag.
    local _grfs = 0
    foreach _f of local firststage {
        if inlist("`_f'", "probit", "fprobit", "poisson") local _grfs = 1
    }

    *--- v0.19: engine selection ----------------------------------------------
    * firststage(biprobit) replaces the two-step control-function estimator by
    * the recursive bivariate probit: joint MLE of depvar and ONE 0/1
    * endogenous regressor under bivariate normality, with the gamma-offset in
    * the outcome equation.  No control function, no Murphy-Topel step.  The
    * control-function engine ("cf") is byte-identical to v0.18.0.
    local engine    "cf"
    local _est_name "cfprobit"
    local _vce_name "Murphy-Topel"
    local _n_bp = 0
    forvalues k = 1/`count2' {
        if "`firststage_`k''" == "biprobit" local ++_n_bp
    }
    if `_n_bp' > 0 {
        if `count2' > 1 {
            _apexogerr "firststage(biprobit) supports exactly ONE endogenous regressor (the model has `count2'): it fits the recursive bivariate probit of `yvar' and a single 0/1 treatment.  Use the control-function first stages (linear/probit/fprobit/poisson) with several endogenous regressors."
        }
        local engine    "biprobit"
        local _est_name "biprobit"
        local _vce_name "MLE"
    }

    *--- Phase 2: parse target() option ---------------------------------------
    * target() can specify either APE targets or an AIE target:
    *
    *   APE targets:  target(d1)        -> APE of d1
    *                 target(d1 d2)     -> APE of both
    *                 target(c.d1)      -> same (prefix is stripped)
    *
    *   AIE target:   target(c.X#c.Z)   -> cts x cts AIE   (cross-partial)
    *                 target(c.X#i.d)   -> cts x binary AIE (firstdiff)
    *                 target(i.X#i.d)   -> binary x binary  (2nd discrete diff)
    *                 prefixes DRIVE the formula choice (overriding value-based
    *                 auto-detection); both interactors must already appear
    *                 somewhere in the regressor set (exog or endog).
    *
    *   Mixing: a single target() can hold multiple APE specs OR exactly one
    *   AIE spec (not both kinds at once; not multiple AIE specs).
    *
    *   Default (target omitted): every endog in varlist2 is an APE target.
    *
    *   aie() and target() are mutually exclusive; aie() is the back-compat
    *   shorthand and gets internally translated to a target() AIE spec below.
    forvalues k = 1/`count2' {
        local target_endog_`k' = 1     // default: every endog is a target
    }
    * AIE-from-target() info, populated only if target() carries a # spec
    local target_is_aie = 0
    local target_aie_spec ""
    local target_aie_X    ""
    local target_aie_Z    ""
    local target_aie_X_pfx ""
    local target_aie_Z_pfx ""

    if "`target'" != "" {
        if "`aie'" != "" {
            _apexogerr "target() and aie() are mutually exclusive.  Use target(c.X#i.d) or aie(X d) but not both."
        }
        * First pass: count AIE specs and APE specs to enforce one-of-each rule.
        local _n_aie_specs = 0
        local _n_ape_specs = 0
        foreach _spec of local target {
            if strpos("`_spec'", "#") local ++_n_aie_specs
            else                       local ++_n_ape_specs
        }
        if `_n_aie_specs' > 1 {
            _apexogerr "target() supports at most one AIE spec (e.g. target(c.X#i.d)); got `_n_aie_specs'."
        }
        if `_n_aie_specs' > 0 & `_n_ape_specs' > 0 {
            _apexogerr "target() cannot mix APE specs (e.g. d1) and AIE specs (e.g. c.X#i.d) in the same call."
        }

        if `_n_aie_specs' == 1 {
            *--- AIE branch via target() ---------------------------------------
            local target_is_aie = 1
            local target_aie_spec : word 1 of `target'
            * Strip any leading equation prefix (Stata's display of e(b) names
            * can include eqn prefixes; user shouldn't normally write these).
            * Split spec on '#'.  Must have exactly 2 components.
            local _spec_raw "`target_aie_spec'"
            local _hash_pos = strpos("`_spec_raw'", "#")
            local _left  = substr("`_spec_raw'", 1, `_hash_pos' - 1)
            local _right = substr("`_spec_raw'", `_hash_pos' + 1, .)
            * Each side may have c./i. prefix.
            local _left_pfx ""
            local _left_bare "`_left'"
            if substr("`_left'", 1, 2) == "i." {
                local _left_pfx "i"
                local _left_bare = substr("`_left'", 3, .)
            }
            else if substr("`_left'", 1, 2) == "c." {
                local _left_pfx "c"
                local _left_bare = substr("`_left'", 3, .)
            }
            local _right_pfx ""
            local _right_bare "`_right'"
            if substr("`_right'", 1, 2) == "i." {
                local _right_pfx "i"
                local _right_bare = substr("`_right'", 3, .)
            }
            else if substr("`_right'", 1, 2) == "c." {
                local _right_pfx "c"
                local _right_bare = substr("`_right'", 3, .)
            }
            * Validate both are in the regressor set (unabbreviate first).
            unab _left_bare  : `_left_bare'
            unab _right_bare : `_right_bare'
            local _allvars `varlist1' `varlist2'
            local _xin : list _left_bare in _allvars
            local _zin : list _right_bare in _allvars
            if !`_xin' _apexogerr "target(): interactor '`_left_bare'' is not in the regressor list (`_allvars')."
            if !`_zin' _apexogerr "target(): interactor '`_right_bare'' is not in the regressor list (`_allvars')."

            local target_aie_X     "`_left_bare'"
            local target_aie_Z     "`_right_bare'"
            local target_aie_X_pfx "`_left_pfx'"
            local target_aie_Z_pfx "`_right_pfx'"
            * Display label stays in fv form for the user.
            local target_aie_spec  "`_left'#`_right'"

            * Disable per-endog APE display (we're in AIE mode).
            forvalues k = 1/`count2' {
                local target_endog_`k' = 0
            }
            local target_list ""
        }
        else {
            *--- APE branch via target() (Phase 1 behavior) --------------------
            forvalues k = 1/`count2' {
                local target_endog_`k' = 0
            }
            local _target_clean ""
            foreach _spec of local target {
                local _bare "`_spec'"
                if substr("`_spec'", 1, 2) == "i." local _bare = substr("`_spec'", 3, .)
                if substr("`_spec'", 1, 2) == "c." local _bare = substr("`_spec'", 3, .)
                capture unab _bare : `_bare'
                local _matched = 0
                local _k = 0
                foreach _e of local varlist2 {
                    local ++_k
                    if "`_e'" == "`_bare'" {
                        local target_endog_`_k' = 1
                        local _matched = 1
                        local _target_clean "`_target_clean' `_bare'"
                        continue, break
                    }
                }
                if !`_matched' {
                    _apexogerr "target() APE spec '`_spec'' (bare: `_bare') is not in the endogenous regressor list (`varlist2')."
                }
            }
            local target_list : list clean _target_clean
        }
    }
    else {
        * Default: all endogs are targets; build display list
        local target_list "`varlist2'"
    }
    * Count active APE targets (irrelevant when target_is_aie=1).
    local n_active_targets = 0
    forvalues k = 1/`count2' {
        if `target_endog_`k'' local ++n_active_targets
    }
    if !`target_is_aie' & `n_active_targets' < 1 {
        _apexogerr "Internal: no active APE targets after parsing target() option."
    }

    *--- vce() is parsed in step (1a) above, before the sample is marked. ------

    if `count2' < 1 {
        _apexogerr "At least one endogenous regressor required."
    }
    if `count_iv' < `count2' {
        _apexogerr "Need at least as many instruments (`count_iv') as endogenous regressors (`count2')."
    }
    * v0.6: multi-endog supported.  Each endogenous regressor gets its own
    * first stage on the same shared instrument set, its own chat0, and its
    * own A_per_obs.  The APE is reported for the FIRST endog (reorder
    * `varlist2' to change which one is the target).

    if "`method'" == "uci" {
        if `countmin' != `count_iv' | `countmax' != `count_iv' {
            _apexogerr "For UCI, gmin() and gmax() each need `count_iv' element(s)."
        }
        local k = 1
        foreach n of numlist `gmin' {
            local gmin`k' = `n'
            local ++k
        }
        local k = 1
        foreach n of numlist `gmax' {
            local gmax`k' = `n'
            local ++k
        }
        forvalues k = 1/`count_iv' {
            if `gmin`k'' > `gmax`k'' {
                _apexogerr "gmin() must not exceed gmax(): instrument `k' has gmin = `gmin`k'' > gmax = `gmax`k''."
            }
        }
        if `grid' < 2 _apexogerr "grid() must be at least 2."
    }

    if "`method'" == "ltz" {
        if "`distribution'" != "" {
            if "`mu'" != "" | "`omega'" != "" {
                _apexogerr "distribution() is mutually exclusive with mu()/omega()."
            }
            * Set ltz_mode = "distribution"; gammaCall1..gammaCallK strings will
            * hold the per-IV random-draw expression to evaluate inside the LTZ loop.
            local ltz_mode "distribution"

            * Tokenise: distribution(<type> p1 p2 ...) where the type is the
            * first word and the remainder are parameters.  Stata's numlist for
            * the parameter portion is dropped here (we want the raw string),
            * so we re-parse from the original `distribution' macro.
            local distribution: subinstr local distribution "," " , ", all
            local distcnt : list sizeof distribution
            local jj = 1
            forvalues j = 1/`distcnt' {
                local dist`jj' : word `j' of `distribution'
                local dist`jj' : subinstr local dist`jj' "," ""
                if "`dist`jj''" != "" local ++jj
            }
            * Now: dist1 = type, dist2... = parameters; jj-1 = total non-empty tokens
            local distname = "`dist1'"
            local expD = 2 + 2*`count_iv'      // tokens expected for 2-param dists
            local expS = 2 + 1*`count_iv'      // tokens expected for 1-param dists
            local cntD = 2*`count_iv'
            local cntS = 1*`count_iv'
            local accept "normal, uniform, chi2, poisson, t, gamma, special"

            if "`distname'" == "normal" {
                if `jj' != `expD' _apexogerr "distribution(normal ...) with `count_iv' IV(s) needs `cntD' parameters (mean and sd per IV)."
                forvalues k = 1/`count_iv' {
                    local kl = `k'*2
                    local kh = `k'*2+1
                    local gammaCall`k' = "rnormal(`dist`kl'', `dist`kh'')"
                }
            }
            else if "`distname'" == "uniform" {
                if `jj' != `expD' _apexogerr "distribution(uniform ...) with `count_iv' IV(s) needs `cntD' parameters (min and max per IV)."
                forvalues k = 1/`count_iv' {
                    local kl = `k'*2
                    local kh = `k'*2+1
                    local gammaCall`k' = "`dist`kl''+(`dist`kh''-`dist`kl'')*runiform()"
                }
            }
            else if "`distname'" == "chi2" {
                if `jj' != `expS' _apexogerr "distribution(chi2 ...) with `count_iv' IV(s) needs `cntS' parameter(s) (df per IV)."
                forvalues k = 1/`count_iv' {
                    local kl = `k'+1
                    if `dist`kl'' < 1 _apexogerr "chi2 df must be >= 1; got `dist`kl''."
                    local gammaCall`k' = "rchi2(`dist`kl'')"
                }
            }
            else if "`distname'" == "poisson" {
                if `jj' != `expS' _apexogerr "distribution(poisson ...) with `count_iv' IV(s) needs `cntS' parameter(s) (mean per IV)."
                forvalues k = 1/`count_iv' {
                    local kl = `k'+1
                    if `dist`kl'' < 1 _apexogerr "poisson mean must be >= 1; got `dist`kl''."
                    local gammaCall`k' = "rpoisson(`dist`kl'')"
                }
            }
            else if "`distname'" == "t" {
                if `jj' != `expS' _apexogerr "distribution(t ...) with `count_iv' IV(s) needs `cntS' parameter(s) (df per IV)."
                forvalues k = 1/`count_iv' {
                    local kl = `k'+1
                    if `dist`kl'' < 1 _apexogerr "t df must be >= 1; got `dist`kl''."
                    local gammaCall`k' = "rt(`dist`kl'')"
                }
            }
            else if "`distname'" == "gamma" {
                if `jj' != `expD' _apexogerr "distribution(gamma ...) with `count_iv' IV(s) needs `cntD' parameters (shape and scale per IV)."
                forvalues k = 1/`count_iv' {
                    local kl = `k'*2
                    local kh = `k'*2+1
                    if `dist`kl''<=0 | `dist`kh''<=0 _apexogerr "gamma shape and scale must be > 0."
                    local gammaCall`k' = "rgamma(`dist`kl'',`dist`kh'')"
                }
            }
            else if "`distname'" == "special" {
                * Each per-IV "parameter" is the NAME of a variable whose
                * empirical distribution we sample from.  We validate names
                * here but defer the donor-matrix construction to the LTZ
                * section (which runs after `marksample touse').
                if `jj' != `expS' _apexogerr "distribution(special ...) needs one variable name per IV (the empirical-distribution donor)."
                local special_donor_vars ""
                forvalues k = 1/`count_iv' {
                    local kl = `k'+1
                    capture confirm variable `dist`kl''
                    if _rc _apexogerr "distribution(special): variable `dist`kl'' not found."
                    local special_donor_vars `special_donor_vars' `dist`kl''
                }
                * gammaCall`k' will be built later (just before the LTZ loop)
                * once the donor matrix exists.
                local special_pending = 1
            }
            else {
                _apexogerr "distribution(): type must be one of: `accept'.  Got '`distname''."
            }
        }
        else {
            * Gaussian (legacy) path
            local ltz_mode "gaussian"
            if `count_mu' != `count_iv' | `count_om' != `count_iv' {
                _apexogerr "For LTZ, mu() and omega() each need `count_iv' element(s)  (or use distribution() instead)."
            }
            local k = 1
            foreach n of numlist `mu' {
                local mu`k' = `n'
                local ++k
            }
            local k = 1
            foreach n of numlist `omega' {
                local omega`k' = `n'
                if `omega`k'' <= 0 _apexogerr "All omega() values must be > 0."
                local ++k
            }
            * Build gammaCall locals so the LTZ loop is unified across both modes
            forvalues k = 1/`count_iv' {
                local _sd_k = sqrt(`omega`k'')
                local gammaCall`k' = "`mu`k'' + `_sd_k'*rnormal()"
            }
        }
        if `iterations' < 10 _apexogerr "iterations() must be >= 10."
        if "`seed'" != "" set seed `seed'
    }

    local dvar : word 1 of `varlist2'

    *==========================================================================
    * (3) Determine the estimation sample
    *==========================================================================
    * (`touse' was fixed in step (1b): if/in + missing values, incl. the
    *  cluster id, so the cleaning and the fits share one sample.)
    qui count if `touse'
    local N = r(N)
    if `N' == 0 _apexogerr "No observations remain after if/in and missing-data filtering."
    * Guard: the outcome must be coded 0/1.  Stata's probit treats any nonzero
    * value as a positive outcome, but the Mata score uses the raw value, so a
    * differently coded outcome would give silently wrong two-step SEs.
    qui count if `touse' & !inlist(`yvar', 0, 1)
    if r(N) > 0 {
        local _nbad = r(N)
        _apexogerr "Dependent variable `yvar' takes values other than 0/1 in `_nbad' observation(s) of the estimation sample.  Recode it first, e.g. gen byte y01 = (`yvar' != 0)."
    }
    if "`vcetype'" == "cluster" {
        mata: st_numscalar("_apexogNCLUST", rows(uniqrows(st_data(., "`clustvar'", "`touse'"))))
    }

    *==========================================================================
    * (3b) Process aie() option, if specified.
    *      aie(X Z) requests the Average Interaction Effect for the c.X##c.Z
    *      interaction.  Either or both interactors may be ENDOGENOUS -- this
    *      is the canonical cfprobit use case (interaction of an endog
    *      treatment with an exog moderator).  We generate a tempvar
    *      XZ = X*Z and add it to varlist1 so the main probit fits with the
    *      interaction.  The user should NOT have written `c.X##c.Z'
    *      separately in the varlist; let apexog add the interaction
    *      term automatically.
    *
    *      IMPORTANT: when X or Z is endogenous, XZ = X*Z is also endogenous,
    *      so it must be EXCLUDED from that endog's first-stage regression
    *      (otherwise the first stage is circular).  The per-endog exog mask
    *      below handles this.
    *
    *      v0.9 supports continuous-by-continuous AIE only.  A binary
    *      interactor would need a discrete-change formulation.
    *==========================================================================
    * Two ways into AIE mode now:
    *   (i)  aie(X Z)            -- back-compat shorthand; auto-detect types
    *   (ii) target(c.X#i.d ...) -- new Phase 2 syntax; prefixes drive types
    * Both routes converge on the same internal state: aie_active=1, aie_X,
    * aie_Z, aie_X_binary, aie_Z_binary, aie_case, and the aie_XZ tempvar in
    * varlist1.  How the type flags are filled differs:
    *
    *   - aie() path:    type flags filled by value-based detection on data
    *   - target() path: type flags filled from c./i. prefix on each spec;
    *                    a missing prefix falls back to value-based detection
    local aie_active = 0
    local aie_xz_is_endog = 0
    if "`aie'" != "" | `target_is_aie' {
        if "`aie'" != "" {
            local n_aie : word count `aie'
            if `n_aie' != 2 _apexogerr "aie() takes exactly 2 variable names; you gave: `aie'"
            local aie_X : word 1 of `aie'
            local aie_Z : word 2 of `aie'
            * Defer to value-based detection for both
            local _aie_X_pfx ""
            local _aie_Z_pfx ""
            local aie_source "aie()"
        }
        else {
            * target() AIE path
            local aie_X     "`target_aie_X'"
            local aie_Z     "`target_aie_Z'"
            local _aie_X_pfx "`target_aie_X_pfx'"
            local _aie_Z_pfx "`target_aie_Z_pfx'"
            local aie_source "target(`target_aie_spec')"
        }
        unab aie_X : `aie_X'
        unab aie_Z : `aie_Z'
        confirm numeric variable `aie_X'
        confirm numeric variable `aie_Z'
        if "`aie_X'" == "`aie_Z'" {
            _apexogerr "The two AIE interactors must be distinct variables; got `aie_X' twice."
        }

        * Validate both interactors are somewhere in the regressor set.
        local _allvars `varlist1' `varlist2'
        local _xin : list aie_X in _allvars
        local _zin : list aie_Z in _allvars
        if !`_xin' _apexogerr "AIE interactor `aie_X' must appear in the regressor list (exog or endog)."
        if !`_zin' _apexogerr "AIE interactor `aie_Z' must appear in the regressor list (exog or endog)."

        * Resolve type per interactor.  Precedence:
        *   1. explicit prefix from target() (i. -> binary, c. -> cts)
        *   2. value-based auto-detection (values strictly in {0,1} -> binary)
        local aie_X_binary = .
        local aie_Z_binary = .
        if "`_aie_X_pfx'" == "i" local aie_X_binary = 1
        if "`_aie_X_pfx'" == "c" local aie_X_binary = 0
        if "`_aie_Z_pfx'" == "i" local aie_Z_binary = 1
        if "`_aie_Z_pfx'" == "c" local aie_Z_binary = 0
        if `aie_X_binary' == . {
            qui count if `touse' & !inlist(`aie_X', 0, 1)
            local aie_X_binary = cond(r(N) == 0, 1, 0)
        }
        if `aie_Z_binary' == . {
            qui count if `touse' & !inlist(`aie_Z', 0, 1)
            local aie_Z_binary = cond(r(N) == 0, 1, 0)
        }

        * Case codes:
        *   0 = continuous-by-continuous          AIE = E[phi*(b_XZ - eta*a*b)]
        *   1 = continuous-by-binary (X cts, Z bin) firstdiff in Z, slope in X
        *   2 = binary-by-continuous (X bin, Z cts) firstdiff in X, slope in Z
        *   3 = binary-by-binary                  AIE = E[2nd discrete diff of Phi]
        local aie_case = 2*`aie_X_binary' + `aie_Z_binary'

        local _aie_case_lbl ""
        if `aie_case' == 0 local _aie_case_lbl "continuous-by-continuous (cross-partial)"
        if `aie_case' == 1 local _aie_case_lbl "continuous-by-binary (firstdiff in `aie_Z')"
        if `aie_case' == 2 local _aie_case_lbl "binary-by-continuous (firstdiff in `aie_X')"
        if `aie_case' == 3 local _aie_case_lbl "binary-by-binary (2nd discrete difference)"

        * v0.13: if this interaction was DECLARED as an endogenous regressor in
        * the varlist (the endogenous-treatment x moderator case), target that
        * existing endog tempvar and let it carry its own control function.  Do
        * NOT create a separate exogenous control.  Matching is on the unordered
        * operand pair.  Otherwise (all-exogenous, or undeclared) fall back to
        * the original behaviour: materialise the product as an exog control.
        if "`aieendog_tempvar'" != ""                                          ///
           & ( ("`aie_X'"=="`aieendog_X'" & "`aie_Z'"=="`aieendog_Z'")         ///
             | ("`aie_X'"=="`aieendog_Z'" & "`aie_Z'"=="`aieendog_X'") ) {
            local aie_xz_is_endog = 1
            local aie_XZ "`aieendog_tempvar'"
        }
        else {
            * Create the interaction tempvar.  For binary-by-binary this is also
            * the `1.X#1.Z' indicator (since X,Z in {0,1}, X*Z=1 iff both=1).
            tempvar aie_XZ
            qui gen double `aie_XZ' = `aie_X' * `aie_Z' if `touse'
            local varlist1 `varlist1' `aie_XZ'
            local count1 = `count1' + 1
        }
        local aie_active = 1

        di as text "  [aie] target = `_aie_case_lbl'  (via `aie_source')"

        * Diagnostic note when an interactor is endogenous
        local _xin_endog : list aie_X in varlist2
        local _zin_endog : list aie_Z in varlist2
        if `_xin_endog' | `_zin_endog' {
            di as text "  [note] AIE: one or both interactors are endogenous; the MT correction propagates first-stage uncertainty into the AIE SE."
        }
    }
    else {
        local aie_X_binary = 0
        local aie_Z_binary = 0
        local aie_case = 0
    }

    if "`engine'" == "biprobit" & `aie_active' {
        _apexogerr "firststage(biprobit) targets the APE of the single 0/1 endogenous regressor; the AIE target (aie() / target(c.X#i.d)) is not available with the bivariate-probit engine."
    }

    *==========================================================================
    * cfinteract: include products of the control functions in the main probit
    * (replicates cfprobit's cfinteract = c.chat1##c.chat2... minus the mains).
    * Auto-on when the AIE interaction is itself an endogenous regressor (the
    * canonical endogenous-treatment x moderator model); also honored as an
    * explicit option for general multi-endog use.  Needs >=2 control functions.
    *==========================================================================
    local do_cfinteract = 0
    if "`cfinteract'" != "" local do_cfinteract = 1
    if `aie_xz_is_endog'    local do_cfinteract = 1
    if `count2' < 2         local do_cfinteract = 0
    if "`cfinteract'" != "" & `count2' < 2 {
        di as text "  [note] cfinteract needs at least two endogenous regressors; ignored."
    }

    *==========================================================================
    * Short labels for matrix column names and e(lb_*)/e(ub_*) scalar names.
    * Stata names are capped at 32 characters, so a long endog/interactor name
    * would otherwise error at the very end of the run.  Short names are
    * unchanged (e(lb_d), ape_d, aie_X_Z, ... exactly as before).
    *==========================================================================
    local _all_labs ""
    local j = 0
    foreach endogvar of local varlist2 {
        local ++j
        local _lab_`j' = substr("`endogvar'", 1, 27)
        local _dup : list posof "`_lab_`j''" in _all_labs
        if `_dup' > 0 local _lab_`j' "endog`j'"
        local _all_labs "`_all_labs' `_lab_`j''"
    }
    if `aie_active' {
        local _labX = substr("`aie_X'", 1, 13)
        local _labZ = substr("`aie_Z'", 1, 13)
    }

    *==========================================================================
    * (4) Run the first stage ONCE (does not depend on gamma).  Save:
    *       - chat0     = control function (residual or generalized residual)
    *       - Wpi_hat   = first-stage linear predictor (used by Mata for the
    *                     per-obs MT weight A_i in nonlinear first stages)
    *       - A_per_obs = per-obs MT weight = -d(chat0)/d(Wpi): 1 for linear,
    *                     chat0*(Wpi+chat0) for probit, the same plus
    *                     d(1-d)*[phi/(Phi(1-Phi))]^2 for fprobit, exp(Wpi)
    *                     for poisson
    *==========================================================================
    * Build the vce() option string that goes into regress/probit/etc.
    if "`vcetype'" == "robust" local vce_opt = "vce(robust)"
    else                        local vce_opt = "vce(cluster `clustvar')"

    * Per-endog first-stage exog list.  Identical to varlist1 EXCEPT for
    * endogs that are AIE interactors -- those exclude the aie_XZ tempvar
    * (since X*d on the RHS of a regression of d would be circular).
    local j = 0
    foreach dj of local varlist2 {
        local ++j
        local exog_for_endog_`j' "`varlist1'"
        if `aie_active' {
            if "`dj'" == "`aie_X'" | "`dj'" == "`aie_Z'" {
                local exog_for_endog_`j' : list exog_for_endog_`j' - aie_XZ
            }
        }
    }

    * For each endogenous regressor, build its own chat0 and A_per_obs by
    * running an independent first-stage regression.  v0.7: each endog can
    * have its own first-stage type via `firststage_<j>'.  v0.9: each endog
    * can have its own first-stage EXOG list via `exog_for_endog_<j>' (used
    * to drop the AIE tempvar when the endog itself is an interactor).
    local chat0_list ""
    local A_list     ""
    * v0.21.2 (display only; no computation uses either):
    *   _binlin        0/1 endogenous regressors that get a linear first stage
    *   _apexogGSMAX   |pi_z|/sigma_v, defined for ONE continuous (not 0/1)
    *                  endogenous regressor with ONE instrument and a linear
    *                  first stage; missing otherwise
    local _binlin ""
    scalar _apexogGSMAX = .
    local j = 0
    foreach dj of local varlist2 {
        local ++j
        tempvar chat0_`j' Wpi_`j' A_`j'
        local fst = "`firststage_`j''"
        * Guard: the first-stage model must suit the endog's support (probit on
        * a continuous variable used to die with r(2000) "outcome does not
        * vary", since probit treats every nonzero value as 1).
        if "`fst'" == "probit" {
            qui count if `touse' & !inlist(`dj', 0, 1)
            if r(N) > 0 _apexogerr "firststage(probit) requires a 0/1 endogenous regressor, but `dj' has other values (use fprobit for a fractional or linear for a continuous variable)."
        }
        else if "`fst'" == "fprobit" {
            qui count if `touse' & (`dj' < 0 | `dj' > 1)
            if r(N) > 0 _apexogerr "firststage(fprobit) requires an endogenous regressor in [0,1], but `dj' has values outside that range."
        }
        else if "`fst'" == "poisson" {
            qui count if `touse' & `dj' < 0
            if r(N) > 0 _apexogerr "firststage(poisson) requires a nonnegative endogenous regressor, but `dj' has negative values."
        }
        if "`fst'" == "biprobit" {
            * v0.19: no control function -- the recursive bivariate probit is
            * fit jointly in _apexogfit at every gamma.  Only the support check.
            qui count if `touse' & !inlist(`dj', 0, 1)
            if r(N) > 0 _apexogerr "firststage(biprobit) requires a 0/1 endogenous regressor, but `dj' has other values."
            continue
        }
        if "`fst'" == "linear" {
            qui regress `dj' `iv_for_endog_`j'' `exog_for_endog_`j'' if `touse', `vce_opt'
            qui predict double `chat0_`j'' if `touse', residuals
            qui predict double `Wpi_`j''   if `touse', xb
            qui gen double `A_`j'' = 1 if `touse'
            * v0.21.2 (display only).  A 0/1 regressor cannot have the normal
            * linear first stage of the control-function model: flag it for the
            * header note (the product regressor of an endogenous interaction
            * is skipped; its treatment variable carries the flag).
            capture assert inlist(`dj', 0, 1) if `touse'
            local _dj01 = (_rc == 0)
            if `_dj01' & "`dj'" != "`aieendog_tempvar'" local _binlin "`_binlin' `dj'"
            * The limit of the structural scale, |pi_z|/sigma_v (help apexog,
            * gunits(structural)): one CONTINUOUS endogenous regressor, one
            * instrument.  Not defined for a 0/1 regressor, for which the linear
            * normal first stage behind the conversion cannot hold.
            local _niv_j : word count `iv_for_endog_`j''
            if `count2' == 1 & `count_iv' == 1 & `_niv_j' == 1 & !`_dj01' {
                local _ivone : word 1 of `iv_for_endog_`j''
                capture scalar _apexogGSMAX = abs(_b[`_ivone'])
                if _rc scalar _apexogGSMAX = .
                qui summarize `chat0_`j'' if `touse'
                if r(sd) > 0 & r(sd) < . scalar _apexogGSMAX = scalar(_apexogGSMAX) / r(sd)
                else                     scalar _apexogGSMAX = .
            }
        }
        else if "`fst'" == "probit" {
            qui probit `dj' `iv_for_endog_`j'' `exog_for_endog_`j'' if `touse', `vce_opt'
            qui predict double `Wpi_`j'' if `touse', xb
            qui gen double `chat0_`j'' = (`dj' - normal(`Wpi_`j'')) * normalden(`Wpi_`j'') ///
                                         / (normal(`Wpi_`j'') * normal(-`Wpi_`j''))     ///
                                         if `touse'
            qui gen double `A_`j'' = `chat0_`j'' * (`Wpi_`j'' + `chat0_`j'') if `touse'
        }
        else if "`fst'" == "fprobit" {
            qui fracreg probit `dj' `iv_for_endog_`j'' `exog_for_endog_`j'' if `touse', `vce_opt'
            qui predict double `Wpi_`j'' if `touse', xb
            qui gen double `chat0_`j'' = (`dj' - normal(`Wpi_`j'')) * normalden(`Wpi_`j'') ///
                                         / (normal(`Wpi_`j'') * normal(-`Wpi_`j''))     ///
                                         if `touse'
            * v0.21.1: -d(chat0)/d(Wpi) for a FRACTIONAL d.  The identity
            * chat0*(Wpi+chat0) of the binary probit needs d in {0,1}; for
            * d in (0,1) the derivative has the extra term d(1-d)*lam^2,
            * lam = phi(Wpi)/[Phi(Wpi)(1-Phi(Wpi))] (zero when d is 0 or 1).
            qui gen double `A_`j'' = `chat0_`j'' * (`Wpi_`j'' + `chat0_`j'')               ///
                + `dj' * (1 - `dj') * (normalden(`Wpi_`j'')                               ///
                                       / (normal(`Wpi_`j'') * normal(-`Wpi_`j'')))^2      ///
                if `touse'
        }
        else if "`fst'" == "poisson" {
            qui poisson `dj' `iv_for_endog_`j'' `exog_for_endog_`j'' if `touse', `vce_opt'
            qui predict double `Wpi_`j'' if `touse', xb
            qui gen double `chat0_`j'' = `dj' - exp(`Wpi_`j'') if `touse'
            qui gen double `A_`j''     = exp(`Wpi_`j'') if `touse'
        }
        * Guard: the first stage must have used the full estimation sample (a
        * probit first stage drops observations under perfect prediction) and
        * the control function must be defined everywhere; otherwise the Mata
        * views (built on the full sample) and the first-stage fit would be
        * misaligned and the two-step variance would be wrong.
        if e(N) < `N' {
            local _nfs = e(N)
            _apexogerr "The `fst' first stage for `dj' used `_nfs' of `N' observations (perfect prediction or dropped groups); apexog cannot proceed with a misaligned first stage."
        }
        qui count if `touse' & missing(`chat0_`j'')
        if r(N) > 0 {
            local _nmis = r(N)
            _apexogerr "The control function for `dj' is undefined in `_nmis' observation(s) (extreme first-stage fitted values); apexog cannot proceed."
        }
        local chat0_list `chat0_list' `chat0_`j''
        local A_list     `A_list'     `A_`j''
    }

    *==========================================================================
    * (4d) cfinteract product terms: all unique pairwise products of the control
    *      functions.  For the canonical 2-endog AIE this is the single product
    *      chat_d * chat_dX.  These enter the main probit AFTER the chat0s (so
    *      e(b) order is [exog endog chat0 cfprod _cons]).
    *==========================================================================
    local cfprod_list ""
    if `do_cfinteract' {
        local _nc : word count `chat0_list'
        local _amax = `_nc' - 1
        forvalues a = 1/`_amax' {
            local _ca : word `a' of `chat0_list'
            local _bmin = `a' + 1
            forvalues b = `_bmin'/`_nc' {
                local _cb : word `b' of `chat0_list'
                tempvar _cfp_`a'_`b'
                qui gen double `_cfp_`a'_`b'' = `_ca' * `_cb' if `touse'
                local cfprod_list `cfprod_list' `_cfp_`a'_`b''
            }
        }
        local _ncf : word count `cfprod_list'
        di as text "  [cfinteract] added `_ncf' control-function product term(s) to the main probit."
    }

    * v0.6: all chat0s enter the main probit as regressors.  The MT correction
    * accumulates contributions from EVERY first stage.  The APE target is the
    * FIRST endogenous regressor (the rest are treated as controls).
    if `count2' > 1 {
        di as text "  [note] `count2' endogenous regressors detected.  All get a CF in the" /*
            */ " main probit; the MT correction accumulates uncertainty from all first stages." /*
            */ "  Bounds are reported for every endog (or those named in target())."
    }

    if "`debug'" != "" {
        di as text "[debug] N = `N',  exog = `varlist1',  endog = `varlist2',  IVs = `varlist_iv'"
        di as text "[debug] First stage(s) = `firststage'.  chat0_list = `chat0_list'"
        forvalues j = 1/`count2' {
            local endog_j : word `j' of `varlist2'
            di as text "[debug] endog `endog_j' IVs: `iv_for_endog_`j''"
        }
    }

    *==========================================================================
    * (4b) Build the IV-mask matrix: J x n_iv with 1 if IV k is in endog j's
    *      first-stage regression, 0 otherwise.  Passed to Mata so the per-
    *      endog MT correction uses only the relevant IVs in M_pi_k.
    *==========================================================================
    matrix _apexogIVMASK = J(`count2', `count_iv', 0)
    local j = 0
    foreach endogvar of local varlist2 {
        local ++j
        local kcol = 0
        foreach ivvar of local varlist_iv {
            local ++kcol
            local _in : list ivvar in iv_for_endog_`j'
            if `_in' matrix _apexogIVMASK[`j', `kcol'] = 1
        }
    }

    *==========================================================================
    * (4c) Build the EXOG-mask matrix: J x count1 with 1 if exog k is in endog
    *      j's first-stage regression, 0 otherwise.  Default: all 1s.  When
    *      aie() is active and endog j is an interactor, the aie_XZ column is
    *      masked out (it would be circular).  Passed to Mata for W_k.
    *==========================================================================
    matrix _apexogEXOGMASK = J(`count2', `count1', 1)
    if `aie_active' {
        * Find the column position of aie_XZ within varlist1
        local _xz_pos = 0
        local _pos    = 0
        foreach v of local varlist1 {
            local ++_pos
            if "`v'" == "`aie_XZ'" {
                local _xz_pos = `_pos'
                continue, break
            }
        }
        if `_xz_pos' > 0 {
            local j = 0
            foreach dj of local varlist2 {
                local ++j
                if "`dj'" == "`aie_X'" | "`dj'" == "`aie_Z'" {
                    matrix _apexogEXOGMASK[`j', `_xz_pos'] = 0
                }
            }
        }
    }

    *==========================================================================
    * (5) Allocate the per-iteration offset variable (filled per gamma)
    *==========================================================================
    tempvar offvar
    qui gen double `offvar' = 0 if `touse'

    * v0.20: which (cleaned) instruments are 0/1 in the estimation sample?
    * Decides how _apexoggpr reports the exact Pr-scale value of a violation
    * (average discrete change for 0/1 instruments, average derivative else).
    local ivbin ""
    foreach z of local varlist_iv {
        capture assert `z' == 0 | `z' == 1 if `touse'
        local ivbin "`ivbin' `=(_rc == 0)'"
    }

    *==========================================================================
    * (6) Baseline fit (gamma = 0)
    *==========================================================================
    qui replace `offvar' = 0 if `touse'
    * v0.19: the biprobit engine warm-starts every later fit from the baseline
    * e(b); make sure no stale matrix from an interrupted run is picked up.
    capture matrix drop _apexogBP_B0
    _apexogfit, touse(`touse')                                                ///
        yvar(`yvar') endoglist(`varlist2') chat0list(`chat0_list')          ///
        offvar(`offvar') aweightlist(`A_list')                              ///
        exoglist(`varlist1') ivlist(`varlist_iv')                           ///
        clustvar(`clustvar') vcetype(`vcetype')                             ///
        aiex(`aie_X') aiez(`aie_Z') aiexz(`aie_XZ') aiecase(`aie_case') cfprodlist(`cfprod_list') engine(`engine')                         ///
        level(`level') `debug'
    * First-endog scalars (back-compat)
    local baseline_ape = r(ape)
    local baseline_se  = r(se)
    local baseline_lo  = r(lo)
    local baseline_hi  = r(hi)
    * Per-endog baseline vectors (J x 1 each)
    matrix _apexogBASE_APE = r(ape_vec)
    matrix _apexogBASE_SE  = r(se_vec)
    matrix _apexogBASE_LO  = r(lo_vec)
    matrix _apexogBASE_HI  = r(hi_vec)

    *==========================================================================
    * (6b) gamma scale.  m_phi = mean(normalden(eta_hat)) at the baseline fit
    *      is the average dPr/dindex slope: 1 latent-index unit of gamma moves
    *      Pr(y=1) by m_phi per unit of Z, on average.  With gunits(pr), all
    *      user-supplied gamma quantities are PROBABILITY POINTS and are
    *      converted to latent units here, AFTER the baseline fit (which is
    *      gamma-free) and BEFORE the UCI/LTZ loops.  The conversion is
    *      first-order (exact as gamma -> 0).  Everything downstream --
    *      offsets, e(grid), graphs, displayed gamma values -- is latent.
    *
    *      With gunits(structural), user gamma quantities are on the
    *      STRUCTURAL latent scale (total error variance 1, the scale a
    *      researcher elicits when reasoning about y* = xb + gamma*z + eps
    *      with sd(eps) = 1) and are divided by sigma_tilde = sd(eps | CF) =
    *      (1 + b_cf' Var(chat) b_cf)^(-1/2), estimated at the baseline fit.
    *      The probit index normalizes the CONDITIONAL error to 1, so index
    *      quantities = structural quantities / sigma_tilde.  Same machinery
    *      as gunits(pr) with a different conversion factor.
    *      v0.19.1: this conversion is derived for a LINEAR first stage.  For
    *      the generalized-residual first stages (probit/fprobit/poisson) the
    *      validation Monte Carlos put the correct
    *      offset at the UNSCALED structural gamma, so the headers warn when
    *      gunits(structural) meets one of them (flag _grfs); the arithmetic
    *      below is deliberately unchanged.
    *==========================================================================
    scalar _apexogMPHI0   = r(mphi)
    scalar _apexogSTILDE0 = r(stilde)
    * Cross-target MT covariance of the baseline estimates (for ereturn post)
    matrix _apexogVTGT0 = r(vtgt)
    if "`engine'" == "biprobit" {
        * baseline endogeneity correlation of the recursive bivariate probit
        capture scalar _apexogBPRHO0 = scalar(_apexogBPRHO)
        if _rc scalar _apexogBPRHO0 = .
    }
    * v0.21.2 (display only): the user's box in STRUCTURAL units, kept for the
    * header warning (gmin/gmax are overwritten with index values just below).
    local _gs_boxmax = .
    if "`gunits'" == "structural" & "`method'" == "uci" & `count_iv' == 1 {
        local _gs_boxmax = max(abs(`gmin1'), abs(`gmax1'))
    }
    if "`gunits'" == "pr" | "`gunits'" == "structural" {
        if "`gunits'" == "pr" {
            if _apexogMPHI0 >= . {
                _apexogerr "gunits(pr) needs m_phi from a successful baseline fit, but the baseline failed."
            }
            local _ginv : di %21.15g 1/scalar(_apexogMPHI0)
        }
        else {
            if _apexogSTILDE0 >= . {
                _apexogerr "gunits(structural) needs sigma_tilde from a successful baseline fit, but the baseline failed."
            }
            local _ginv : di %21.15g 1/scalar(_apexogSTILDE0)
        }
        if "`method'" == "uci" {
            forvalues k = 1/`count_iv' {
                local gmin`k' = `gmin`k'' * `_ginv'
                local gmax`k' = `gmax`k'' * `_ginv'
            }
        }
        if "`method'" == "ltz" & "`ltz_mode'" == "gaussian" {
            * convert the prior parameters and rebuild the draw expressions
            forvalues k = 1/`count_iv' {
                local mu`k'    = `mu`k''    * `_ginv'
                local omega`k' = `omega`k'' * (`_ginv')^2
                local _sd_k = sqrt(`omega`k'')
                local gammaCall`k' = "`mu`k'' + `_sd_k'*rnormal()"
            }
        }
        if "`method'" == "ltz" & "`ltz_mode'" == "distribution" & "`special_pending'" != "1" {
            * non-scale parameters (chi2 df, t df, ...) cannot be rescaled;
            * rescaling the DRAW is equivalent and universal
            forvalues k = 1/`count_iv' {
                local gammaCall`k' = "(`gammaCall`k'')*(`_ginv')"
            }
        }
        * distribution(special): gammaCall is rebuilt inside the LTZ block
        * (after touse exists); the wrap is applied there.
    }

    *==========================================================================
    * (7) UCI: grid over [gmin_k, gmax_k] for each instrument k
    *==========================================================================
    if "`method'" == "uci" {
        local total_pts = `grid'^`count_iv'
        * Guard the curse of dimensionality: the UCI grid is grid^(#instruments)
        * and is stored in a Stata matrix capped at c(max_matsize) rows (11000
        * in Stata/SE).  Past that the J() allocation dies with a cryptic
        * matsize error.  Catch it early with an actionable message: coarsen
        * the grid, drop instruments, or use the draw-based LTZ (no grid).
        if `total_pts' > c(max_matsize) {
            di as error "UCI grid has `total_pts' points (grid(`grid')^`count_iv' instruments), " ///
                "above the matrix limit c(max_matsize) = `c(max_matsize)'."
            di as error "This is the CHR UCI curse of dimensionality.  Options: lower grid(), " ///
                "reduce the number of instruments, or use method ltz (draws, no grid)."
            exit 908
        }
        di as text _newline "Running UCI grid: `grid' point(s) per instrument, " /*
            */ "`count_iv' instrument(s), `total_pts' total grid point(s)."

        * n_targets = how many "things" we report bounds for.
        *   APE mode: one APE per endog (n_targets = count2)
        *   AIE mode: a single AIE bound       (n_targets = 1)
        if `aie_active' local n_targets = 1
        else            local n_targets = `count2'

        * Target names for labeling output, ereturn, etc.
        if `aie_active' local target_labels "aie"
        else            local target_labels "`varlist2'"

        * Grid matrix: K gamma columns then 3 columns per target (ape/lo/hi or aie/lo/hi).
        matrix _apexogGRID = J(`total_pts', `count_iv' + 3*`n_targets', .)
        local colnames ""
        forvalues k = 1/`count_iv' {
            local colnames "`colnames' g`k'"
        }
        forvalues j = 1/`n_targets' {
            if `aie_active' local colnames "`colnames' aie_`_labX'_`_labZ' lo_`_labX'_`_labZ' hi_`_labX'_`_labZ'"
            else            local colnames "`colnames' ape_`_lab_`j'' lo_`_lab_`j'' hi_`_lab_`j''"
        }
        matrix colnames _apexogGRID = `colnames'
        * v0.20: exact Pr-scale value of each instrument's violation per row
        matrix _apexogGRIDPR = J(`total_pts', `count_iv', .)
        local colnames ""
        forvalues k = 1/`count_iv' {
            local colnames "`colnames' gpr`k'"
        }
        matrix colnames _apexogGRIDPR = `colnames'

        * Per-target union trackers
        local lo_union = .
        local hi_union = .
        forvalues j = 1/`n_targets' {
            local lo_union_`j' = .
            local hi_union_`j' = .
        }

        forvalues iter = 1/`total_pts' {
            * Decode iter into a tuple of grid indices, one per IV
            local R = `iter' - 1
            forvalues k = `count_iv'(-1)1 {
                local idx`k' = floor(`R' / (`grid'^(`k'-1)))
                local R = `R' - (`grid'^(`k'-1)) * `idx`k''
                local g`k' = `gmin`k'' + ((`gmax`k'' - `gmin`k'') / (`grid' - 1)) * `idx`k''
            }

            * Build offset variable: gen offvar = sum_k g_k * Z_k
            qui replace `offvar' = 0 if `touse'
            local k = 1
            local _gams ""
            foreach z of local varlist_iv {
                qui replace `offvar' = `offvar' + `g`k''*`z' if `touse'
                local _gams "`_gams' `g`k''"
                local ++k
            }

            _apexogfit, touse(`touse')                                            ///
                yvar(`yvar') endoglist(`varlist2') chat0list(`chat0_list')      ///
                offvar(`offvar') aweightlist(`A_list')                          ///
                exoglist(`varlist1') ivlist(`varlist_iv')                       ///
                clustvar(`clustvar') vcetype(`vcetype')                         ///
                aiex(`aie_X') aiez(`aie_Z') aiexz(`aie_XZ') aiecase(`aie_case') cfprodlist(`cfprod_list') engine(`engine')                     ///
                gammas(`_gams') ivbin(`ivbin')                                  ///
                level(`level') `debug'
            matrix _apexogITER_APE = r(ape_vec)
            matrix _apexogITER_LO  = r(lo_vec)
            matrix _apexogITER_HI  = r(hi_vec)
            capture matrix _apexogITER_GPR = r(gpr_vec)
            local _gpr_ok = (_rc == 0)

            forvalues k = 1/`count_iv' {
                matrix _apexogGRID[`iter', `k'] = `g`k''
                if `_gpr_ok' matrix _apexogGRIDPR[`iter', `k'] = _apexogITER_GPR[`k', 1]
            }
            forvalues j = 1/`n_targets' {
                local col_ape = `count_iv' + 3*(`j'-1) + 1
                local col_lo  = `count_iv' + 3*(`j'-1) + 2
                local col_hi  = `count_iv' + 3*(`j'-1) + 3
                local ape_ij  = _apexogITER_APE[`j', 1]
                local lo_ij   = _apexogITER_LO[`j', 1]
                local hi_ij   = _apexogITER_HI[`j', 1]
                matrix _apexogGRID[`iter', `col_ape'] = `ape_ij'
                matrix _apexogGRID[`iter', `col_lo']  = `lo_ij'
                matrix _apexogGRID[`iter', `col_hi']  = `hi_ij'
                local lo_union_`j' = min(`lo_union_`j'', `lo_ij')
                local hi_union_`j' = max(`hi_union_`j'', `hi_ij')
            }
            * Back-compat: lo/hi_union track FIRST target
            local lo_union = `lo_union_1'
            local hi_union = `hi_union_1'
        }

        *--- Display ---------------------------------------------------------
        if "`header'" == "" {
            di as text _newline "{hline 78}"
            if `aie_active' {
                di as text "Plausibly-exogenous AIE bounds (CHR 2012, UCI method)"
            }
            else {
                di as text "Plausibly-exogenous APE bounds (CHR 2012, UCI method)"
            }
            local _blopt ""
            local _binlin : list clean _binlin
            if "`_binlin'" != "" local _blopt "binlin(`_binlin')"
            if "`_binlin'" != "" & `aie_active' local _blopt "`_blopt' aie"
            _apexoghdr, engine(`engine') firststage(`firststage') n(`N') `_blopt'
            if `aie_active' {
                if "`target_aie_spec'" != "" {
                    di as text "Target:    Average Interaction Effect (`target_aie_spec')"
                }
                else {
                    di as text "Target:    Average Interaction Effect (`aie_X' ## `aie_Z')"
                }
            }
            else {
                di as text "Endogenous: `varlist2'"
                if "`target'" != "" {
                    di as text "Target:    APE of {`target_list'} (selected via target() option)"
                }
                else {
                    di as text "Target:    APE of every endogenous regressor (default)"
                }
            }
            di as text "{hline 78}"
            local lvl = `level' * 100
            di as text "Baseline (gamma = 0):"
            if `aie_active' {
                local b_ape_1 = _apexogBASE_APE[1, 1]
                local b_se_1  = _apexogBASE_SE[1, 1]
                local b_lo_1  = _apexogBASE_LO[1, 1]
                local b_hi_1  = _apexogBASE_HI[1, 1]
                di as text "    AIE(`aie_X' ## `aie_Z') = " %9.6f `b_ape_1' "  (SE " %7.6f `b_se_1' ")  [" %9.6f `b_lo_1' ", " %9.6f `b_hi_1' "]"
            }
            else {
                local j = 0
                foreach endogvar of local varlist2 {
                    local ++j
                    if !`target_endog_`j'' continue
                    local b_ape_j = _apexogBASE_APE[`j', 1]
                    local b_se_j  = _apexogBASE_SE[`j', 1]
                    local b_lo_j  = _apexogBASE_LO[`j', 1]
                    local b_hi_j  = _apexogBASE_HI[`j', 1]
                    di as text "    APE(`endogvar')              = " %9.6f `b_ape_j' "  (SE " %7.6f `b_se_j' ")  [" %9.6f `b_lo_j' ", " %9.6f `b_hi_j' "]"
                }
            }
            di as text " "
            di as text "Scale:     1 latent-index unit of gamma = " %6.4f scalar(_apexogMPHI0) ///
                " pr. points per unit of Z (avg., baseline)"
            di as text "           1 structural-latent unit of gamma = " ///
                %6.4f 1/scalar(_apexogSTILDE0) " latent-index units" ///
                "  (sigma_tilde = " %6.4f scalar(_apexogSTILDE0) ")"
            if "`engine'" == "biprobit" {
                di as text "           biprobit: the outcome-equation error has unit variance, so the"
                di as text "           structural-latent and index scales coincide (sigma_tilde = 1)."
                di as text "           rho(eps_y, eps_d) at the baseline fit = " %7.4f scalar(_apexogBPRHO0)
            }
            if "`gunits'" == "pr" {
                di as text "           gunits(pr): user gamma values were probability points;" ///
                    " divided by m_phi.  Values below are latent units."
                di as text "           (first order; exact Pr-scale equivalents shown alongside)"
            }
            if "`gunits'" == "structural" {
                di as text "           gunits(structural): user gamma values were structural-latent" ///
                    " units; divided by sigma_tilde.  Values below are latent-index units."
                if `_grfs' {
                    * v0.19.1: the conversion is a linear-first-stage result
                    di as text "Warning:   sigma_tilde conversion derived for a LINEAR first stage (see help)."
                }
                * v0.21.2: baseline conversion, and the limit of the structural scale
                local _binlin : list clean _binlin
                if "`_binlin'" != "" {
                    di as text "Warning:   sigma_tilde conversion not derived for a 0/1 regressor (see Note)."
                }
                else if "`engine'" == "cf" & !`_grfs' {
                    di as text "           The conversion uses the baseline sigma_tilde, exact only at gamma = 0."
                }
                if scalar(_apexogGSMAX) < . {
                    di as text "           A structural value has a unique index-scale counterpart only below"
                    di as text "           |pi_z|/sigma_v = " %6.4f scalar(_apexogGSMAX) " in absolute value (see help)."
                    if `_gs_boxmax' < . {
                        if `_gs_boxmax' >= scalar(_apexogGSMAX) {
                            di as text "Warning:   the box reaches " %6.4f `_gs_boxmax' " structural units, at or beyond that limit."
                        }
                    }
                }
            }
            di as text " "
            di as text "UCI union over gamma grid (`grid' point(s) per IV):"
            forvalues k = 1/`count_iv' {
                local zname : word `k' of `ivlabels'
                * v0.20: the same range on the exact Pr scale (over the refits)
                mata: st_local("_prlo", strofreal(min(select(st_matrix("_apexogGRIDPR")[., `k'], st_matrix("_apexogGRIDPR")[., `k'] :< .)), "%21.15g"))
                mata: st_local("_prhi", strofreal(max(select(st_matrix("_apexogGRIDPR")[., `k'], st_matrix("_apexogGRIDPR")[., `k'] :< .)), "%21.15g"))
                di as text "    gamma_`zname' range  : [" %7.4f `gmin`k'' ", " %7.4f `gmax`k'' "]" ///
                    "   exact Pr scale: [" %7.4f `_prlo' ", " %7.4f `_prhi' "]"
            }
            di as text " "
            if `aie_active' {
                di as text "    AIE(`aie_X' ## `aie_Z') bounds = [" %9.6f `lo_union_1' ", " %9.6f `hi_union_1' "]"
            }
            else {
                local j = 0
                foreach endogvar of local varlist2 {
                    local ++j
                    if !`target_endog_`j'' continue
                    di as text "    APE(`endogvar') bounds       = [" %9.6f `lo_union_`j'' ", " %9.6f `hi_union_`j'' "]"
                }
            }
            di as text "{hline 78}"
        }

        *--- Optional graph: APE/AIE vs gamma over the UCI grid ----------------
        * Do this BEFORE the `ereturn matrix' so _apexogGRID still exists.
        if "`graph'" != "" & `count_iv' == 1 {
            preserve
            clear
            svmat _apexogGRID, names(col)
            svmat _apexogGRIDPR, names(col)
            local zname : word 1 of `ivlabels'
            * v0.20: with gunits(pr) the x-axis is the exact Pr scale of the refits
            if "`gunits'" == "pr" {
                local _gx  "gpr1"
                local _gxt "gamma_`zname' (exact Pr points)"
            }
            else {
                local _gx  "g1"
                local _gxt "gamma_`zname'"
            }
            sort `_gx'
            if `aie_active' {
                * (v0.18 fix: the svmat columns are aie_*, not ape_*)
                local _ape_lab "aie_`_labX'_`_labZ'"
                local _lo_lab  "lo_`_labX'_`_labZ'"
                local _hi_lab  "hi_`_labX'_`_labZ'"
                local _title   "UCI AIE bounds vs assumed gamma_`zname'"
                local _ytitle  "AIE (probability scale)"
            }
            else {
                local _ape_lab "ape_`_lab_1'"
                local _lo_lab  "lo_`_lab_1'"
                local _hi_lab  "hi_`_lab_1'"
                local _title   "UCI APE bounds vs assumed gamma_`zname'"
                local _ytitle  "APE (probability scale)"
            }
            local _lvlg = `level'*100
            if `"`graphoptions'"' == "" {
                local graphoptions ///
                    title("`_title'") ///
                    ytitle("`_ytitle'")                               ///
                    xtitle("`_gxt'")                                   ///
                    legend(order(1 "Lower `_lvlg'% CI" 2 "Upper `_lvlg'% CI" 3 "Point estimate"))
            }
            twoway (line `_lo_lab' `_gx', lpattern(dash) lcolor(black))   ///
                   (line `_hi_lab' `_gx', lpattern(dash) lcolor(black))   ///
                   (line `_ape_lab' `_gx', lcolor(blue)),                 ///
                   yline(0, lpattern(dot) lcolor(gray))                  ///
                   `graphoptions'
            restore
        }
        else if "`graph'" != "" & `count_iv' > 1 {
            di as text "  [note] graph is only supported with a single instrument; skipping."
        }

        * Post the baseline target estimates as e(b)/e(V) with e(sample), so
        * `if e(sample)', test/lincom, estimates store, and replay all work.
        * e(b) = baseline (gamma=0) APE per endog (or the 1x1 AIE); e(V) = the
        * Murphy-Topel delta-method covariance of those targets (cross-target
        * terms from the shared V_beta).  Falls back to a bare ereturn clear
        * if the baseline fit failed.
        local _postok = !matmissing(_apexogBASE_APE) & !matmissing(_apexogVTGT0)
        tempname _bT _VT
        if `_postok' {
            matrix `_bT' = _apexogBASE_APE'
            matrix `_VT' = _apexogVTGT0
            local _tgtnames ""
            if `aie_active' local _tgtnames "aie"
            else {
                local j = 0
                foreach endogvar of local varlist2 {
                    local ++j
                    if strpos("`endogvar'", "__") == 1 local _tgtnames "`_tgtnames' endog`j'"
                    else                               local _tgtnames "`_tgtnames' `endogvar'"
                }
            }
            matrix colnames `_bT' = `_tgtnames'
            matrix colnames `_VT' = `_tgtnames'
            matrix rownames `_VT' = `_tgtnames'
            ereturn post `_bT' `_VT', esample(`touse') obs(`N')
        }
        else {
            ereturn clear
            ereturn scalar N         = `N'
        }
        ereturn scalar level         = `level'
        ereturn scalar mphi          = scalar(_apexogMPHI0)
        ereturn scalar stilde        = scalar(_apexogSTILDE0)
        ereturn scalar gs_max        = scalar(_apexogGSMAX)
        if "`engine'" == "biprobit" ereturn scalar rho = scalar(_apexogBPRHO0)
        ereturn local  gunits        "`gunits'"
        ereturn scalar baseline_ape  = `baseline_ape'
        ereturn scalar baseline_se   = `baseline_se'
        ereturn scalar baseline_lo   = `baseline_lo'
        ereturn scalar baseline_hi   = `baseline_hi'
        if `aie_active' {
            * AIE mode: single bound stored under e(lb_aie) / e(ub_aie)
            ereturn scalar lb_aie = `lo_union_1'
            ereturn scalar ub_aie = `hi_union_1'
        }
        else {
            * Back-compat: first endog scalars.  Only set when target() was NOT
            * used (Phase 1: respect target() exclusions; legacy callers that
            * don't use target() get the back-compat scalars as before).
            if "`target'" == "" {
                ereturn scalar lb_`_lab_1'   = `lo_union'
                ereturn scalar ub_`_lab_1'   = `hi_union'
            }
            * Per-endog scalars -- only for targets (Phase 1: skip non-targets)
            local j = 0
            foreach endogvar of local varlist2 {
                local ++j
                if !`target_endog_`j'' continue
                ereturn scalar lb_`_lab_`j'' = `lo_union_`j''
                ereturn scalar ub_`_lab_`j'' = `hi_union_`j''
            }
        }
        ereturn matrix grid               = _apexogGRID
        ereturn matrix grid_pr            = _apexogGRIDPR
        ereturn matrix baseline_ape_vec   = _apexogBASE_APE
        ereturn matrix baseline_se_vec    = _apexogBASE_SE
        ereturn matrix baseline_lo_vec    = _apexogBASE_LO
        ereturn matrix baseline_hi_vec    = _apexogBASE_HI
        if `aie_active' {
            ereturn local target          "aie"
            ereturn local aie_vars        "`aie_X' `aie_Z'"
            if "`target_aie_spec'" != "" ereturn local aie_spec "`target_aie_spec'"
        }
        else {
            ereturn local target          "ape"
            ereturn local target_list     "`target_list'"
        }
        ereturn local cmd            "apexog"
        ereturn local cmdline        `"apexog `_cmdline_full'"'
        ereturn local method         "uci"
        ereturn local estimator      "`_est_name'"
        ereturn local firststage     "`firststage'"
        ereturn local vce            "`_vce_name'"
        ereturn local depvar         "`yvar'"
        ereturn local endog          "`varlist2'"
        ereturn local exog           "`varlist1'"
        ereturn local instruments    "`ivlabels'"
        if "`vcetype'" == "cluster" {
            ereturn scalar N_clust = scalar(_apexogNCLUST)
            ereturn local  clustvar "`clustvar_name'"
        }
    }

    *==========================================================================
    * (8) LTZ: simulate gamma ~ N(mu, diag(omega)), refit per draw
    *==========================================================================
    if "`method'" == "ltz" {
        di as text _newline "Running LTZ: `iterations' simulation draw(s) " /*
            */ "from `ltz_mode' prior over `count_iv' instrument violation(s)."

        * AIE mode: single target.  APE mode: per-endog.
        if `aie_active' local n_targets = 1
        else            local n_targets = `count2'

        matrix _apexogAPESIM = J(`iterations', `n_targets', .)
        matrix _apexogSESIM  = J(`iterations', `n_targets', .)
        local sim_colnames ""
        if `aie_active' {
            local sim_colnames "aie_`_labX'_`_labZ'"
        }
        else {
            forvalues j = 1/`count2' {
                local sim_colnames "`sim_colnames' ape_`_lab_`j''"
            }
        }
        matrix colnames _apexogAPESIM = `sim_colnames'
        matrix colnames _apexogSESIM  = `sim_colnames'
        * v0.20: the latent draws and their exact Pr-scale values
        matrix _apexogGSIM   = J(`iterations', `count_iv', .)
        matrix _apexogGPRSIM = J(`iterations', `count_iv', .)
        local _gcn ""
        local _gprcn ""
        forvalues k = 1/`count_iv' {
            local _gcn   "`_gcn' g`k'"
            local _gprcn "`_gprcn' gpr`k'"
        }
        matrix colnames _apexogGSIM   = `_gcn'
        matrix colnames _apexogGPRSIM = `_gprcn'

        * Build the donor matrix for distribution(special) now that touse exists.
        if "`special_pending'" == "1" {
            tempname _apexogDONOR
            qui mkmat `special_donor_vars' if `touse', matrix(`_apexogDONOR') nomissing
            local _apexog_nrows = rowsof(`_apexogDONOR')
            forvalues k = 1/`count_iv' {
                local gammaCall`k' = "`_apexogDONOR'[ceil(runiform()*`_apexog_nrows'), `k']"
                if "`gunits'" == "pr" {
                    local _ginv : di %21.15g 1/scalar(_apexogMPHI0)
                    local gammaCall`k' = "(`gammaCall`k'')*(`_ginv')"
                }
                else if "`gunits'" == "structural" {
                    local _ginv : di %21.15g 1/scalar(_apexogSTILDE0)
                    local gammaCall`k' = "(`gammaCall`k'')*(`_ginv')"
                }
            }
        }

        forvalues iter = 1/`iterations' {
            * Draw gamma_k from the per-IV gammaCall expression
            * (gammaCall`k' was built above for both Gaussian and non-Gaussian
            * priors; for Gaussian it's "mu + sd*rnormal()", for distribution()
            * it's the appropriate Stata RNG expression).
            qui replace `offvar' = 0 if `touse'
            local k = 1
            local _gams ""
            foreach z of local varlist_iv {
                local g_draw = `gammaCall`k''
                qui replace `offvar' = `offvar' + `g_draw'*`z' if `touse'
                matrix _apexogGSIM[`iter', `k'] = `g_draw'
                local _gams "`_gams' `g_draw'"
                local ++k
            }

            _apexogfit, touse(`touse')                                            ///
                yvar(`yvar') endoglist(`varlist2') chat0list(`chat0_list')      ///
                offvar(`offvar') aweightlist(`A_list')                          ///
                exoglist(`varlist1') ivlist(`varlist_iv')                       ///
                clustvar(`clustvar') vcetype(`vcetype')                         ///
                aiex(`aie_X') aiez(`aie_Z') aiexz(`aie_XZ') aiecase(`aie_case') cfprodlist(`cfprod_list') engine(`engine')                     ///
                gammas(`_gams') ivbin(`ivbin')                                  ///
                level(`level') `debug'

            matrix _apexogITER_APE = r(ape_vec)
            matrix _apexogITER_SE  = r(se_vec)
            forvalues j = 1/`n_targets' {
                matrix _apexogAPESIM[`iter', `j'] = _apexogITER_APE[`j', 1]
                matrix _apexogSESIM[`iter', `j']  = _apexogITER_SE[`j', 1]
            }
            capture matrix _apexogITER_GPR = r(gpr_vec)
            if _rc == 0 {
                forvalues k = 1/`count_iv' {
                    matrix _apexogGPRSIM[`iter', `k'] = _apexogITER_GPR[`k', 1]
                }
            }

            if mod(`iter', 100) == 0 {
                di as text "  ... `iter' / `iterations' draws done"
            }
        }

        mata: st_local("_n_missing", strofreal(sum(rowmissing(st_matrix("_apexogAPESIM")) :| rowmissing(st_matrix("_apexogSESIM")))))
        if `_n_missing' > 0 {
            di as text "  [note] `_n_missing' of `iterations' LTZ draws failed (probit did not converge or MT SE missing); excluded from percentile interval."
        }

        mata: _apexogsim_pct(`level')
        matrix colnames _apexogSIMCNV = `sim_colnames'

        if "`header'" == "" {
            di as text _newline "{hline 78}"
            if `aie_active' {
                di as text "Plausibly-exogenous AIE bounds (CHR 2012, LTZ method, `ltz_mode' prior)"
            }
            else {
                di as text "Plausibly-exogenous APE bounds (CHR 2012, LTZ method, `ltz_mode' prior)"
            }
            local _blopt ""
            local _binlin : list clean _binlin
            if "`_binlin'" != "" local _blopt "binlin(`_binlin')"
            if "`_binlin'" != "" & `aie_active' local _blopt "`_blopt' aie"
            _apexoghdr, engine(`engine') firststage(`firststage') n(`N') `_blopt'
            if `aie_active' {
                if "`target_aie_spec'" != "" {
                    di as text "Target:    Average Interaction Effect (`target_aie_spec')"
                }
                else {
                    di as text "Target:    Average Interaction Effect (`aie_X' ## `aie_Z')"
                }
            }
            else {
                di as text "Endogenous: `varlist2'"
                if "`target'" != "" {
                    di as text "Target:    APE of {`target_list'} (selected via target() option)"
                }
                else {
                    di as text "Target:    APE of every endogenous regressor (default)"
                }
            }
            di as text "{hline 78}"
            di as text "Baseline (gamma = 0):"
            if `aie_active' {
                local b_ape_1 = _apexogBASE_APE[1, 1]
                local b_se_1  = _apexogBASE_SE[1, 1]
                di as text "    AIE(`aie_X' ## `aie_Z') = " %9.6f `b_ape_1' "  (SE " %7.6f `b_se_1' ")"
            }
            else {
                local j = 0
                foreach endogvar of local varlist2 {
                    local ++j
                    if !`target_endog_`j'' continue
                    local b_ape_j = _apexogBASE_APE[`j', 1]
                    local b_se_j  = _apexogBASE_SE[`j', 1]
                    di as text "    APE(`endogvar')              = " %9.6f `b_ape_j' "  (SE " %7.6f `b_se_j' ")"
                }
            }
            di as text " "
            di as text "Scale:     1 latent-index unit of gamma = " %6.4f scalar(_apexogMPHI0) ///
                " pr. points per unit of Z (avg., baseline)"
            di as text "           1 structural-latent unit of gamma = " ///
                %6.4f 1/scalar(_apexogSTILDE0) " latent-index units" ///
                "  (sigma_tilde = " %6.4f scalar(_apexogSTILDE0) ")"
            if "`engine'" == "biprobit" {
                di as text "           biprobit: the outcome-equation error has unit variance, so the"
                di as text "           structural-latent and index scales coincide (sigma_tilde = 1)."
                di as text "           rho(eps_y, eps_d) at the baseline fit = " %7.4f scalar(_apexogBPRHO0)
            }
            if "`gunits'" == "pr" {
                di as text "           gunits(pr): user gamma values were probability points;" ///
                    " divided by m_phi.  Values below are latent units."
                di as text "           (first order; exact Pr-scale equivalents shown alongside)"
            }
            if "`gunits'" == "structural" {
                di as text "           gunits(structural): user gamma values were structural-latent" ///
                    " units; divided by sigma_tilde.  Values below are latent-index units."
                if `_grfs' {
                    * v0.19.1: the conversion is a linear-first-stage result
                    di as text "Warning:   sigma_tilde conversion derived for a LINEAR first stage (see help)."
                }
                * v0.21.2: baseline conversion, and the limit of the structural scale
                local _binlin : list clean _binlin
                if "`_binlin'" != "" {
                    di as text "Warning:   sigma_tilde conversion not derived for a 0/1 regressor (see Note)."
                }
                else if "`engine'" == "cf" & !`_grfs' {
                    di as text "           The conversion uses the baseline sigma_tilde, exact only at gamma = 0."
                }
                if scalar(_apexogGSMAX) < . {
                    di as text "           A structural value has a unique index-scale counterpart only below"
                    di as text "           |pi_z|/sigma_v = " %6.4f scalar(_apexogGSMAX) " in absolute value (see help)."
                    * share of the prior draws (index units * sigma_tilde = structural units) at or beyond it
                    mata: st_local("_gs_share", strofreal(mean(abs(st_matrix("_apexogGSIM")[., 1]) :* st_numscalar("_apexogSTILDE0") :>= st_numscalar("_apexogGSMAX")), "%21.15g"))
                    if `_gs_share' > 0 {
                        di as text "Warning:   " %4.1f 100*`_gs_share' "% of the prior draws are at or beyond that limit."
                    }
                }
            }
            di as text " "
            di as text "LTZ prior:"
            forvalues k = 1/`count_iv' {
                local zname : word `k' of `ivlabels'
                if "`ltz_mode'" == "gaussian" {
                    di as text "    gamma_`zname' ~ N(" %7.4f `mu`k'' ", " %7.4f `omega`k'' ")"
                }
                else {
                    di as text "    gamma_`zname' ~ `distname'   draw expr: `gammaCall`k''"
                }
                * v0.20: the draws on the exact Pr scale (over the refits)
                mata: _apexog_v = select(st_matrix("_apexogGPRSIM")[., `k'], st_matrix("_apexogGPRSIM")[., `k'] :< .)
                mata: st_local("_prm", strofreal(mean(_apexog_v), "%21.15g"))
                mata: st_local("_prs", strofreal(sqrt(variance(_apexog_v)), "%21.15g"))
                mata: mata drop _apexog_v
                di as text "      exact Pr scale of the draws: mean " %7.4f `_prm' ", sd " %7.4f `_prs'
            }
            di as text "    iterations         = `iterations'"
            di as text " "
            local lvl = `level' * 100
            if `aie_active' {
                di as text "LTZ `lvl'% interval (sampling noise convolved with prior; CHR 2012 eq. 7):"
                scalar _apexogLO = _apexogLO_PCT[1, 1]
                scalar _apexogHI = _apexogHI_PCT[1, 1]
                di as text "    AIE(`aie_X' ## `aie_Z') `lvl'% CI = [" %9.6f _apexogLO ", " %9.6f _apexogHI "]"
                di as text " "
                di as text "Prior-only spread of the point AIEs (sensitivity band, no sampling noise):"
                scalar _apexogLO = _apexogLO_PRIOR[1, 1]
                scalar _apexogHI = _apexogHI_PRIOR[1, 1]
                di as text "    AIE(`aie_X' ## `aie_Z') `lvl'% range = [" %9.6f _apexogLO ", " %9.6f _apexogHI "]"
            }
            else {
                di as text "LTZ `lvl'% intervals (sampling noise convolved with prior; CHR 2012 eq. 7):"
                local j = 0
                foreach endogvar of local varlist2 {
                    local ++j
                    if !`target_endog_`j'' continue
                    scalar _apexogLO = _apexogLO_PCT[`j', 1]
                    scalar _apexogHI = _apexogHI_PCT[`j', 1]
                    di as text "    APE(`endogvar') `lvl'% CI = [" %9.6f _apexogLO ", " %9.6f _apexogHI "]"
                }
                di as text " "
                di as text "Prior-only spread of the point APEs (sensitivity band, no sampling noise):"
                local j = 0
                foreach endogvar of local varlist2 {
                    local ++j
                    if !`target_endog_`j'' continue
                    scalar _apexogLO = _apexogLO_PRIOR[`j', 1]
                    scalar _apexogHI = _apexogHI_PRIOR[`j', 1]
                    di as text "    APE(`endogvar') `lvl'% range = [" %9.6f _apexogLO ", " %9.6f _apexogHI "]"
                }
            }
            di as text "{hline 78}"
        }

        *--- Optional delta sweep: rerun the LTZ simulation for a sequence of
        *    prior widths (graphomega(); one value per step, broadcast to all
        *    instruments; Gaussian priors only) -- the CHR "bounds vs delta"
        *    figure.  Runs AFTER the main results, so the headline above is
        *    unchanged; uses dedicated _apexogSW* matrices so nothing the
        *    main path stored is clobbered.  Tracks the FIRST target (the
        *    AIE, or the first endog's APE).
        local n_sweep : word count `graphomega'
        if `n_sweep' > 0 {
            if "`ltz_mode'" != "gaussian" {
                _apexogerr "graphomega() requires the Gaussian mu()/omega() prior (not distribution())."
            }
            local n_smu : word count `graphmu'
            if `n_smu' == 0 {
                forvalues s = 1/`n_sweep' {
                    local graphmu "`graphmu' 0"
                }
            }
            else if `n_smu' != `n_sweep' {
                _apexogerr "graphmu() must have as many values as graphomega()."
            }
            local n_sd : word count `graphdelta'
            if `n_sd' > 0 & `n_sd' != `n_sweep' {
                _apexogerr "graphdelta() must have as many values as graphomega()."
            }
            matrix _apexogSWEEP = J(`n_sweep', 6, .)
            matrix colnames _apexogSWEEP = delta mu lb ub lb_prior ub_prior
            di as text _newline "Delta sweep: `n_sweep' prior width(s) x `iterations' draws each (first target)."
            forvalues s = 1/`n_sweep' {
                local om_s : word `s' of `graphomega'
                local mu_s : word `s' of `graphmu'
                if `n_sd' > 0 local de_s : word `s' of `graphdelta'
                else          local de_s = `om_s'
                if `om_s' <= 0 _apexogerr "graphomega() values must be > 0."
                * same unit conversion as the main prior
                local om_use = `om_s'
                local mu_use = `mu_s'
                if "`gunits'" == "pr" {
                    local _ginvs : di %21.15g 1/scalar(_apexogMPHI0)
                    local om_use = `om_s' * (`_ginvs')^2
                    local mu_use = `mu_s' * `_ginvs'
                }
                else if "`gunits'" == "structural" {
                    local _ginvs : di %21.15g 1/scalar(_apexogSTILDE0)
                    local om_use = `om_s' * (`_ginvs')^2
                    local mu_use = `mu_s' * `_ginvs'
                }
                local sd_s = sqrt(`om_use')
                matrix _apexogSWAPE = J(`iterations', 1, .)
                matrix _apexogSWSE  = J(`iterations', 1, .)
                forvalues iter = 1/`iterations' {
                    qui replace `offvar' = 0 if `touse'
                    foreach z of local varlist_iv {
                        local g_draw = `mu_use' + `sd_s'*rnormal()
                        qui replace `offvar' = `offvar' + `g_draw'*`z' if `touse'
                    }
                    _apexogfit, touse(`touse')                                            ///
                        yvar(`yvar') endoglist(`varlist2') chat0list(`chat0_list')      ///
                        offvar(`offvar') aweightlist(`A_list')                          ///
                        exoglist(`varlist1') ivlist(`varlist_iv')                       ///
                        clustvar(`clustvar') vcetype(`vcetype')                         ///
                        aiex(`aie_X') aiez(`aie_Z') aiexz(`aie_XZ') aiecase(`aie_case') cfprodlist(`cfprod_list') engine(`engine') ///
                        level(`level')
                    matrix _apexogITER_APE = r(ape_vec)
                    matrix _apexogITER_SE  = r(se_vec)
                    matrix _apexogSWAPE[`iter', 1] = _apexogITER_APE[1, 1]
                    matrix _apexogSWSE[`iter', 1]  = _apexogITER_SE[1, 1]
                }
                mata: _apexogsweep_pct(`level')
                matrix _apexogSWEEP[`s', 1] = `de_s'
                matrix _apexogSWEEP[`s', 2] = `mu_s'
                matrix _apexogSWEEP[`s', 3] = scalar(_apexogSW_lo)
                matrix _apexogSWEEP[`s', 4] = scalar(_apexogSW_hi)
                matrix _apexogSWEEP[`s', 5] = scalar(_apexogSW_lop)
                matrix _apexogSWEEP[`s', 6] = scalar(_apexogSW_hip)
                di as text "  sweep `s'/`n_sweep': omega = " %-9.0g `om_s' ///
                    "  CI [" %9.6f scalar(_apexogSW_lo) ", " %9.6f scalar(_apexogSW_hi) "]"
            }
        }

        *--- Optional graph.  With graphomega(): the CHR "bounds vs delta"
        *    sweep plot.  Otherwise: histogram of the CONVOLVED LTZ draws
        *    (APE_i + SE_i*eps_i -- the draws behind the headline CI).
        *    For multi-endog APE, both track the FIRST endog only.
        if "`graph'" != "" & `n_sweep' > 0 {
            preserve
            clear
            svmat _apexogSWEEP, names(col)
            local _xunits ""
            if "`gunits'" == "pr"         local _xunits " (pr units)"
            if "`gunits'" == "structural" local _xunits " (structural units)"
            if `"`graphoptions'"' == "" {
                local graphoptions ///
                    title("LTZ bounds vs prior width")                          ///
                    subtitle("`iterations' draws per width; first target")      ///
                    xtitle("delta = prior variance`_xunits'")                   ///
                    ytitle("APE / AIE (probability scale)")                     ///
                    legend(order(1 "LTZ CI" 3 "prior-only band"))               ///
                    note("Solid: convolved LTZ CI (CHR 2012 eq. 7).  Dashed: prior-only sensitivity band." ///
                         "Dotted horizontal line: baseline (gamma = 0) estimate.")
            }
            sort delta
            twoway (line lb delta, lcolor(blue))                        ///
                   (line ub delta, lcolor(blue))                        ///
                   (line lb_prior delta, lpattern(dash) lcolor(gray))   ///
                   (line ub_prior delta, lpattern(dash) lcolor(gray)),  ///
                   yline(`baseline_ape', lpattern(dot) lcolor(black))   ///
                   `graphoptions'
            restore
        }
        else if "`graph'" != "" {
            preserve
            clear
            svmat _apexogSIMCNV, names(col)
            if `aie_active' {
                local _histvar "aie_`_labX'_`_labZ'"
                local _title   "LTZ AIE distribution across `iterations' draws"
                local _subt    "Interaction: `aie_X' ## `aie_Z'"
                local _xt      "AIE (probability scale)"
            }
            else {
                local zfirst : word 1 of `varlist2'
                local _histvar "ape_`_lab_1'"
                local _title   "LTZ APE distribution across `iterations' draws"
                local _subt    "First endog: `zfirst'"
                local _xt      "APE (probability scale)"
            }
            if `"`graphoptions'"' == "" {
                local graphoptions ///
                    title("`_title'") ///
                    subtitle("`_subt'")                                    ///
                    xtitle("`_xt'")                                         ///
                    ytitle("Density")                                       ///
                    note("Draws convolve sampling noise with the prior (CHR 2012 eq. 7)." ///
                         "Dashed lines: `=`level'*100'% percentile bounds; solid: baseline (gamma=0)")
            }
            twoway (histogram `_histvar', frequency)                   ///
                   ,                                                   ///
                   xline(`lo_pct',       lpattern(dash) lcolor(red))   ///
                   xline(`hi_pct',       lpattern(dash) lcolor(red))   ///
                   xline(`baseline_ape', lpattern(solid) lcolor(blue)) ///
                   `graphoptions'
            restore
        }

        * Post the baseline target estimates as e(b)/e(V) with e(sample) (see
        * the UCI block for the rationale); falls back to ereturn clear if the
        * baseline fit failed.
        local _postok = !matmissing(_apexogBASE_APE) & !matmissing(_apexogVTGT0)
        tempname _bT _VT
        if `_postok' {
            matrix `_bT' = _apexogBASE_APE'
            matrix `_VT' = _apexogVTGT0
            local _tgtnames ""
            if `aie_active' local _tgtnames "aie"
            else {
                local j = 0
                foreach endogvar of local varlist2 {
                    local ++j
                    if strpos("`endogvar'", "__") == 1 local _tgtnames "`_tgtnames' endog`j'"
                    else                               local _tgtnames "`_tgtnames' `endogvar'"
                }
            }
            matrix colnames `_bT' = `_tgtnames'
            matrix colnames `_VT' = `_tgtnames'
            matrix rownames `_VT' = `_tgtnames'
            ereturn post `_bT' `_VT', esample(`touse') obs(`N')
        }
        else {
            ereturn clear
            ereturn scalar N         = `N'
        }
        ereturn scalar level         = `level'
        ereturn scalar iterations    = `iterations'
        ereturn scalar baseline_ape  = `baseline_ape'
        ereturn scalar baseline_se   = `baseline_se'
        ereturn scalar mphi          = scalar(_apexogMPHI0)
        ereturn scalar stilde        = scalar(_apexogSTILDE0)
        ereturn scalar gs_max        = scalar(_apexogGSMAX)
        if "`engine'" == "biprobit" ereturn scalar rho = scalar(_apexogBPRHO0)
        ereturn local  gunits        "`gunits'"
        if `n_sweep' > 0 {
            ereturn matrix sweep = _apexogSWEEP
        }
        if `aie_active' {
            * AIE mode: single interval stored under e(lb_aie) / e(ub_aie)
            scalar _apexogLO_temp = _apexogLO_PCT[1, 1]
            scalar _apexogHI_temp = _apexogHI_PCT[1, 1]
            ereturn scalar lb_aie = _apexogLO_temp
            ereturn scalar ub_aie = _apexogHI_temp
        }
        else {
            * Back-compat: first endog scalars.  Only set when target() was NOT
            * used (Phase 1: respect target() exclusions).
            if "`target'" == "" {
                ereturn scalar lb_`_lab_1'   = `lo_pct'
                ereturn scalar ub_`_lab_1'   = `hi_pct'
            }
            * Per-endog scalars -- only for targets (Phase 1: skip non-targets)
            local j = 0
            foreach endogvar of local varlist2 {
                local ++j
                if !`target_endog_`j'' continue
                scalar _apexogLO_temp = _apexogLO_PCT[`j', 1]
                scalar _apexogHI_temp = _apexogHI_PCT[`j', 1]
                ereturn scalar lb_`_lab_`j'' = _apexogLO_temp
                ereturn scalar ub_`_lab_`j'' = _apexogHI_temp
            }
        }
        ereturn matrix ape_sim       = _apexogAPESIM
        ereturn matrix se_sim        = _apexogSESIM
        ereturn matrix gamma_sim     = _apexogGSIM
        ereturn matrix gamma_pr_sim  = _apexogGPRSIM
        ereturn matrix ape_sim_cnv   = _apexogSIMCNV
        ereturn matrix lo_pct_vec    = _apexogLO_PCT
        ereturn matrix hi_pct_vec    = _apexogHI_PCT
        ereturn matrix lo_prior_vec  = _apexogLO_PRIOR
        ereturn matrix hi_prior_vec  = _apexogHI_PRIOR
        ereturn matrix baseline_ape_vec = _apexogBASE_APE
        ereturn matrix baseline_se_vec  = _apexogBASE_SE
        if `aie_active' {
            ereturn local target      "aie"
            ereturn local aie_vars    "`aie_X' `aie_Z'"
        }
        else {
            ereturn local target      "ape"
            ereturn local target_list "`target_list'"
        }
        ereturn local cmd            "apexog"
        ereturn local cmdline        `"apexog `_cmdline_full'"'
        ereturn local method         "ltz"
        ereturn local estimator      "`_est_name'"
        ereturn local firststage     "`firststage'"
        ereturn local vce            "`_vce_name'"
        ereturn local depvar         "`yvar'"
        ereturn local endog          "`varlist2'"
        ereturn local exog           "`varlist1'"
        ereturn local instruments    "`ivlabels'"
        if "`vcetype'" == "cluster" {
            ereturn scalar N_clust = scalar(_apexogNCLUST)
            ereturn local  clustvar "`clustvar_name'"
        }
    }

    *==========================================================================
    * (9) Housekeeping: drop the internal _apexog* matrices and scalars left in
    *     the session (all results live in e()); user objects are untouched.
    *==========================================================================
    _apexogclean
end

*==============================================================================
* _apexogfit: workhorse.  At a given gamma-offset (already baked into `offvar'):
*   (v0.19) engine(biprobit): fit `biprobit (y = X d, offset(offvar)) (d = Z X)'
*   and compute the ASF discrete-change APE + MLE delta-method SE in Mata
*   (_apexogfit_bp_ape); same r() contract as below.  Otherwise:
*   1. Fit `probit y X d chat0, offset(offvar) vce(robust)'
*   2. Pass results to Mata: compute Murphy-Topel-corrected V_beta, then
*      APE delta-method point estimate and SE.
*   3. Return ape, se, lo, hi via r().
*==============================================================================
program _apexogfit, rclass
    syntax,                                                                 ///
        touse(string)                                                       ///
        yvar(string) endoglist(string) offvar(string)                       ///
        exoglist(string) ivlist(string)                                     ///
        vcetype(string)                                                     ///
        level(real)                                                         ///
        [chat0list(string) aweightlist(string)                              ///
         clustvar(string) aiex(string) aiez(string) aiexz(string)           ///
         aiecase(integer 0) cfprodlist(string) engine(string)               ///
         gammas(string) ivbin(string) debug]

    if "`debug'" != "" {
        qui sum `offvar' if `touse'
        di as text "[debug _apexogfit] offvar `offvar' mean=" %9.6f r(mean) " sd=" %9.6f r(sd)
        di as text "[debug _apexogfit] endoglist = `endoglist'   chat0list = `chat0list'"
    }

    * Build the vce string from vcetype + clustvar
    if "`vcetype'" == "cluster" local _vce = "vce(cluster `clustvar')"
    else                         local _vce = "vce(robust)"

    *--- v0.19: recursive bivariate probit engine ------------------------------
    * firststage(biprobit): fit the outcome and the 0/1 treatment JOINTLY,
    *     biprobit (y = exog d, offset(off)) (d = iv exog)
    * with the gamma-offset in the outcome equation.  The APE (ASF discrete
    * change) and its delta-method SE come from Mata (_apexogfit_bp_ape); there
    * is no control function and no Murphy-Topel step.  Fits after the baseline
    * are warm-started from the baseline e(b) (_apexogBP_B0: created by the
    * first successful fit of the run, dropped by _apexogclean at exit).
    if "`engine'" == "biprobit" {
        local _n_endog : word count `endoglist'
        local _from ""
        capture confirm matrix _apexogBP_B0
        if _rc == 0 local _from "from(_apexogBP_B0, copy)"
        capture {
            if "`debug'" != "" {
                biprobit (`yvar' = `exoglist' `endoglist', offset(`offvar'))      ///
                         (`endoglist' = `ivlist' `exoglist') if `touse', `_vce' `_from'
            }
            else {
                quietly biprobit (`yvar' = `exoglist' `endoglist', offset(`offvar')) ///
                         (`endoglist' = `ivlist' `exoglist') if `touse', `_vce' `_from'
            }
        }
        local _rc = _rc
        local _bad ""
        if `_rc' != 0 local _bad "biprobit failed (rc=`_rc')"
        else {
            qui count if `touse'
            if e(converged) == 0 local _bad "biprobit did not converge"
            else if e(N) < r(N)  local _bad "biprobit used `e(N)' of `r(N)' obs"
            else if [#1]_b[`endoglist'] == 0 & [#1]_se[`endoglist'] == 0 {
                local _bad "`endoglist' omitted from the outcome equation (collinear)"
            }
        }
        if "`_bad'" != "" {
            qui sum `offvar' if `touse'
            di as text "  [warn] `_bad' for offvar mean=" %9.6f r(mean) "; returning missing for this fit."
        matrix _apexogAPE_vec = J(`_n_endog', 1, .)
        matrix _apexogSE_vec  = J(`_n_endog', 1, .)
        matrix _apexogLO_vec  = J(`_n_endog', 1, .)
        matrix _apexogHI_vec  = J(`_n_endog', 1, .)
        return matrix ape_vec = _apexogAPE_vec
        return matrix se_vec  = _apexogSE_vec
        return matrix lo_vec  = _apexogLO_vec
        return matrix hi_vec  = _apexogHI_vec
        return scalar ape     = .
        return scalar se      = .
        return scalar lo      = .
        return scalar hi      = .
        return scalar mphi    = .
        return scalar stilde  = .
        matrix _apexogVTGT    = J(`_n_endog', `_n_endog', .)
        return matrix vtgt    = _apexogVTGT
        exit
        }
        matrix _apexogPR_B = e(b)
        matrix _apexogPR_V = e(V)
        if "`_from'" == "" matrix _apexogBP_B0 = e(b)
        scalar _apexogBPRHO = e(rho)
        capture noisily ///
        mata: _apexogfit_bp_ape("`endoglist'", "`offvar'", "`exoglist'",          ///
                                "_apexogAPE_vec", "_apexogSE_vec",                  ///
                                "_apexogLO_vec",  "_apexogHI_vec", `level')
        if _rc != 0 {
            di as error "[_apexogfit] Mata biprobit APE failed (rc=`_rc')."
        matrix _apexogAPE_vec = J(`_n_endog', 1, .)
        matrix _apexogSE_vec  = J(`_n_endog', 1, .)
        matrix _apexogLO_vec  = J(`_n_endog', 1, .)
        matrix _apexogHI_vec  = J(`_n_endog', 1, .)
        return matrix ape_vec = _apexogAPE_vec
        return matrix se_vec  = _apexogSE_vec
        return matrix lo_vec  = _apexogLO_vec
        return matrix hi_vec  = _apexogHI_vec
        return scalar ape     = .
        return scalar se      = .
        return scalar lo      = .
        return scalar hi      = .
        return scalar mphi    = .
        return scalar stilde  = .
        matrix _apexogVTGT    = J(`_n_endog', `_n_endog', .)
        return matrix vtgt    = _apexogVTGT
        exit
        }
        * v0.20: exact Pr-scale value of each instrument's violation at this
        * refit (before any `return': _apexoggpr is r-class).
        if "`gammas'" != "" {
            _apexoggpr, touse(`touse') ivlist(`ivlist') gammas(`gammas') ivbin(`ivbin') eq(xb1)
            matrix _apexogGPR_vec = r(gpr)
        }
        return scalar ape    = _apexogAPE_vec[1, 1]
        return scalar se     = _apexogSE_vec[1, 1]
        return scalar lo     = _apexogLO_vec[1, 1]
        return scalar hi     = _apexogHI_vec[1, 1]
        return scalar mphi   = scalar(_apexogMPHI)
        return scalar stilde = scalar(_apexogSTILDE)
        if "`gammas'" != "" return matrix gpr_vec = _apexogGPR_vec
        if "`debug'" != "" di as text "[debug _apexogfit] biprobit APE=" %12.7f _apexogAPE_vec[1,1] " SE=" %12.7f _apexogSE_vec[1,1] " rho=" %8.4f scalar(_apexogBPRHO)
        return matrix vtgt    = _apexogVTGT
        return matrix ape_vec = _apexogAPE_vec
        return matrix se_vec  = _apexogSE_vec
        return matrix lo_vec  = _apexogLO_vec
        return matrix hi_vec  = _apexogHI_vec
        exit
    }

    * Run probit with the offset.  Use `capture noisily quietly' pattern.
    capture {
        if "`debug'" != "" {
            probit `yvar' `exoglist' `endoglist' `chat0list' `cfprodlist' if `touse', ///
                offset(`offvar') `_vce'
        }
        else {
            quietly probit `yvar' `exoglist' `endoglist' `chat0list' `cfprodlist' if `touse', ///
                offset(`offvar') `_vce'
        }
    }
    local _rc = _rc

    * Number of endogs (for empty-fallback matrices on failure)
    local _n_endog : word count `endoglist'

    if `_rc' != 0 {
        qui sum `offvar' if `touse'
        di as text "  [warn] probit failed (rc=`_rc') for offvar mean=" %9.6f r(mean) "; returning ."
        matrix _apexogAPE_vec = J(`_n_endog', 1, .)
        matrix _apexogSE_vec  = J(`_n_endog', 1, .)
        matrix _apexogLO_vec  = J(`_n_endog', 1, .)
        matrix _apexogHI_vec  = J(`_n_endog', 1, .)
        return matrix ape_vec = _apexogAPE_vec
        return matrix se_vec  = _apexogSE_vec
        return matrix lo_vec  = _apexogLO_vec
        return matrix hi_vec  = _apexogHI_vec
        return scalar ape     = .
        return scalar se      = .
        return scalar lo      = .
        return scalar hi      = .
        return scalar mphi    = .
        return scalar stilde  = .
        matrix _apexogVTGT    = J(`_n_endog', `_n_endog', .)
        return matrix vtgt    = _apexogVTGT
        exit
    }

    * Guard 1: the probit must have used the FULL touse sample.  Under perfect
    * prediction/separation, probit silently drops observations, which would
    * misalign the Mata views (built on the full touse) with the fitted
    * coefficients.
    * Guard 2: no endogenous regressor, control function, or CF product may
    * have been omitted as collinear (b = 0 AND se = 0) -- the APE/AIE and the
    * MT cross-Hessian are undefined for an omitted load-bearing column.
    * (Exog controls MAY carry omitted factor base levels; invsym handles
    * those directions, as validated on the factor-control real-data specs.)
    local _bad ""
    qui count if `touse'
    if e(N) < r(N) {
        local _bad "probit used `e(N)' of `r(N)' obs (perfect prediction?)"
    }
    if "`_bad'" == "" {
        foreach _v in `endoglist' `chat0list' `cfprodlist' {
            if _b[`_v'] == 0 & _se[`_v'] == 0 {
                local _bad "`_v' omitted from the probit (collinear)"
            }
        }
    }
    if "`_bad'" != "" {
        di as text "  [warn] `_bad'; returning missing for this fit."
        matrix _apexogAPE_vec = J(`_n_endog', 1, .)
        matrix _apexogSE_vec  = J(`_n_endog', 1, .)
        matrix _apexogLO_vec  = J(`_n_endog', 1, .)
        matrix _apexogHI_vec  = J(`_n_endog', 1, .)
        return matrix ape_vec = _apexogAPE_vec
        return matrix se_vec  = _apexogSE_vec
        return matrix lo_vec  = _apexogLO_vec
        return matrix hi_vec  = _apexogHI_vec
        return scalar ape     = .
        return scalar se      = .
        return scalar lo      = .
        return scalar hi      = .
        return scalar mphi    = .
        return scalar stilde  = .
        matrix _apexogVTGT    = J(`_n_endog', `_n_endog', .)
        return matrix vtgt    = _apexogVTGT
        exit
    }

    * Save probit results
    matrix _apexogPR_B = e(b)
    matrix _apexogPR_V = e(V)

    * Mata: compute MT-corrected V_beta and per-endog APE delta-method SEs.
    * The four matrices below are filled by Mata and returned to the caller.
    if "`debug'" != "" di as text "[debug _apexogfit] entering Mata MT+APE..."
    capture noisily ///
    mata: _apexogfit_apese_mt("`yvar'", "`endoglist'", "`chat0list'",         ///
                            "`offvar'",                                     ///
                            "`aweightlist'", "`clustvar'",                  ///
                            "`exoglist'", "`ivlist'",                       ///
                            "`aiex'", "`aiez'", "`aiexz'", `aiecase',       ///
                            "`cfprodlist'",                                 ///
                            "_apexogAPE_vec", "_apexogSE_vec",                  ///
                            "_apexogLO_vec",  "_apexogHI_vec", `level')
    if _rc != 0 {
        di as error "[_apexogfit] Mata MT+APE failed (rc=`_rc')."
        matrix _apexogAPE_vec = J(`_n_endog', 1, .)
        matrix _apexogSE_vec  = J(`_n_endog', 1, .)
        matrix _apexogLO_vec  = J(`_n_endog', 1, .)
        matrix _apexogHI_vec  = J(`_n_endog', 1, .)
        return matrix ape_vec = _apexogAPE_vec
        return matrix se_vec  = _apexogSE_vec
        return matrix lo_vec  = _apexogLO_vec
        return matrix hi_vec  = _apexogHI_vec
        return scalar ape     = .
        return scalar se      = .
        return scalar lo      = .
        return scalar hi      = .
        return scalar mphi    = .
        return scalar stilde  = .
        matrix _apexogVTGT    = J(`_n_endog', `_n_endog', .)
        return matrix vtgt    = _apexogVTGT
        exit
    }

    * v0.20: exact Pr-scale value of each instrument's violation at this
    * refit (before any `return': _apexoggpr is r-class).
    if "`gammas'" != "" {
        _apexoggpr, touse(`touse') ivlist(`ivlist') gammas(`gammas') ivbin(`ivbin') eq(xb)
        matrix _apexogGPR_vec = r(gpr)
    }

    * Back-compat: also expose the FIRST endog's APE/SE as scalars for
    * existing callers that index r(ape), r(se), r(lo), r(hi).
    if "`gammas'" != "" return matrix gpr_vec = _apexogGPR_vec
    return scalar ape = _apexogAPE_vec[1, 1]
    return scalar se  = _apexogSE_vec[1, 1]
    return scalar lo  = _apexogLO_vec[1, 1]
    return scalar hi  = _apexogHI_vec[1, 1]
    return scalar mphi   = scalar(_apexogMPHI)
    return scalar stilde = scalar(_apexogSTILDE)
    * (v0.18: the debug line must come BEFORE `return matrix', which MOVES the
    *  matrices into r() -- reading them afterwards gave r(111).)
    if "`debug'" != "" di as text "[debug _apexogfit] Mata returned: target APE=" %12.7f _apexogAPE_vec[1,1] " SE=" %12.7f _apexogSE_vec[1,1]
    return matrix vtgt = _apexogVTGT
    return matrix ape_vec = _apexogAPE_vec
    return matrix se_vec  = _apexogSE_vec
    return matrix lo_vec  = _apexogLO_vec
    return matrix hi_vec  = _apexogHI_vec
end

*==============================================================================
* _apexoggpr (v0.20): the EXACT probability-scale value of each instrument's
* violation at the current refit.  Called by _apexogfit right after a
* successful fit, while e() still holds the probit (eq(xb)) or the bivariate
* probit (eq(xb1)) with the offset sum_k gamma_k z_k.  eta = the fitted outcome
* index INCLUDING the offset and the control functions; for instrument k:
*   binary z_k (0/1 in the estimation sample):
*       mean[ Phi(eta - g_k z_k + g_k) - Phi(eta - g_k z_k) ]
*     = the average change in Pr(y=1) when z_k switches 0 -> 1, holding the
*       treatment, the control functions and the other instruments' offsets
*       at their observed values (cf. the discrete-change APE);
*   otherwise:  mean[ phi(eta) ] * g_k   (the average derivative).
* The first-order value gunits(pr) uses is m_phi(baseline) * g_k.  Uses no
* random numbers, so seeded LTZ draws are unaffected.
*==============================================================================
program _apexoggpr, rclass
    syntax, touse(string) ivlist(string) gammas(string) ivbin(string) eq(string)
    tempvar eta w
    tempname G
    if "`eq'" == "xb1" qui predict double `eta' if `touse', xb1
    else               qui predict double `eta' if `touse', xb
    local K : word count `ivlist'
    matrix `G' = J(`K', 1, .)
    qui gen double `w' = . if `touse'
    forvalues k = 1/`K' {
        local z : word `k' of `ivlist'
        local g : word `k' of `gammas'
        local b : word `k' of `ivbin'
        if `b' {
            qui replace `w' = normal(`eta' - (`g')*`z' + (`g')) - normal(`eta' - (`g')*`z') if `touse'
        }
        else {
            qui replace `w' = normalden(`eta') * (`g') if `touse'
        }
        qui summarize `w' if `touse', meanonly
        matrix `G'[`k', 1] = r(mean)
    }
    return matrix gpr = `G'
end

*==============================================================================
* _apexogerr: pretty error message
*==============================================================================
program _apexogerr
    di as error `"`0'"'
    exit 198
end

*==============================================================================
* _apexoghdr (v0.19): the "Estimator:" header line, shared by the UCI and LTZ
* displays, plus the disclosure note for generalized-residual first stages.
*==============================================================================
program _apexoghdr
    syntax, engine(string) firststage(string) n(string) [binlin(string) AIE]
    if "`engine'" == "biprobit" {
        di as text "Estimator: biprobit (recursive bivariate probit, MLE)  |  SE: delta method  |  N = `n'"
        di as text "           outcome and 0/1 treatment fit jointly under bivariate normality;"
        di as text "           gamma is an offset in the outcome equation; no generated regressor."
    }
    else {
        * v0.21.2: "two-step" = the Murphy-Topel parameter covariance plus the
        * direct first-stage term of the effect (the SE computed since 0.21.0)
        di as text "Estimator: cfprobit  |  First stage: `firststage'  |  SE: two-step  |  N = `n'"
        if "`binlin'" != "" {
            * v0.21.2: a 0/1 endogenous regressor with a linear first stage
            di as text "Note:      0/1 endogenous regressor with a linear first stage (`binlin'):"
            di as text "           a binary treatment cannot have the normal linear first stage the"
            di as text "           control-function model assumes, so the effect is a descriptive"
            if "`aie'" != "" {
                * the exact engine and the probit first stage are not available
                * for a 0/1 treatment inside an interaction
                di as text "           fitted contrast without a causal interpretation (see help)."
            }
            else {
                di as text "           fitted contrast without a causal interpretation.  For a causal"
                di as text "           reading use firststage(biprobit) or firststage(probit); see help."
            }
        }
        local _gr = 0
        foreach _f of local firststage {
            if inlist("`_f'", "probit", "fprobit", "poisson") local _gr = 1
        }
        if `_gr' {
            di as text "Note:      the generalized-residual control function (probit/fprobit/"
            di as text "           poisson first stage) is consistent under Wooldridge's sufficiency"
            di as text "           assumption D(eps|z,d) = D(eps|GR) and only approximate under joint"
            di as text "           normality of the structural errors (Wooldridge 2014, 2015; [R]"
            di as text "           cfprobit).  For a 0/1 treatment firststage(biprobit) is exact;"
            di as text "           see help apexog, Methods, for the large-N bias table."
        }
    }
end

*==============================================================================
* _apexogclean: drop the internal _apexog* matrices and scalars.  All results
* live in e() by the time this runs; nothing the user owns matches the prefix.
*==============================================================================
program _apexogclean
    local _mats : all matrices "_apexog*"
    foreach _x of local _mats {
        capture matrix drop `_x'
    }
    local _scals : all scalars "_apexog*"
    foreach _x of local _scals {
        capture scalar drop `_x'
    }
end

*==============================================================================
* Mata library
*==============================================================================
mata:

// ----------------------------------------------------------------------------
// _apexogsim_pct: per-target percentile intervals for the LTZ simulation.
// _apexogAPESIM (N_iter x J) holds the per-draw POINT APEs/AIEs; _apexogSESIM
// (N_iter x J) the per-draw Murphy-Topel SEs.
//
// The HEADLINE interval is the CHR (2012) eq. (7) LTZ confidence interval:
// sampling noise convolved with the gamma prior.  Each refit at gamma_i gives
// the bias-adjusted point estimate (the nonlinear analog of betahat - A*gamma_i),
// so convolving per draw
//        APE*_i = APE_i + SE_i * eps_i ,   eps_i ~ N(0,1)
// and taking percentiles of {APE*_i} reproduces CHR's pivot interval
// [APE0 - c_{1-a/2}(eta), APE0 - c_{a/2}(eta)] exactly (symmetric noise), and
// reproduces plausexog's analytic V + A*Omega*A' interval in the linear model
// (checked by simulation).
//
// The eps draws happen HERE, after the gamma loop, so the gamma-draw sequence
// under seed() is identical to pre-0.14 versions.  A full N_iter x J eps
// matrix is drawn regardless of failed (missing) draws, so the RNG state
// consumed is reproducible.  For J > 1 the eps draws are independent across
// targets: exact for the per-target MARGINAL intervals reported (joint
// inference across targets would need the cross-APE covariance, not returned).
//
// Writes:
//   _apexogLO_PCT  / _apexogHI_PCT   (J x 1)  convolved percentiles (headline)
//   _apexogLO_PRIOR/ _apexogHI_PRIOR (J x 1)  prior-only point-APE percentiles
//   _apexogSIMCNV  (N_iter x J)               the convolved draws
//   locals lo_pct / hi_pct                    first target, convolved
// ----------------------------------------------------------------------------
void _apexogsim_pct(real scalar level)
{
    real matrix      S, SE, C
    real colvector   col_clean, lo_vec, hi_vec, lo_raw, hi_raw
    real scalar      n_iter, J, j, n_valid, lo_pos, hi_pos, alpha
    S      = st_matrix("_apexogAPESIM")
    SE     = st_matrix("_apexogSESIM")
    n_iter = rows(S)
    J      = cols(S)
    alpha  = 1 - level
    // convolve: missing APE_i or SE_i propagates to a missing convolved draw
    C      = S + SE :* rnormal(n_iter, J, 0, 1)
    lo_vec = J(J, 1, .)
    hi_vec = J(J, 1, .)
    lo_raw = J(J, 1, .)
    hi_raw = J(J, 1, .)
    for (j = 1; j <= J; j++) {
        // headline: convolved draws
        col_clean = select(C[., j], !missing(C[., j]))
        n_valid   = rows(col_clean)
        if (n_valid > 0) {
            col_clean = sort(col_clean, 1)
            lo_pos    = max((1, round(alpha/2 * n_valid)))
            hi_pos    = min((n_valid, round((1 - alpha/2) * n_valid)))
            lo_vec[j] = col_clean[lo_pos]
            hi_vec[j] = col_clean[hi_pos]
        }
        // secondary: prior-only spread of the point APEs (pre-0.14 headline)
        col_clean = select(S[., j], !missing(S[., j]))
        n_valid   = rows(col_clean)
        if (n_valid > 0) {
            col_clean = sort(col_clean, 1)
            lo_pos    = max((1, round(alpha/2 * n_valid)))
            hi_pos    = min((n_valid, round((1 - alpha/2) * n_valid)))
            lo_raw[j] = col_clean[lo_pos]
            hi_raw[j] = col_clean[hi_pos]
        }
    }
    st_matrix("_apexogLO_PCT",   lo_vec)
    st_matrix("_apexogHI_PCT",   hi_vec)
    st_matrix("_apexogLO_PRIOR", lo_raw)
    st_matrix("_apexogHI_PRIOR", hi_raw)
    st_matrix("_apexogSIMCNV",   C)
    // Back-compat: locals for the first target (now the convolved interval)
    st_local("lo_pct", strofreal(lo_vec[1], "%21.15g"))
    st_local("hi_pct", strofreal(hi_vec[1], "%21.15g"))
}

// ----------------------------------------------------------------------------
// _apexogsweep_pct: percentile intervals for ONE delta-sweep step.  Reads the
// n_iter x 1 matrices _apexogSWAPE (point APEs/AIEs) and _apexogSWSE (MT SEs),
// convolves exactly as _apexogsim_pct does, and writes four scalars:
// _apexogSW_lo/_apexogSW_hi (convolved CI) and _apexogSW_lop/_apexogSW_hip
// (prior-only band).  Uses its own matrix names so the main LTZ results are
// never clobbered.
// ----------------------------------------------------------------------------
void _apexogsweep_pct(real scalar level)
{
    real colvector ape, se, T, s
    real matrix    R
    real scalar    n, lo_pos, hi_pos, a
    R   = st_matrix("_apexogSWAPE"), st_matrix("_apexogSWSE")
    R   = select(R, !rowmissing(R))
    n   = rows(R)
    a   = 1 - level
    if (n == 0) {
        st_numscalar("_apexogSW_lo",  .)
        st_numscalar("_apexogSW_hi",  .)
        st_numscalar("_apexogSW_lop", .)
        st_numscalar("_apexogSW_hip", .)
        return
    }
    ape = R[., 1]
    se  = R[., 2]
    T   = ape + se :* rnormal(n, 1, 0, 1)
    lo_pos = max((1, round(a/2 * n)))
    hi_pos = min((n, round((1 - a/2) * n)))
    s = sort(T, 1)
    st_numscalar("_apexogSW_lo",  s[lo_pos])
    st_numscalar("_apexogSW_hi",  s[hi_pos])
    s = sort(ape, 1)
    st_numscalar("_apexogSW_lop", s[lo_pos])
    st_numscalar("_apexogSW_hip", s[hi_pos])
}

// ----------------------------------------------------------------------------
// _apexogfit_meat (v0.21.0): outer product of per-observation contributions,
// summed within clusters first when clustered == 1 (sort_idx / info from
// order() and panelsetup() on the cluster id).
// ----------------------------------------------------------------------------
real matrix _apexogfit_meat(real matrix M, real matrix sort_idx,
                            real matrix info, real scalar clustered)
{
    real matrix Ms, Mc
    real scalar ci, m
    if (!clustered) return(cross(M, M))
    m  = cols(M)
    Ms = M[sort_idx, .]
    Mc = J(rows(info), m, 0)
    for (ci = 1; ci <= rows(info); ci++) {
        Mc[ci, .] = colsum(Ms[|info[ci,1], 1 \ info[ci,2], m|])
    }
    return(cross(Mc, Mc))
}

// ----------------------------------------------------------------------------
// _apexogfit_apese_mt: Murphy-Topel-corrected V_beta + APE delta-method SE.
//
// Unified for linear and probit first stages via a per-obs weight A_i:
//   linear:  A_i = 1                                  =>  w_dchat_i = -A_i = -1
//   probit:  A_i = chat0_i * (W_i'*pi + chat0_i)      =>  w_dchat_i = -A_i
//   fprobit: A_i = chat0_i * (W_i'*pi + chat0_i)
//                  + d_i(1-d_i)*[phi/(Phi(1-Phi))]^2   =>  w_dchat_i = -A_i
//   poisson: A_i = exp(W_i'*pi)                       =>  w_dchat_i = -A_i
// where w_dchat_i = d(chat0_i) / d(W_i'*pi) -- the per-obs derivative of the
// control function with respect to the first-stage linear predictor.
//
// Model:
//   First stage:    chat0_i is a known function f(d_i, W_i'*pi).  W = [Z, X, 1].
//   Main probit:    P(y_i=1) = Phi(eta_i)
//                   eta_i = X_i' b_X + d_i*b_d + chat0_i*b_chat + b_const + offvar_i
//                   x_main_i = [X_i, d_i, chat0_i, 1]
//
//   Probit score:   s_{2,i} = psi_i * x_main_i,
//                   where psi_i = (y_i - Phi(eta_i)) * phi(eta_i) / [Phi*(1-Phi)]
//   Probit Hess:    H_b = sum_i psi_i*(eta_i + psi_i) * x_main_i x_main_i'
//
//   Cross-Hess J_{2,pi} (the unified form, valid for any first stage):
//     J_{2,pi} = -b_chat * sum_i psi_i(eta+psi_i) * w_dchat_i * x_main_i W_i'
//                + e_{chat0} * sum_i psi_i * w_dchat_i * W_i'
//
//   First-stage bread:  M_pi = (W' * diag(A) * W)^{-1}
//   First-stage IF:      psi_{pi,i} = M_pi * W_i * chat0_i
//      (linear: M_pi = (W'W)^{-1}, identical to OLS;
//       probit: M_pi = (-H_1)^{-1} where H_1 is the first-stage probit Hessian)
//
//   Total score:    s_total_i = s_{2,i} + J_{2,pi} * psi_{pi,i}
//   MT VCE:         V_beta_MT = H_b^{-1} * (sum_i s_total_i s_total_i') * H_b^{-1}
//
//   APE:            APE = mean(phi(eta_i)) * b_d
//   APE Jacobian J_beta (length cols(x_main); aligned to probit e(b) order):
//     d/db_const = b_d * mean(-eta*phi)
//     d/db_X_v   = b_d * mean(-eta*phi * X_v)
//     d/db_d     = b_d * mean(-eta*phi * d) + mean(phi)
//     d/db_chat  = b_d * mean(-eta*phi * chat0)
//
//   v0.21.0 -- full two-step influence function of the target.  The effect
//   q_i depends on pi-hat not only through the probit coefficients (the
//   Murphy-Topel channel above) but DIRECTLY through chat0_i in the index:
//     J_pi_k  = mean( Dq_i * g_k,i * w_dchat_k,i * W_k,i' ),
//               Dq_i = d q_i / d(index shift), g_k,i = d eta_i / d chat_k,i
//     l_i     = psi_2,i' J_beta' + sum_k psi_pi_k,i' J_pi_k'
//               psi_2,i = H_b^{-1} s_total_i   (the MT influence of beta-hat)
//     SE      = sqrt( sum_i l_i^2 )   (cluster-summed l_i under vce(cluster))
//   Through v0.20.1 only the first term was used (so does cfprobit + margins,
//   which treats the stored control function as a fixed regressor); the
//   omission is about 1% of the SE in the validation designs.  Checked by
//   _sanity/59_direct_channel.do against a formula-free oracle.
//   level x 100 % CI: APE +/- invnormal(1 - alpha/2) * SE
// ----------------------------------------------------------------------------
// Multi-endog, multi-target (v0.8).
// x_main = [X, d_1..d_J, chat0_1..chat0_J, _cons].
// Returns APE and SE for EVERY endog in `endognames' (each one is the target
// in its own delta-method computation, with V_beta_MT shared).  Results
// written into Stata matrices via the m*name arguments (J x 1 each).
void _apexogfit_apese_mt(string scalar yvar,
                       string scalar endognames,   string scalar chatnames,
                       string scalar offvarname,
                       string scalar Anames,       string scalar clustvarname,
                       string scalar exognames,    string scalar ivnames,
                       string scalar aiex_name,    string scalar aiez_name,
                       string scalar aiexz_name,   real scalar   aie_case,
                       string scalar cfprodnames,
                       string scalar mname_ape,    string scalar mname_se,
                       string scalar mname_lo,     string scalar mname_hi,
                       real scalar   level)
{
    string rowvector evars, zvars, endogs, chats, Avars, cfprods
    real matrix      X_exog, Z_iv, D, CHAT, AA, CFPROD, x_main, W
    real matrix      H_b, U_2, V_beta_MT, meat, s_2, s_total
    real matrix      M_pi_k, psi_pi_k, A1_k, J2pi_k
    real matrix      S_clust, cluster_ids, info
    real matrix      sort_idx, s_sorted, cluster_sorted
    real colvector   y, offset_vals, A_k, chat_k, eta, phi_eta, Phi_eta, Phim_eta
    real colvector   psi, psi_eta_psi, w_dchat_k, A2_k, ones_col, mean_term
    real colvector   eta_1, eta_0, Phi_1, Phi_0, phi_1, phi_0, diff_phi
    real colvector   g_k, cf_a, cf_b
    real rowvector   b_pr, J_beta, J_pi
    real matrix      Jall, Vtgt, PSI1, GWD, psi2, L
    real colvector   Dq
    real scalar      clustered
    real scalar      n_exog, n_iv, n_endog, n_obs, nx, nw, j, k, ci, n_cf
    real scalar      cf_start_idx, cc, aa, bb, m_partner, cf_col_k
    real scalar      d_start_idx, chat_start_idx, cons_idx
    real scalar      d_target_idx, chat_target_idx
    real scalar      b_chat_k, b_d_target
    real scalar      mean_phi, ape, se, se_sq, z_crit
    real scalar      binary_d_target
    string scalar    touse

    evars   = (exognames == "") ? J(1, 0, "") : tokens(exognames)
    zvars   = tokens(ivnames)
    endogs  = tokens(endognames)
    chats   = tokens(chatnames)
    Avars   = tokens(Anames)
    n_exog  = cols(evars)
    n_iv    = cols(zvars)
    n_endog = cols(endogs)
    touse   = st_local("touse")

    if (cols(chats) != n_endog | cols(Avars) != n_endog) {
        _error("_apexogfit_apese_mt: chat0list and aweightlist must each have one element per endogenous regressor.")
    }

    // ---------- pull data ----------
    st_view(y,           ., yvar,       touse)
    st_view(D,           ., endogs,     touse)        // N x n_endog
    st_view(CHAT,        ., chats,      touse)        // N x n_endog
    st_view(AA,          ., Avars,      touse)        // N x n_endog
    st_view(offset_vals, ., offvarname, touse)
    if (n_exog > 0) st_view(X_exog, ., evars, touse)
    if (n_iv > 0)   st_view(Z_iv,   ., zvars, touse)
    n_obs    = rows(y)
    ones_col = J(n_obs, 1, 1)

    // cfinteract product columns (chat_i * chat_j), if any
    cfprods = (cfprodnames == "") ? J(1, 0, "") : tokens(cfprodnames)
    n_cf    = cols(cfprods)
    if (n_cf > 0) st_view(CFPROD, ., cfprods, touse)

    // x_main column order MUST match Stata's probit:
    //   y EXOG d_1..d_J chat0_1..chat0_J [cfprod_1..cfprod_C] _cons
    // e(b) order: [X_1..X_p, d_1..d_J, chat0_1..chat0_J, cfprod_1..C, _cons]
    if (n_exog > 0) x_main = (X_exog, D, CHAT)
    else            x_main = (D, CHAT)
    if (n_cf > 0)   x_main = (x_main, CFPROD)
    x_main = (x_main, ones_col)
    nx = cols(x_main)

    // Global first-stage regressor matrix (union of all IVs across clauses
    // plus the shared exog and constant):
    //   W_full = [Z_1..Z_m, X_1..X_p, 1]
    // For each endog k, the per-endog W_k is a SUBMATRIX selected by an IV-mask
    // (read from the Stata matrix _apexogIVMASK, J x n_iv) AND an EXOG-mask (read
    // from _apexogEXOGMASK, J x n_exog).  The masks say which global IVs and
    // which global exog columns each endog's first-stage regression actually
    // uses.  EXOGMASK is all-1 in normal cases; in aie() with an endog
    // interactor, the aie_XZ column is masked out for that endog (otherwise
    // the first stage would be circular: d on the LHS, X*d on the RHS).
    if (n_iv > 0 & n_exog > 0)       W = (Z_iv, X_exog, ones_col)
    else if (n_iv > 0 & n_exog == 0) W = (Z_iv, ones_col)
    else if (n_iv == 0 & n_exog > 0) W = (X_exog, ones_col)
    else                             W = ones_col
    nw = cols(W)
    real matrix IVMASK, EXOGMASK
    if (n_iv > 0) IVMASK   = st_matrix("_apexogIVMASK")
    else          IVMASK   = J(n_endog, 0, .)
    if (n_exog > 0) EXOGMASK = st_matrix("_apexogEXOGMASK")
    else            EXOGMASK = J(n_endog, 0, .)

    // Column indices into x_main
    d_start_idx     = n_exog + 1
    chat_start_idx  = n_exog + n_endog + 1
    cf_start_idx    = chat_start_idx + n_endog       // first cfprod column (if any)
    cons_idx        = nx
    d_target_idx    = d_start_idx                   // first endog = the target
    chat_target_idx = chat_start_idx                // its CF column

    // Map each cfinteract product column c -> the endog pair (a,b) it multiplies.
    // The wrapper builds ALL unique pairwise products of the control functions in
    // lexicographic order, so column c of CFPROD = chat_{cf_a[c]} * chat_{cf_b[c]}.
    cf_a = J(n_cf, 1, .)
    cf_b = J(n_cf, 1, .)
    if (n_cf > 0) {
        cc = 0
        for (aa = 1; aa <= n_endog - 1; aa++) {
            for (bb = aa + 1; bb <= n_endog; bb++) {
                cc++
                if (cc <= n_cf) {
                    cf_a[cc] = aa
                    cf_b[cc] = bb
                }
            }
        }
    }

    // ---------- pull probit coefficients ----------
    b_pr = st_matrix("_apexogPR_B")
    if (cols(b_pr) != nx) {
        _error("_apexogfit_apese_mt: probit e(b) cols (" + strofreal(cols(b_pr)) +
               ") != expected nx (" + strofreal(nx) + ").  Check that the probit was fit on " +
               "[exog endog_1..endog_J chat0_1..chat0_J + constant].")
    }
    b_d_target = b_pr[d_target_idx]

    // ---------- compute eta, psi, psi(eta+psi) ----------
    eta      = x_main * b_pr' + offset_vals
    phi_eta  = normalden(eta)
    Phi_eta  = normal(eta)
    Phim_eta = normal(-eta)
    // average dPr/dindex slope: converts latent-index units to probability
    // points (used by gunits(pr) and reported as an elicitation aid)
    st_numscalar("_apexogMPHI", mean(phi_eta))
    // sigma_tilde: probit-index-to-STRUCTURAL-latent scale ratio.  In the
    // structural equation y* = xb + gamma*z + eps with TOTAL error variance
    // normalized to 1, the probit instead normalizes the CONDITIONAL (on the
    // control functions) error sd to 1, so index units = structural units /
    // sigma_tilde, sigma_tilde = sd(eps | CF).  From the CF loadings b_cf
    // and the sample covariance of the control functions:
    //   sigma_tilde = (1 + b_cf' Var(chat) b_cf)^(-1/2)
    // which reduces to (1 + theta_cf^2)^(-1/2) for one CF with unit-variance
    // first-stage residual.  Derived for, and exact with, LINEAR first
    // stages (chat estimates the structural first-stage error).  For the
    // generalized-residual first stages the same formula is evaluated on
    // the generalized residual, but it is NOT the index-scale conversion
    // there: the validation Monte Carlos put the
    // correct offset at the unscaled structural gamma, so v0.19.1 warns
    // when gunits(structural) is combined with such a first stage (the
    // value is still reported in e(stilde), as a baseline-fit quantity).
    // cfinteract CF-product loadings are ignored (first-order in the CF
    // coefficients).  Used by gunits(structural) and reported as an
    // elicitation aid.
    real colvector b_cf_col
    b_cf_col = b_pr[|chat_start_idx \ chat_start_idx + n_endog - 1|]'
    st_numscalar("_apexogSTILDE",
        1/sqrt(1 + b_cf_col' * variance(CHAT) * b_cf_col))
    psi      = (y :- Phi_eta) :* phi_eta :/ (Phi_eta :* Phim_eta)
    psi_eta_psi = psi :* (eta :+ psi)

    // ---------- probit Hessian H_b ----------
    H_b = cross(x_main, psi_eta_psi, x_main)
    U_2 = invsym(H_b)

    // ---------- accumulate first-stage contributions ----------
    // For each endog k:
    //   chat_k = CHAT[., k];  A_k = AA[., k];  w_dchat_k = -A_k
    //   M_pi_k    = (W' diag(A_k) W)^{-1}
    //   psi_pi_k  = (W :* chat_k) * M_pi_k'        // N x nw
    //   A1_k = cross(x_main, psi_eta_psi :* w_dchat_k, W)
    //   A2_k = cross(W,        psi          :* w_dchat_k)
    //   J2pi_k = -b_chat_k :* A1_k
    //   J2pi_k[chat_idx_k, .] += A2_k'
    //   s_total += psi_pi_k * J2pi_k'
    s_2     = x_main :* psi                          // N x nx
    s_total = s_2
    real matrix W_k
    real vector active_z, active_x, active_cols_k
    real scalar nz_k, nx_k
    for (k = 1; k <= n_endog; k++) {
        // Build per-endog W_k by selecting columns of W_full according to masks.
        if (n_iv > 0) active_z = selectindex(IVMASK[k, .])
        else          active_z = J(1, 0, .)
        nz_k = cols(active_z)
        if (n_exog > 0) active_x = selectindex(EXOGMASK[k, .])
        else            active_x = J(1, 0, .)
        nx_k = cols(active_x)
        // active_cols_k = [active Z indices, active exog indices (offset by n_iv), constant idx]
        active_cols_k = J(1, 0, .)
        if (nz_k > 0) active_cols_k = active_z
        if (nx_k > 0) active_cols_k = (active_cols_k, n_iv :+ active_x)
        active_cols_k = (active_cols_k, nw)         // the constant column
        W_k = W[., active_cols_k]

        chat_k    = CHAT[., k]
        A_k       = AA[., k]
        w_dchat_k = -A_k
        b_chat_k  = b_pr[chat_start_idx + k - 1]

        M_pi_k   = invsym(cross(W_k, A_k, W_k))
        psi_pi_k = (W_k :* chat_k) * M_pi_k'        // N x nw_k

        // g_k = d(eta)/d(linear-pred-k) coefficient.  Normally just b_chat_k.
        // With cfinteract, every product chat_k*chat_m also depends on pi_k
        // (through chat_k), so eta gains b_cf*chat_m*d(chat_k); hence
        // g_k,i = b_chat_k + sum_{products p with k} b_cf_p * chat_m,i.
        // Reduces to the scalar b_chat_k when n_cf == 0 (backward compatible).
        g_k = J(n_obs, 1, b_chat_k)
        for (cc = 1; cc <= n_cf; cc++) {
            m_partner = .
            if      (cf_a[cc] == k) m_partner = cf_b[cc]
            else if (cf_b[cc] == k) m_partner = cf_a[cc]
            if (m_partner != .) {
                cf_col_k = cf_start_idx + cc - 1
                g_k = g_k + b_pr[cf_col_k] :* CHAT[., m_partner]
            }
        }

        A1_k = cross(x_main, psi_eta_psi :* g_k :* w_dchat_k, W_k)   // nx x nw_k
        A2_k = cross(W_k,      psi          :* w_dchat_k)            // nw_k x 1
        J2pi_k = -A1_k
        J2pi_k[chat_start_idx + k - 1, .] = J2pi_k[chat_start_idx + k - 1, .] + A2_k'
        // d(cfprod_c)/d(pi_k) = chat_m * d(chat_k)/d(pi_k) adds an A2-type row at
        // each product column the endog k participates in.
        for (cc = 1; cc <= n_cf; cc++) {
            m_partner = .
            if      (cf_a[cc] == k) m_partner = cf_b[cc]
            else if (cf_b[cc] == k) m_partner = cf_a[cc]
            if (m_partner != .) {
                cf_col_k = cf_start_idx + cc - 1
                J2pi_k[cf_col_k, .] = J2pi_k[cf_col_k, .] + cross(W_k, psi :* CHAT[., m_partner] :* w_dchat_k)'
            }
        }

        s_total = s_total + psi_pi_k * J2pi_k'       // N x nx (product)

        // v0.21.0: keep the first-stage influence functions and the per-obs
        // derivative of the index w.r.t. pi_k (d eta_i/d pi_k' = g_k,i *
        // w_dchat_k,i * W_k,i') for the direct term of the target's variance.
        if (k == 1) {
            PSI1 = psi_pi_k
            GWD  = (g_k :* w_dchat_k) :* W_k
        }
        else {
            PSI1 = (PSI1, psi_pi_k)
            GWD  = (GWD, (g_k :* w_dchat_k) :* W_k)
        }
    }

    // ---------- meat: cluster-robust or sandwich ----------
    clustered = (clustvarname != "")
    if (clustered) {
        st_view(cluster_ids, ., clustvarname, touse)
        sort_idx       = order(cluster_ids, 1)
        cluster_sorted = cluster_ids[sort_idx, .]
        info = panelsetup(cluster_sorted, 1)
    }
    else {
        sort_idx = J(0, 0, .)
        info     = J(0, 0, .)
    }
    meat      = _apexogfit_meat(s_total, sort_idx, info, clustered)
    V_beta_MT = U_2 * meat * U_2
    // v0.21.0: per-obs influence of the probit coefficients (N x nx); the
    // target's full influence is psi2 * J_beta' + PSI1 * J_pi'.
    psi2 = s_total * U_2

    // ---------- TARGET FUNCTIONAL: APE (per endog) OR AIE (single) ----------
    real colvector ape_vec, se_vec, lo_vec, hi_vec
    real scalar    t, b_d_t, target_idx_t, chat_target_idx_t, binary_d_t
    real scalar    se_t, se_sq_t
    z_crit  = invnormal(1 - (1 - level)/2)

    if (aiex_name != "") {
        // ====== AIE branch: three cases ======
        // aie_case: 0 = cts x cts, 1 = cts x bin (Z bin), 2 = bin x cts (X bin),
        //           3 = bin x bin
        //
        // Shared setup: locate X, Z in either evars or endogs.  x_main column
        // order = [evars, endogs, chats, _cons], so:
        //   exog var at evars position k -> x_main column k
        //   endog var at endogs position k -> x_main column n_exog + k
        // The XZ tempvar is in evars when the interaction is an exogenous
        // control (aie() all-exog case) OR in endogs when it was declared as an
        // endogenous regressor (the endogenous-treatment x moderator case).
        real scalar X_xmcol, Z_xmcol, XZ_idx, XZ_endogk, k_search
        real colvector X_val, Z_val, XZ_val
        real rowvector J_beta_aie
        real scalar b_X_coef, b_Z_coef, b_XZ_coef, AIE_pt, se_aie, se_sq_aie
        X_xmcol = .; Z_xmcol = .; XZ_idx = .; XZ_endogk = .
        X_val   = J(n_obs, 1, .)
        Z_val   = J(n_obs, 1, .)
        for (k_search = 1; k_search <= n_exog; k_search++) {
            if (evars[k_search] == aiex_name) {
                X_xmcol = k_search
                X_val   = X_exog[., k_search]
            }
            if (evars[k_search] == aiez_name) {
                Z_xmcol = k_search
                Z_val   = X_exog[., k_search]
            }
            if (evars[k_search] == aiexz_name) XZ_idx = k_search
        }
        for (k_search = 1; k_search <= n_endog; k_search++) {
            if (endogs[k_search] == aiex_name) {
                X_xmcol = n_exog + k_search
                X_val   = D[., k_search]
            }
            if (endogs[k_search] == aiez_name) {
                Z_xmcol = n_exog + k_search
                Z_val   = D[., k_search]
            }
            if (endogs[k_search] == aiexz_name) {
                XZ_idx    = n_exog + k_search
                XZ_endogk = k_search
            }
        }
        if (X_xmcol == . | Z_xmcol == . | XZ_idx == .) {
            _error("_apexogfit_apese_mt: aie() variables not found in exog/endog list.")
        }
        if (XZ_endogk != .) XZ_val = D[., XZ_endogk]
        else                XZ_val = X_exog[., XZ_idx]
        b_X_coef  = b_pr[X_xmcol]
        b_Z_coef  = b_pr[Z_xmcol]
        b_XZ_coef = b_pr[XZ_idx]

        // Branch on case.  Each branch fills AIE_pt and J_beta_aie (1 x nx).
        J_beta_aie = J(1, nx, 0)

        if (aie_case == 0) {
            // -------- CASE 0: continuous-by-continuous (cross-partial) --------
            //   AIE_i = phi(eta) * (b_XZ - eta*a*b),
            //   a = b_X + b_XZ * Z,  b = b_Z + b_XZ * X
            real colvector a_aie, b_aie, w_aie, T_aie, AIE_i_c, phi_T_aie
            a_aie  = b_X_coef :+ b_XZ_coef :* Z_val
            b_aie  = b_Z_coef :+ b_XZ_coef :* X_val
            w_aie  = b_XZ_coef :- eta :* a_aie :* b_aie
            AIE_i_c = phi_eta :* w_aie
            AIE_pt  = mean(AIE_i_c)

            T_aie     = -eta :* w_aie :- a_aie :* b_aie
            phi_T_aie = phi_eta :* T_aie
            Dq        = phi_T_aie                      // d AIE_i / d(index shift)
            for (j = 1; j <= nx; j++) {
                J_beta_aie[j] = mean(phi_T_aie :* x_main[., j])
            }
            J_beta_aie[X_xmcol] = J_beta_aie[X_xmcol] - mean(phi_eta :* eta :* b_aie)
            J_beta_aie[Z_xmcol] = J_beta_aie[Z_xmcol] - mean(phi_eta :* eta :* a_aie)
            J_beta_aie[XZ_idx]  = J_beta_aie[XZ_idx] + mean(phi_eta)             ///
                                - mean(phi_eta :* eta :* Z_val :* b_aie)         ///
                                - mean(phi_eta :* eta :* a_aie :* X_val)
        }
        else if (aie_case == 1 | aie_case == 2) {
            // -------- CASE 1 or 2: cts-by-binary (or binary-by-cts) ----------
            // Without loss of generality, let CTS = the continuous interactor
            // and BIN = the binary one.  Define
            //   eta^B1 = eta + (1 - B_i) * (b_BIN + b_XZ * C_i)
            //   eta^B0 = eta -       B_i * (b_BIN + b_XZ * C_i)
            //   AIE_i  = phi(eta^B1) * (b_CTS + b_XZ) - phi(eta^B0) * b_CTS
            //          [continuous-derivative differenced across the binary]
            // (The CTS derivative is taken at C=C_i; the BIN difference is
            //  D(B=1) - D(B=0) of that derivative.)
            real scalar CTS_xmcol, BIN_xmcol, b_CTS, b_BIN
            real colvector CTS_val, BIN_val, eta_b1, eta_b0
            real colvector phi_b1, phi_b0, AIE_i_cb
            real colvector A_cb, B_cb, w_pen, w_act
            if (aie_case == 1) {                 // X cts, Z bin
                CTS_xmcol = X_xmcol; BIN_xmcol = Z_xmcol
                CTS_val   = X_val  ; BIN_val   = Z_val
                b_CTS     = b_X_coef; b_BIN    = b_Z_coef
            }
            else {                                // X bin, Z cts
                CTS_xmcol = Z_xmcol; BIN_xmcol = X_xmcol
                CTS_val   = Z_val  ; BIN_val   = X_val
                b_CTS     = b_Z_coef; b_BIN    = b_X_coef
            }
            // w_pen = b_BIN + b_XZ * CTS_val  (the shift to eta when B flips)
            w_pen = b_BIN :+ b_XZ_coef :* CTS_val
            eta_b1 = eta + (1 :- BIN_val) :* w_pen
            eta_b0 = eta -        BIN_val  :* w_pen
            phi_b1 = normalden(eta_b1)
            phi_b0 = normalden(eta_b0)
            // act = b_CTS + b_XZ * B (the active derivative coefficient)
            //   at B=1: act = b_CTS + b_XZ
            //   at B=0: act = b_CTS
            AIE_i_cb = phi_b1 :* (b_CTS + b_XZ_coef) :- phi_b0 :* b_CTS
            AIE_pt    = mean(AIE_i_cb)

            // Jacobian.  Build "generic" piece phi_b1*(...)*x_main[k] etc., and
            // add corrections for CTS_xmcol, BIN_xmcol, XZ_idx.
            //
            // ∂AIE_i / ∂beta_k =
            //     -eta_b1*phi_b1 * (b_CTS+b_XZ) * (∂eta_b1/∂beta_k)
            //     + phi_b1 * (∂(b_CTS+b_XZ)/∂beta_k)
            //     - (-eta_b0*phi_b0 * b_CTS * (∂eta_b0/∂beta_k))
            //     - phi_b0 * (∂b_CTS/∂beta_k)
            //
            // For most k: ∂eta_b1/∂beta_k = ∂eta_b0/∂beta_k = x_main[k]
            //             and ∂b_CTS/∂beta_k = ∂b_XZ/∂beta_k = 0
            // So generic[k] = x_main[k] * (-eta_b1*phi_b1*(b_CTS+b_XZ) + eta_b0*phi_b0*b_CTS)
            //
            // Special k:
            //  CTS_xmcol: ∂b_CTS/∂beta = 1; generic still holds for eta-derivs.
            //             Add: + phi_b1 - phi_b0  (from the explicit b_CTS terms)
            //  BIN_xmcol: ∂eta_b1/∂b_BIN = ∂eta/∂b_BIN + (1-B) = B + (1-B) = 1
            //             ∂eta_b0/∂b_BIN = ∂eta/∂b_BIN - B    = B - B     = 0
            //             So the b_BIN entry is:
            //               -1 * eta_b1 * phi_b1 * (b_CTS+b_XZ)
            //             (no contribution from the eta_b0 term, no explicit
            //              b_CTS/b_XZ change)
            //  XZ_idx: ∂eta_b1/∂b_XZ = X_i*Z_i + (1-B)*CTS = CTS_i*BIN_i + (1-BIN)*CTS = CTS
            //          ∂eta_b0/∂b_XZ = X_i*Z_i - B*CTS = 0
            //          ∂b_XZ/∂b_XZ = 1 (only affects the b_CTS+b_XZ piece in eta_b1's coefficient)
            //          So:
            //            generic_xz_using_x_main_XZ = -eta_b1*phi_b1*(b_CTS+b_XZ)*X_i*Z_i + eta_b0*phi_b0*b_CTS*X_i*Z_i
            //            correction = phi_b1 * 1                                            (from ∂b_XZ/∂b_XZ in phi_b1 coeff)
            //                       + [-eta_b1*phi_b1*(b_CTS+b_XZ)*(CTS - X*Z)             // ∂eta_b1 diff
            //                          + eta_b0*phi_b0*b_CTS * X*Z]                        // ∂eta_b0 diff
            //
            // For clean accounting, recompute the J_beta entries for the three
            // special columns explicitly (overriding any earlier generic value).
            real colvector G_i
            G_i = -eta_b1 :* phi_b1 :* (b_CTS + b_XZ_coef) :+ eta_b0 :* phi_b0 :* b_CTS
            Dq  = G_i                                  // d AIE_i / d(index shift)
            for (j = 1; j <= nx; j++) {
                J_beta_aie[j] = mean(G_i :* x_main[., j])
            }
            // CTS correction: add phi_b1 - phi_b0 (the explicit b_CTS bump).
            J_beta_aie[CTS_xmcol] = J_beta_aie[CTS_xmcol] + mean(phi_b1) - mean(phi_b0)
            // BIN correction: the b_BIN entry has a non-generic d(eta_b1)/d(b_BIN)=1
            // and d(eta_b0)/d(b_BIN)=0, so the generic value (which used
            // x_main[BIN_xmcol]=BIN_val) is replaced.
            // Note: must SUBTRACT what we already put in (mean(G_i * BIN_val))
            // and ADD the correct mean(G_only_via_eta_b1).  Since the b_BIN
            // entry only contributes via the eta_b1 piece (with derivative 1)
            // -- NOT via eta_b0 -- the correct value is:
            //   mean( -eta_b1 * phi_b1 * (b_CTS+b_XZ) * 1 )
            //  = mean( -eta_b1 * phi_b1 * (b_CTS+b_XZ) ).
            J_beta_aie[BIN_xmcol] = mean(-eta_b1 :* phi_b1 :* (b_CTS + b_XZ_coef))
            // XZ_idx correction (replace generic):
            //   ∂AIE_i/∂b_XZ = -eta_b1*phi_b1*(b_CTS+b_XZ) * CTS_i + phi_b1
            //                   (the b_XZ term in the b_CTS+b_XZ multiplier kicks in
            //                    via the explicit coefficient deriv)
            //                  - 0   (eta_b0 derivative is 0)
            J_beta_aie[XZ_idx] = mean(-eta_b1 :* phi_b1 :* (b_CTS + b_XZ_coef) :* CTS_val) ///
                                + mean(phi_b1)
        }
        else if (aie_case == 3) {
            // -------- CASE 3: binary-by-binary (2nd discrete difference) ------
            //   AIE_i = Phi(eta|X=1,Z=1) - Phi(eta|X=1,Z=0)
            //         - Phi(eta|X=0,Z=1) + Phi(eta|X=0,Z=0)
            // Build the four counterfactual etas:
            //   eta^{xz} = eta + (x - X_i)*b_X + (z - Z_i)*b_Z + (x*z - X_i*Z_i)*b_XZ
            real colvector eta_11, eta_10, eta_01, eta_00
            real colvector Phi_11, Phi_10, Phi_01, Phi_00
            real colvector phi_11, phi_10, phi_01, phi_00
            real colvector AIE_i_bb
            real colvector d_eta_X1Z1, d_eta_X1Z0, d_eta_X0Z1, d_eta_X0Z0  // per-obs deta/dbeta_k for generic k
            real colvector S_xx, S_xz, S_zz  // accumulators for X, Z, XZ specials

            eta_11 = eta + (1 :- X_val):*b_X_coef + (1 :- Z_val):*b_Z_coef + (1 :- X_val:*Z_val):*b_XZ_coef
            eta_10 = eta + (1 :- X_val):*b_X_coef + (0 :- Z_val):*b_Z_coef + (0 :- X_val:*Z_val):*b_XZ_coef
            eta_01 = eta + (0 :- X_val):*b_X_coef + (1 :- Z_val):*b_Z_coef + (0 :- X_val:*Z_val):*b_XZ_coef
            eta_00 = eta + (0 :- X_val):*b_X_coef + (0 :- Z_val):*b_Z_coef + (0 :- X_val:*Z_val):*b_XZ_coef
            Phi_11 = normal(eta_11)
            Phi_10 = normal(eta_10)
            Phi_01 = normal(eta_01)
            Phi_00 = normal(eta_00)
            phi_11 = normalden(eta_11)
            phi_10 = normalden(eta_10)
            phi_01 = normalden(eta_01)
            phi_00 = normalden(eta_00)
            AIE_i_bb = Phi_11 :- Phi_10 :- Phi_01 :+ Phi_00
            AIE_pt   = mean(AIE_i_bb)

            // Jacobian.  For each k that affects eta via x_main[., k] in the
            // SAME way regardless of (x,z) (i.e., not one of X_xmcol/Z_xmcol/XZ_idx):
            //   ∂eta^{xz}/∂beta_k = x_main_i[k]   (constant across x,z)
            //   ∂AIE_i/∂beta_k = x_main_i[k] * (phi_11 - phi_10 - phi_01 + phi_00)
            real colvector D_phi
            D_phi = phi_11 :- phi_10 :- phi_01 :+ phi_00
            Dq    = D_phi                              // d AIE_i / d(index shift)
            for (j = 1; j <= nx; j++) {
                J_beta_aie[j] = mean(D_phi :* x_main[., j])
            }
            // X_xmcol correction: ∂eta^{xz}/∂b_X = (x - X_i)
            //   So at (1,1): (1-X), at (1,0): (1-X), at (0,1): -X, at (0,0): -X
            //   ∂AIE_i/∂b_X = (1-X)*(phi_11 - phi_10) - X*(phi_01 - phi_00)
            // The generic above used X_val for x_main[X_xmcol] which is wrong; replace.
            J_beta_aie[X_xmcol] = mean( (1 :- X_val) :* (phi_11 - phi_10) ) ///
                                 - mean( X_val :* (phi_01 - phi_00) )
            // Z_xmcol correction: ∂eta^{xz}/∂b_Z = (z - Z_i)
            //   ∂AIE_i/∂b_Z = (1-Z)*(phi_11 - phi_01) - Z*(phi_10 - phi_00)
            J_beta_aie[Z_xmcol] = mean( (1 :- Z_val) :* (phi_11 - phi_01) ) ///
                                 - mean( Z_val :* (phi_10 - phi_00) )
            // XZ_idx correction: ∂eta^{xz}/∂b_XZ = (x*z - X*Z)
            //   at (1,1): (1 - X*Z), at (1,0): -X*Z, at (0,1): -X*Z, at (0,0): -X*Z
            //   ∂AIE_i/∂b_XZ = (1 - X*Z)*phi_11 + X*Z*(-phi_10 - phi_01 + phi_00)
            //                = phi_11 - X*Z*(phi_11 + phi_10 + phi_01 - phi_00)
            //   Wait, redo: ∂AIE_i/∂b_XZ
            //     = (1-X*Z)*phi_11 - (-X*Z)*phi_10 - (-X*Z)*phi_01 + (-X*Z)*phi_00
            //     = (1-X*Z)*phi_11 + X*Z*phi_10 + X*Z*phi_01 - X*Z*phi_00
            //     = phi_11 + X*Z*(phi_10 + phi_01 - phi_00 - phi_11)
            //     = phi_11 - X*Z*(phi_11 - phi_10 - phi_01 + phi_00)
            //     = phi_11 - X*Z*D_phi
            J_beta_aie[XZ_idx] = mean(phi_11) - mean(XZ_val :* D_phi)
        }
        else {
            _error("_apexogfit_apese_mt: unknown aie_case = " + strofreal(aie_case))
        }

        // Full two-step influence function and its (cluster) variance
        // (v0.21.0; through v0.20.1: J_beta_aie * V_beta_MT * J_beta_aie').
        J_pi      = mean(Dq :* GWD)                     // 1 x sum_k nw_k
        L         = psi2 * J_beta_aie' + PSI1 * J_pi'   // N x 1
        se_sq_aie = _apexogfit_meat(L, sort_idx, info, clustered)[1, 1]
        if (se_sq_aie < 0) {
            printf("{err}_apexogfit_apese_mt: J*V*J' = %g < 0 for AIE.\n", se_sq_aie)
            se_aie = .
        }
        else se_aie = sqrt(se_sq_aie)

        // Output: 1x1 matrices (the AIE is a single quantity, not per-endog)
        ape_vec = AIE_pt
        se_vec  = se_aie
        lo_vec  = AIE_pt - z_crit*se_aie
        hi_vec  = AIE_pt + z_crit*se_aie

        st_matrix(mname_ape, (AIE_pt))
        st_matrix(mname_se,  (se_aie))
        st_matrix(mname_lo,  (AIE_pt - z_crit*se_aie))
        st_matrix(mname_hi,  (AIE_pt + z_crit*se_aie))
        st_matrixrowstripe(mname_ape, ("", "AIE"))
        st_matrixrowstripe(mname_se,  ("", "AIE"))
        st_matrixrowstripe(mname_lo,  ("", "AIE"))
        st_matrixrowstripe(mname_hi,  ("", "AIE"))
        // 1x1 target covariance (for ereturn post)
        st_matrix("_apexogVTGT", J(1, 1, se_aie^2))
        return
    }

    // ====== APE branch (existing): per-endog loop ======
    ape_vec = J(n_endog, 1, .)
    se_vec  = J(n_endog, 1, .)
    lo_vec  = J(n_endog, 1, .)
    hi_vec  = J(n_endog, 1, .)
    Jall    = J(nx, n_endog, 0)      // stacked Jacobians for the cross-target cov
    L       = J(n_obs, n_endog, 0)   // per-target full influence functions (v0.21.0)

    for (t = 1; t <= n_endog; t++) {
        target_idx_t      = d_start_idx + t - 1
        chat_target_idx_t = chat_start_idx + t - 1     // (kept for clarity)
        b_d_t             = b_pr[target_idx_t]
        binary_d_t        = all(D[., t] :== 0 :| D[., t] :== 1)

        J_beta = J(1, nx, 0)

        if (binary_d_t) {
            // Discrete-change APE on d_target = d_t.
            eta_1    = eta + b_d_t :* (1 :- D[., t])
            eta_0    = eta - b_d_t :*  D[., t]
            Phi_1    = normal(eta_1)
            Phi_0    = normal(eta_0)
            phi_1    = normalden(eta_1)
            phi_0    = normalden(eta_0)
            diff_phi = phi_1 - phi_0
            ape_vec[t] = mean(Phi_1 - Phi_0)
            Dq         = diff_phi                      // d q_i / d(index shift)

            // Jacobian:
            //   d/d(b_const)        = mean(phi_1 - phi_0)
            //   d/d(b_X_v)          = mean(diff_phi * X_v)
            //   d/d(b_d_t)          = mean(phi_1)                 (asymmetric)
            //   d/d(b_d_k, k!=t)    = mean(diff_phi * d_k)
            //   d/d(b_chat_k, all k)= mean(diff_phi * chat0_k)
            //   d/d(b_cfprod_c)     = mean(diff_phi * cfprod_c)   (v0.20.1)
            for (j = 1; j <= n_exog; j++) {
                J_beta[j] = mean(diff_phi :* X_exog[., j])
            }
            for (k = 1; k <= n_endog; k++) {
                if (k == t) J_beta[d_start_idx + k - 1] = mean(phi_1)
                else        J_beta[d_start_idx + k - 1] = mean(diff_phi :* D[., k])
            }
            for (k = 1; k <= n_endog; k++) {
                J_beta[chat_start_idx + k - 1] = mean(diff_phi :* CHAT[., k])
            }
            // v0.20.1: cfinteract product coefficients (left at 0 before)
            for (cc = 1; cc <= n_cf; cc++) {
                J_beta[cf_start_idx + cc - 1] = mean(diff_phi :* CFPROD[., cc])
            }
            J_beta[cons_idx] = mean(diff_phi)
        }
        else {
            // Slope APE: mean(phi(eta)) * b_d_t
            mean_phi   = mean(phi_eta)
            ape_vec[t] = mean_phi * b_d_t
            mean_term  = -eta :* phi_eta
            Dq         = b_d_t :* mean_term            // d q_i / d(index shift)

            for (j = 1; j <= n_exog; j++) {
                J_beta[j] = b_d_t * mean(mean_term :* X_exog[., j])
            }
            for (k = 1; k <= n_endog; k++) {
                if (k == t)
                    J_beta[d_start_idx + k - 1] = b_d_t * mean(mean_term :* D[., k]) + mean_phi
                else
                    J_beta[d_start_idx + k - 1] = b_d_t * mean(mean_term :* D[., k])
            }
            for (k = 1; k <= n_endog; k++) {
                J_beta[chat_start_idx + k - 1] = b_d_t * mean(mean_term :* CHAT[., k])
            }
            // v0.20.1: cfinteract product coefficients (left at 0 before)
            for (cc = 1; cc <= n_cf; cc++) {
                J_beta[cf_start_idx + cc - 1] = b_d_t * mean(mean_term :* CFPROD[., cc])
            }
            J_beta[cons_idx] = b_d_t * mean(mean_term)
        }

        Jall[., t] = J_beta'

        // Full two-step influence function of target t (v0.21.0; through
        // v0.20.1 the SE was sqrt(J_beta * V_beta_MT * J_beta')).
        J_pi    = mean(Dq :* GWD)                      // 1 x sum_k nw_k
        L[., t] = psi2 * J_beta' + PSI1 * J_pi'
    }

    // Full cross-target covariance = (cluster) outer product of the stacked
    // influence functions.  Diagonal = se_vec^2; off-diagonals give the
    // covariance between the targets' APEs (used for ereturn post).
    Vtgt = _apexogfit_meat(L, sort_idx, info, clustered)
    Vtgt = (Vtgt + Vtgt') / 2
    st_matrix("_apexogVTGT", Vtgt)
    for (t = 1; t <= n_endog; t++) {
        se_sq_t = Vtgt[t, t]
        if (se_sq_t < 0) {
            printf("{err}_apexogfit_apese_mt: variance = %g < 0 for target %g.\n", se_sq_t, t)
            se_t = .
        }
        else se_t = sqrt(se_sq_t)
        se_vec[t] = se_t
        lo_vec[t] = ape_vec[t] - z_crit*se_t
        hi_vec[t] = ape_vec[t] + z_crit*se_t
    }

    // ---------- Return result vectors as Stata matrices ----------
    st_matrix(mname_ape, ape_vec)
    st_matrix(mname_se,  se_vec)
    st_matrix(mname_lo,  lo_vec)
    st_matrix(mname_hi,  hi_vec)
    st_matrixrowstripe(mname_ape, (J(n_endog, 1, ""), endogs'))
    st_matrixrowstripe(mname_se,  (J(n_endog, 1, ""), endogs'))
    st_matrixrowstripe(mname_lo,  (J(n_endog, 1, ""), endogs'))
    st_matrixrowstripe(mname_hi,  (J(n_endog, 1, ""), endogs'))
}

// ----------------------------------------------------------------------------
// _apexogfit_bp_ape (v0.19): APE + delta-method SE for the recursive bivariate
// probit engine (firststage(biprobit)).  The caller has just fit
//     biprobit (y = X d, offset(off)) (d = Z X)
// and stored e(b)/e(V) in _apexogPR_B/_apexogPR_V.  e(b) is ordered
//     [eq1: X_1..X_p, d, _cons | eq2: Z.., X.., _cons | athrho],
// and the first-equation block is located from the column stripe (so a depvar
// or endog with any name works).  With the gamma-offset baked into the outcome
// equation the target is the average-structural-function discrete change
//     APE = mean( Phi(eta_1) - Phi(eta_0) ),
//     eta = X b_X + d b_d + b_0 + off,  eta_1 = eta + b_d (1-d),  eta_0 = eta - b_d d.
// The outcome-equation error is standard normal, so Phi(.) IS the ASF: nothing
// is averaged out, and there is no control function.  Jacobian wrt the FULL
// parameter vector (zero for the treatment equation and athrho):
//     d/d b_X_v = mean((phi_1 - phi_0) X_v),   d/d b_d = mean(phi_1),
//     d/d b_0   = mean(phi_1 - phi_0).
// SE = sqrt(J V J') with V = biprobit's e(V): a plain MLE delta method, since
// d is modelled jointly (no generated regressor, no Murphy-Topel step).
// Also writes _apexogMPHI = mean(phi(eta)) (the gunits(pr) conversion) and
// _apexogSTILDE = 1: the outcome-equation error has unit variance, so the
// structural-latent and index scales coincide for this engine.
// ----------------------------------------------------------------------------
void _apexogfit_bp_ape(string scalar dname,     string scalar offvarname,
                       string scalar exognames,
                       string scalar mname_ape, string scalar mname_se,
                       string scalar mname_lo,  string scalar mname_hi,
                       real scalar   level)
{
    string rowvector evars
    string matrix    cstripe
    real matrix      X_exog, x1, V
    real colvector   d, off, eta, eta_1, eta_0, Phi_1, Phi_0, phi_1, phi_0, diff_phi
    real rowvector   b, b1, Jb
    real scalar      n_obs, n1, K, nxc, j, ape, se_sq, se, z_crit, d_idx, cons_idx
    string scalar    touse, eq1

    touse = st_local("touse")
    evars = (exognames == "") ? J(1, 0, "") : tokens(exognames)
    st_view(d,   ., dname,      touse)
    st_view(off, ., offvarname, touse)
    n_obs = rows(d)
    if (cols(evars) > 0) {
        st_view(X_exog, ., evars, touse)
        x1 = (X_exog, d, J(n_obs, 1, 1))
    }
    else x1 = (d, J(n_obs, 1, 1))
    nxc = cols(x1)

    b = st_matrix("_apexogPR_B")
    V = st_matrix("_apexogPR_V")
    K = cols(b)
    cstripe = st_matrixcolstripe("_apexogPR_B")
    eq1 = cstripe[1, 1]
    n1  = sum(cstripe[., 1] :== eq1)
    if (n1 != nxc) {
        _error("_apexogfit_bp_ape: biprobit outcome equation has " + strofreal(n1) +
               " columns; expected " + strofreal(nxc) + " ([exog d _cons]).")
    }
    b1       = b[|1 \ n1|]
    d_idx    = nxc - 1
    cons_idx = nxc

    eta = x1 * b1' + off
    st_numscalar("_apexogMPHI",   mean(normalden(eta)))
    st_numscalar("_apexogSTILDE", 1)

    eta_1    = eta + b1[d_idx] :* (1 :- d)
    eta_0    = eta - b1[d_idx] :* d
    Phi_1    = normal(eta_1)
    Phi_0    = normal(eta_0)
    phi_1    = normalden(eta_1)
    phi_0    = normalden(eta_0)
    diff_phi = phi_1 - phi_0
    ape      = mean(Phi_1 - Phi_0)

    Jb = J(1, K, 0)
    for (j = 1; j <= nxc - 2; j++) Jb[j] = mean(diff_phi :* X_exog[., j])
    Jb[d_idx]    = mean(phi_1)
    Jb[cons_idx] = mean(diff_phi)
    se_sq = (Jb * V * Jb')[1, 1]
    if (se_sq < 0) {
        printf("{err}_apexogfit_bp_ape: J*V*J' = %g < 0.\n", se_sq)
        se = .
    }
    else se = sqrt(se_sq)
    z_crit = invnormal(1 - (1 - level)/2)

    st_matrix(mname_ape, (ape))
    st_matrix(mname_se,  (se))
    st_matrix(mname_lo,  (ape - z_crit*se))
    st_matrix(mname_hi,  (ape + z_crit*se))
    st_matrixrowstripe(mname_ape, ("", dname))
    st_matrixrowstripe(mname_se,  ("", dname))
    st_matrixrowstripe(mname_lo,  ("", dname))
    st_matrixrowstripe(mname_hi,  ("", dname))
    st_matrix("_apexogVTGT", J(1, 1, se^2))
}

end
