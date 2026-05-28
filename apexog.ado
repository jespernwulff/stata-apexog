*! apexog: APE/AIE bounds for cfprobit under plausibly exogenous IV (CHR 2012)
*! Version 0.12.0  2026-05-27
*! Targets the Average Partial Effect (APE) -- or, with target(c.X#i.d) or the
*! back-compat aie() shorthand, the Average Interaction Effect (AIE) in the
*! style of ginteff (Radean 2023) -- in a control-function probit model,
*! using the Conley/Hansen/Rossi (2012) UCI and LTZ frameworks.
*!
*! v0.12.0 changelog (vs v0.11.0):
*!   - RENAMED: plausexog_ape -> apexog.  Internal helpers _papefit,
*!     _papefit_apese_mt, _papeerr, and all _pape* matrices/scalars
*!     consistently renamed to _apexog* (clean prefix; matches the command).
*!     v0.11.0 source preserved in ./legacy/v0.11.0/ per the project's
*!     "copy not move" rule.  No behavioral change -- this version is
*!     functionally identical to v0.11.0.
*!
*! v0.11.0 changelog (vs v0.10.0) -- "Phase 2 of the Stata-native-syntax refactor":
*!   - target() now accepts Stata factor-variable interaction specs to declare
*!     an AIE target with EXPLICIT type prefixes that drive the formula choice:
*!         target(c.X#c.Z) -> cts x cts AIE   (cross-partial)
*!         target(c.X#i.d) -> cts x binary AIE (firstdiff in d)
*!         target(i.X#i.d) -> binary x binary AIE (2nd discrete difference)
*!     The prefixes OVERRIDE any value-based auto-detection (so the user can,
*!     e.g., force a discrete-change formulation for an indicator variable
*!     that happens to take values in {0,1} even if treated as continuous).
*!   - target() and aie() are mutually exclusive: target() supersedes aie() for
*!     anyone using the new fv-style syntax.  aie() is kept as a back-compat
*!     shorthand: aie(X Z) is internally translated to target(<X>#<Z>) with
*!     value-based prefix auto-detection (so existing scripts still work).
*!   - Display + ereturn additions: the "Target:" header shows the fv-style spec
*!     when AIE is requested via target(); new e(aie_spec) = "c.X#i.d" gives
*!     the parsed spec.
*!   - Internals retained from v0.10: the wrapper still creates an aie_XZ
*!     tempvar = X*Z and routes through the existing Mata AIE branches.  The
*!     "use a native fv interaction column in the probit" refactor is deferred
*!     to v0.12 (Phase 2b) -- it's purely an internal/representation change
*!     and does NOT affect the user-visible output.
*!
*! v0.10.0 changelog (vs v0.9.0) -- "Phase 1 of the Stata-native-syntax refactor":
*!   - New target() option for selecting which endogenous regressor's APE to
*!     report.  Default (omitted) = APE for every endog (current behavior).
*!     target(d1)         -> only the APE of d1 is displayed/ereturned
*!     target(d1 d2)      -> the APEs of d1 and d2
*!   The interaction-AIE workflow remains via the aie() option in Phase 1;
*!   target() will be extended to AIE specs (e.g. target(c.X#i.d)) in Phase 2,
*!   when the varlist gains factor-variable-interaction support.
*!   - New i./c. prefix support on endogenous-clause variables:
*!     (c.d = z)  -> treats d as continuous; defaults firststage(linear)
*!     (i.d = z)  -> treats d as binary;     defaults firststage(probit)
*!   The prefix only sets the DEFAULT first-stage type; firststage() still
*!   overrides.  Internally the prefix is stripped, so the main probit fits
*!   on bare `d' (binary-ness is then auto-detected from the values).
*!   - One-line "Target:" header in the output makes it explicit which
*!     functional is being reported.
*!   - Backup of v0.9.0 lives in ./legacy/v0.9.0/ (copy, not move).
*!
*! v0.9.0 changelog (vs v0.8.0):
*!   - New aie(<X> <Z>) option: switch the sensitivity target from the APE to
*!     the Average Interaction Effect for the interaction of X and Z.  The
*!     type of each interactor (continuous vs binary) is auto-detected and
*!     the appropriate AIE formula is used:
*!         cts x cts -- cross-partial: AIE = E[phi(eta)*(b_XZ - eta*a*b)]
*!         cts x bin -- firstdiff in the binary of the cts-derivative
*!         bin x bin -- 2nd discrete difference of Phi
*!     Multi-level factor interactors (3+ unique non-{0,1} values) are not
*!     supported in v0.9.
*!   - EITHER OR BOTH interactors may be endogenous -- this is the canonical
*!     cfprobit use case (interaction of an endog treatment with a moderator).
*!     apexog generates XZ = X*Z as a tempvar and adds it to the exog
*!     list for the main probit.  When an interactor is endogenous, the
*!     wrapper automatically EXCLUDES XZ = X*d from that endog's first-stage
*!     regression (otherwise the first stage would be circular: d on the LHS,
*!     X*d on the RHS).  A per-endog exog-mask matrix (_apexogEXOGMASK) carries
*!     this through to Mata so W_k matches what was actually fit.
*!   - Each AIE case has its own closed-form Jacobian wrt the full beta vector
*!     (with special corrections at the b_X, b_Z, b_XZ entries).  The Jacobian
*!     is passed through the same Murphy-Topel V_beta to give the delta-method
*!     SE, which then incorporates first-stage generated-regressor uncertainty.
*!   - Verified to machine precision against hand-coded analytic references in
*!     all four configurations: exog-by-exog cts (sanity 24), endog-by-exog cts
*!     (sanity 26), cts-by-binary (sanity 27), and binary-by-binary (sanity 28).
*!   - Display, ereturn, and graphing all branch cleanly on aie_active.  New
*!     ereturn locals: e(target)="aie", e(aie_vars)="<X> <Z>"; new scalars:
*!     e(lb_aie), e(ub_aie); baseline_ape_vec is 1x1 (the AIE).
*!
*! v0.8.0 changelog (vs v0.7.0):
*!   - Multi-clause syntax: (d1 = z1 z2) (d2 = z3 z4) so each endog can have
*!     its own instrument set.  Validated against cfprobit's per-clause syntax
*!     to machine precision on APEs and 0.1% on SEs.
*!   - APE bounds reported for EVERY endogenous regressor (not just the first).
*!     New ereturn matrices: e(baseline_ape_vec), e(baseline_se_vec), etc.
*!     Per-endog scalars: e(lb_<name>), e(ub_<name>) for each endog.
*!     UCI grid matrix has per-endog columns (ape_<name>, lo_<name>, hi_<name>).
*!     LTZ simulation matrix has one column per endog.
*!   - Mata refactor: per-endog W_k built from a Stata IV-mask matrix
*!     (J x n_iv) selecting active first-stage IVs.
*!
*! v0.7.0 changelog (vs v0.6.0):
*!   - Non-Gaussian LTZ priors via distribution(): normal, uniform, chi2, t,
*!     gamma, poisson, special.  Mirrors plausexog.ado's machinery.
*!   - firststage() now accepts either ONE type (applied to all endogs) OR
*!     a list (one per endog, in order).  Validated against cfprobit's
*!     per-clause first-stage syntax (machine-precision agreement).
*!   - Added SSC package metadata (apexog.pkg, stata.toc).
*!
*! v0.6.0 changelog (vs v0.5.0):
*!   - Multiple endogenous regressors.  Each gets its own first stage; all
*!     chat0s enter the main probit; MT correction accumulates first-stage
*!     uncertainty from every endog.  The APE target is the FIRST endog
*!     (reorder the varlist to change which one gets the SE).
*!     Validated: 2-endog APE matches cfprobit + margins to 1.5e-11.
*!
*! v0.5.0 changelog (vs v0.4.0):
*!   - vce(cluster <varname>): cluster-robust SEs throughout (first stage,
*!     main probit, AND the MT meat — sum of cluster-summed scores).
*!     Validated against Stata's native probit cluster SE.
*!   - graph option: line plot of UCI APE bounds vs gamma, or histogram
*!     of LTZ APE draws.  Optional graphoptions() forwarded to twoway.
*!   - apexog.sthlp written for SSC-ready documentation.
*!   - Investigated cfprobit's `margins, dydx(d)' for fprobit: confirmed that
*!     cfprobit uses the binary-d CF formula in its `predict, pr' even when
*!     d is fractional (see _remake_cfs.ado in Stata's source).  This is a
*!     cfprobit quirk; our APE uses the proper generalized residual.
*!
*! v0.4.0 changelog (vs v0.3.0):
*!   - Added firststage(fprobit) for fractional endogenous regressors in [0,1].
*!     The Bernoulli QMLE used by fracreg probit has the same score as probit,
*!     so chat0 and A_i are computed identically; only the Stata first-stage
*!     command differs.
*!   - Added firststage(poisson) for count endogenous regressors.  chat0 =
*!     d - exp(Wpi) (raw residual on the natural scale); A_i = exp(Wpi) = mu.
*!     The unified MT formula handles this with zero Mata changes.
*!
*! v0.3.0 changelog (vs v0.2.0):
*!   - Added firststage(probit) for binary endogenous regressors.
*!   - Auto-detect binary d -> discrete-change APE; else slope APE.
*!   - Unified MT Mata via per-obs weight A_i (probit: chat*(Wpi+chat); linear: 1).
*!
*! v0.2.0 changelog (vs v0.1.0):
*!   - Replaced fragile in-loop `gmm' call with `probit ... offset(...)' plus
*!     Murphy-Topel (two-step sandwich) SE computed in Mata.
*!
*! Supports:
*!   - method: uci, ltz
*!   - target: APE (default) for each endogenous regressor, OR AIE via
*!     aie(<X> <Z>) for a continuous-by-continuous interaction
*!   - LTZ priors: Gaussian via mu()/omega(), OR distribution(normal | uniform
*!     | chi2 | t | gamma | poisson | special) for non-Gaussian priors
*!   - 1+ endogenous regressors (APE reported for EACH one)
*!   - Multi-clause syntax: (d1 = z1 z2) (d2 = z3 z4) with per-endog IV sets
*!   - firststage: one type (linear | probit | fprobit | poisson) applied to
*!     all endogs, OR a list (one per endog) for per-endog typing
*!   - vce: robust (default), cluster <varname>
*!     (for bootstrap: wrap apexog in Stata's `bootstrap' prefix)
*!   - graph option for UCI grid plot and LTZ histogram (first endog)
*!
*! Known limitation:
*!   - cfprobit-style cfopts inside parens (e.g. `(d1 = z1, probit)`) are NOT
*!     supported because Stata's `syntax' command treats commas inside
*!     `anything()' as outer option separators.  Use the firststage(<t1>...)
*!     option instead for per-endog typing.

cap program drop apexog
cap program drop _apexogfit
cap program drop _apexogerr

program apexog, eclass
    version 18
    #delimit ;
    syntax anything(name=0 id="variable list")
    [if] [in]
    [,
        ESTimator(string)
        FIRststage(string)
        grid(integer 11)
        gmin(numlist)
        gmax(numlist)
        level(real 0.95)
        mu(numlist min=1)
        omega(numlist min=1)
        DISTribution(string)
        seed(numlist min=1 max=1)
        iterations(integer 1000)
        vce(string)
        noheader
        debug
        GRAph
        GRAPHOPTions(string asis)
        AIE(string)
        TARget(string)
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
    if `level' <= 0 | `level' >= 1 {
        _apexogerr "Confidence level must be in (0,1); you gave level(`level')."
    }

    local yvar `1'
    macro shift

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

        * Phase 1: parse i./c. prefixes on endog vars BEFORE fvexpand.  The
        * prefix only declares a TYPE (used to default first-stage type).
        * Strip it so downstream code (which expects bare varnames) is
        * unchanged; the type is stored in `endog_prefix_<n_clauses>_<j>'.
        local _stripped_endog ""
        local _pos = 0
        foreach _tok of local endog_`n_clauses' {
            local ++_pos
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
        local endog_`n_clauses' : list clean _stripped_endog

        * Expand factor-variable lists (now on the bare names)
        fvexpand `endog_`n_clauses''
        local endog_`n_clauses' `r(varlist)'
        fvexpand `iv_`n_clauses''
        local iv_`n_clauses' `r(varlist)'

        local all_endog `all_endog' `endog_`n_clauses''
        local all_iv    `all_iv'    `iv_`n_clauses''
    }
    if `n_clauses' < 1 {
        _apexogerr "No clauses found.  Use: METHOD yvar [exog] (endog = iv) [, opts]"
    }

    * Expand exog
    fvexpand `varlist1'
    local varlist1 `r(varlist)'

    * Build varlist2 (all endogs in order) and varlist_iv (UNIQUE IVs in order)
    local varlist2 : list uniq all_endog
    local varlist_iv : list uniq all_iv

    local count1   : word count `varlist1'
    local count2   : word count `varlist2'
    local count_iv : word count `varlist_iv'

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
    local countmin : word count `gmin'
    local countmax : word count `gmax'
    local count_mu : word count `mu'
    local count_om : word count `omega'

    *==========================================================================
    * (2) Validate options and v0.2 scope restrictions
    *==========================================================================
    if "`estimator'" == "" local estimator "cfprobit"
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
            if !inlist("`fst_k'", "linear", "probit", "fprobit", "poisson") {
                _apexogerr "Unknown firststage type '`fst_k''.  Valid: linear, probit, fprobit, poisson."
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
                if inlist("`_w'", "linear", "probit", "fprobit", "poisson") {
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
            * Validate both are in the regressor set.
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

    *--- Parse vce() ----------------------------------------------------------
    * Accept: vce(robust) [default], vce(cluster <varname>)
    if "`vce'" == "" local vce "robust"
    local vcefirst  : word 1 of `vce'
    local vcesecond : word 2 of `vce'
    if "`vcefirst'" == "robust" {
        local vcetype  = "robust"
        local clustvar = ""
    }
    else if "`vcefirst'" == "cluster" {
        if "`vcesecond'" == "" {
            _apexogerr "vce(cluster ...) requires a cluster variable name."
        }
        confirm variable `vcesecond'
        local vcetype  = "cluster"
        local clustvar = "`vcesecond'"
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
    marksample touse
    markout `touse' `yvar' `varlist1' `dvar' `varlist_iv'
    qui count if `touse'
    local N = r(N)
    if `N' == 0 _apexogerr "No observations remain after if/in and missing-data filtering."

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
        confirm numeric variable `aie_X'
        confirm numeric variable `aie_Z'

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

        * Create the interaction tempvar.  For binary-by-binary this is also the
        * `1.X#1.Z' indicator (since X,Z in {0,1}, X*Z = 1 iff both = 1).
        * The first-stage exog list (`exog_for_endog_<j>') is built below and
        * EXCLUDES this tempvar for any endog that is an interactor.
        tempvar aie_XZ
        qui gen double `aie_XZ' = `aie_X' * `aie_Z' if `touse'
        local varlist1 `varlist1' `aie_XZ'
        local count1 = `count1' + 1
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

    *==========================================================================
    * (4) Run the first stage ONCE (does not depend on gamma).  Save:
    *       - chat0     = control function (residual or generalized residual)
    *       - Wpi_hat   = first-stage linear predictor (used by Mata for the
    *                     per-obs MT weight A_i in nonlinear first stages)
    *       - A_per_obs = per-obs MT weight: 1 for linear, chat0*(Wpi+chat0)
    *                     for probit
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
    local j = 0
    foreach dj of local varlist2 {
        local ++j
        tempvar chat0_`j' Wpi_`j' A_`j'
        local fst = "`firststage_`j''"
        if "`fst'" == "linear" {
            qui regress `dj' `iv_for_endog_`j'' `exog_for_endog_`j'' if `touse', `vce_opt'
            qui predict double `chat0_`j'' if `touse', residuals
            qui predict double `Wpi_`j''   if `touse', xb
            qui gen double `A_`j'' = 1 if `touse'
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
            qui gen double `A_`j'' = `chat0_`j'' * (`Wpi_`j'' + `chat0_`j'') if `touse'
        }
        else if "`fst'" == "poisson" {
            qui poisson `dj' `iv_for_endog_`j'' `exog_for_endog_`j'' if `touse', `vce_opt'
            qui predict double `Wpi_`j'' if `touse', xb
            qui gen double `chat0_`j'' = `dj' - exp(`Wpi_`j'') if `touse'
            qui gen double `A_`j''     = exp(`Wpi_`j'') if `touse'
        }
        local chat0_list `chat0_list' `chat0_`j''
        local A_list     `A_list'     `A_`j''
    }

    * v0.6: all chat0s enter the main probit as regressors.  The MT correction
    * accumulates contributions from EVERY first stage.  The APE target is the
    * FIRST endogenous regressor (the rest are treated as controls).
    if `count2' > 1 {
        di as text "  [note] `count2' endogenous regressors detected.  All get a CF in the" /*
            */ " main probit; the MT correction accumulates uncertainty from all first stages." /*
            */ "  APE and SE are reported for the FIRST endog (`dvar')."
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

    *==========================================================================
    * (6) Baseline fit (gamma = 0)
    *==========================================================================
    qui replace `offvar' = 0 if `touse'
    _apexogfit, touse(`touse')                                                ///
        yvar(`yvar') endoglist(`varlist2') chat0list(`chat0_list')          ///
        offvar(`offvar') aweightlist(`A_list')                              ///
        exoglist(`varlist1') ivlist(`varlist_iv')                           ///
        clustvar(`clustvar') vcetype(`vcetype')                             ///
        aiex(`aie_X') aiez(`aie_Z') aiexz(`aie_XZ') aiecase(`aie_case')                         ///
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
    * (7) UCI: grid over [gmin_k, gmax_k] for each instrument k
    *==========================================================================
    if "`method'" == "uci" {
        local total_pts = `grid'^`count_iv'
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
        foreach tlabel of local target_labels {
            if `aie_active' local colnames "`colnames' aie_`aie_X'_`aie_Z' lo_`aie_X'_`aie_Z' hi_`aie_X'_`aie_Z'"
            else            local colnames "`colnames' ape_`tlabel' lo_`tlabel' hi_`tlabel'"
        }
        matrix colnames _apexogGRID = `colnames'

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
            foreach z of local varlist_iv {
                qui replace `offvar' = `offvar' + `g`k''*`z' if `touse'
                local ++k
            }

            _apexogfit, touse(`touse')                                            ///
                yvar(`yvar') endoglist(`varlist2') chat0list(`chat0_list')      ///
                offvar(`offvar') aweightlist(`A_list')                          ///
                exoglist(`varlist1') ivlist(`varlist_iv')                       ///
                clustvar(`clustvar') vcetype(`vcetype')                         ///
                aiex(`aie_X') aiez(`aie_Z') aiexz(`aie_XZ') aiecase(`aie_case')                     ///
                level(`level') `debug'
            matrix _apexogITER_APE = r(ape_vec)
            matrix _apexogITER_LO  = r(lo_vec)
            matrix _apexogITER_HI  = r(hi_vec)

            forvalues k = 1/`count_iv' {
                matrix _apexogGRID[`iter', `k'] = `g`k''
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
            di as text "Estimator: cfprobit  |  First stage: `firststage'  |  SE: Murphy-Topel  |  N = `N'"
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
            di as text "UCI union over gamma grid (`grid' point(s) per IV):"
            forvalues k = 1/`count_iv' {
                local zname : word `k' of `varlist_iv'
                di as text "    gamma_`zname' range  : [" %7.4f `gmin`k'' ", " %7.4f `gmax`k'' "]"
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

        ereturn clear
        ereturn scalar N             = `N'
        ereturn scalar level         = `level'
        *--- Optional graph: APE/AIE vs gamma over the UCI grid ----------------
        * Do this BEFORE the `ereturn matrix' so _apexogGRID still exists.
        if "`graph'" != "" & `count_iv' == 1 {
            preserve
            clear
            svmat _apexogGRID, names(col)
            sort g1
            local zname : word 1 of `varlist_iv'
            if `aie_active' {
                local _ape_lab "ape_`aie_X'_`aie_Z'"
                local _lo_lab  "lo_`aie_X'_`aie_Z'"
                local _hi_lab  "hi_`aie_X'_`aie_Z'"
                local _title   "UCI AIE bounds vs assumed gamma_`zname'"
                local _ytitle  "AIE (probability scale)"
            }
            else {
                local zfirst : word 1 of `varlist2'
                local _ape_lab "ape_`zfirst'"
                local _lo_lab  "lo_`zfirst'"
                local _hi_lab  "hi_`zfirst'"
                local _title   "UCI APE bounds vs assumed gamma_`zname'"
                local _ytitle  "APE (probability scale)"
            }
            if `"`graphoptions'"' == "" {
                local graphoptions ///
                    title("`_title'") ///
                    ytitle("`_ytitle'")                               ///
                    xtitle("gamma_`zname'")                            ///
                    legend(order(1 "Lower 95% CI" 2 "Upper 95% CI" 3 "Point estimate"))
            }
            twoway (line `_lo_lab' g1, lpattern(dash) lcolor(black))   ///
                   (line `_hi_lab' g1, lpattern(dash) lcolor(black))   ///
                   (line `_ape_lab' g1, lcolor(blue)),                 ///
                   yline(0, lpattern(dot) lcolor(gray))                  ///
                   `graphoptions'
            restore
        }
        else if "`graph'" != "" & `count_iv' > 1 {
            di as text "  [note] graph is only supported with a single instrument; skipping."
        }

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
                ereturn scalar lb_`dvar'     = `lo_union'
                ereturn scalar ub_`dvar'     = `hi_union'
            }
            * Per-endog scalars -- only for targets (Phase 1: skip non-targets)
            local j = 0
            foreach endogvar of local varlist2 {
                local ++j
                if !`target_endog_`j'' continue
                ereturn scalar lb_`endogvar' = `lo_union_`j''
                ereturn scalar ub_`endogvar' = `hi_union_`j''
            }
        }
        ereturn matrix grid               = _apexogGRID
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
        ereturn local cmdline        `"apexog `0_raw'"'
        ereturn local method         "uci"
        ereturn local estimator      "cfprobit"
        ereturn local firststage     "`firststage'"
        ereturn local vce            "Murphy-Topel"
        ereturn local depvar         "`yvar'"
        ereturn local endog          "`varlist2'"
        ereturn local exog           "`varlist1'"
        ereturn local instruments    "`varlist_iv'"
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
        local sim_colnames ""
        if `aie_active' {
            local sim_colnames "aie_`aie_X'_`aie_Z'"
        }
        else {
            foreach endogvar of local varlist2 {
                local sim_colnames "`sim_colnames' ape_`endogvar'"
            }
        }
        matrix colnames _apexogAPESIM = `sim_colnames'

        * Build the donor matrix for distribution(special) now that touse exists.
        if "`special_pending'" == "1" {
            tempname _apexogDONOR
            qui mkmat `special_donor_vars' if `touse', matrix(`_apexogDONOR') nomissing
            local _apexog_nrows = rowsof(`_apexogDONOR')
            forvalues k = 1/`count_iv' {
                local gammaCall`k' = "`_apexogDONOR'[ceil(runiform()*`_apexog_nrows'), `k']"
            }
        }

        forvalues iter = 1/`iterations' {
            * Draw gamma_k from the per-IV gammaCall expression
            * (gammaCall`k' was built above for both Gaussian and non-Gaussian
            * priors; for Gaussian it's "mu + sd*rnormal()", for distribution()
            * it's the appropriate Stata RNG expression).
            qui replace `offvar' = 0 if `touse'
            local k = 1
            foreach z of local varlist_iv {
                local g_draw = `gammaCall`k''
                qui replace `offvar' = `offvar' + `g_draw'*`z' if `touse'
                local ++k
            }

            _apexogfit, touse(`touse')                                            ///
                yvar(`yvar') endoglist(`varlist2') chat0list(`chat0_list')      ///
                offvar(`offvar') aweightlist(`A_list')                          ///
                exoglist(`varlist1') ivlist(`varlist_iv')                       ///
                clustvar(`clustvar') vcetype(`vcetype')                         ///
                aiex(`aie_X') aiez(`aie_Z') aiexz(`aie_XZ') aiecase(`aie_case')                     ///
                level(`level') `debug'

            matrix _apexogITER_APE = r(ape_vec)
            forvalues j = 1/`n_targets' {
                matrix _apexogAPESIM[`iter', `j'] = _apexogITER_APE[`j', 1]
            }

            if mod(`iter', 100) == 0 {
                di as text "  ... `iter' / `iterations' draws done"
            }
        }

        mata: st_local("_n_missing", strofreal(sum(rowmissing(st_matrix("_apexogAPESIM")))))
        if `_n_missing' > 0 {
            di as text "  [note] `_n_missing' of `iterations' LTZ draws failed (probit did not converge); excluded from percentile interval."
        }

        mata: _apexogsim_pct(`level')

        if "`header'" == "" {
            di as text _newline "{hline 78}"
            if `aie_active' {
                di as text "Plausibly-exogenous AIE bounds (CHR 2012, LTZ method, `ltz_mode' prior)"
            }
            else {
                di as text "Plausibly-exogenous APE bounds (CHR 2012, LTZ method, `ltz_mode' prior)"
            }
            di as text "Estimator: cfprobit  |  First stage: `firststage'  |  SE: Murphy-Topel  |  N = `N'"
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
            di as text "LTZ prior:"
            forvalues k = 1/`count_iv' {
                local zname : word `k' of `varlist_iv'
                if "`ltz_mode'" == "gaussian" {
                    di as text "    gamma_`zname' ~ N(" %7.4f `mu`k'' ", " %7.4f `omega`k'' ")"
                }
                else {
                    di as text "    gamma_`zname' ~ `distname'   draw expr: `gammaCall`k''"
                }
            }
            di as text "    iterations         = `iterations'"
            di as text " "
            local lvl = `level' * 100
            if `aie_active' {
                di as text "AIE distribution across draws (percentile intervals):"
                scalar _apexogLO = _apexogLO_PCT[1, 1]
                scalar _apexogHI = _apexogHI_PCT[1, 1]
                di as text "    AIE(`aie_X' ## `aie_Z') `lvl'% interval = [" %9.6f _apexogLO ", " %9.6f _apexogHI "]"
            }
            else {
                di as text "APE distribution across draws (percentile intervals):"
                local j = 0
                foreach endogvar of local varlist2 {
                    local ++j
                    if !`target_endog_`j'' continue
                    scalar _apexogLO = _apexogLO_PCT[`j', 1]
                    scalar _apexogHI = _apexogHI_PCT[`j', 1]
                    di as text "    APE(`endogvar') `lvl'% interval = [" %9.6f _apexogLO ", " %9.6f _apexogHI "]"
                }
            }
            di as text "{hline 78}"
        }

        ereturn clear
        ereturn scalar N             = `N'
        ereturn scalar level         = `level'
        ereturn scalar iterations    = `iterations'
        ereturn scalar baseline_ape  = `baseline_ape'
        ereturn scalar baseline_se   = `baseline_se'
        *--- Optional graph: histogram of the LTZ APE/AIE simulation draws ----
        * For multi-endog APE, plot the FIRST endog's distribution only.
        if "`graph'" != "" {
            preserve
            clear
            svmat _apexogAPESIM, names(col)
            if `aie_active' {
                local _histvar "aie_`aie_X'_`aie_Z'"
                local _title   "LTZ AIE distribution across `iterations' draws"
                local _subt    "Interaction: `aie_X' ## `aie_Z'"
                local _xt      "AIE (probability scale)"
            }
            else {
                local zfirst : word 1 of `varlist2'
                local _histvar "ape_`zfirst'"
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
                    note("Dashed lines: `=`level'*100'% percentile bounds; solid: baseline (gamma=0)")
            }
            twoway (histogram `_histvar', frequency)                   ///
                   ,                                                   ///
                   xline(`lo_pct',       lpattern(dash) lcolor(red))   ///
                   xline(`hi_pct',       lpattern(dash) lcolor(red))   ///
                   xline(`baseline_ape', lpattern(solid) lcolor(blue)) ///
                   `graphoptions'
            restore
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
                ereturn scalar lb_`dvar'     = `lo_pct'
                ereturn scalar ub_`dvar'     = `hi_pct'
            }
            * Per-endog scalars -- only for targets (Phase 1: skip non-targets)
            local j = 0
            foreach endogvar of local varlist2 {
                local ++j
                if !`target_endog_`j'' continue
                scalar _apexogLO_temp = _apexogLO_PCT[`j', 1]
                scalar _apexogHI_temp = _apexogHI_PCT[`j', 1]
                ereturn scalar lb_`endogvar' = _apexogLO_temp
                ereturn scalar ub_`endogvar' = _apexogHI_temp
            }
        }
        ereturn matrix ape_sim       = _apexogAPESIM
        ereturn matrix lo_pct_vec    = _apexogLO_PCT
        ereturn matrix hi_pct_vec    = _apexogHI_PCT
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
        ereturn local cmdline        `"apexog `0_raw'"'
        ereturn local method         "ltz"
        ereturn local estimator      "cfprobit"
        ereturn local firststage     "`firststage'"
        ereturn local vce            "Murphy-Topel"
        ereturn local depvar         "`yvar'"
        ereturn local endog          "`varlist2'"
        ereturn local exog           "`varlist1'"
        ereturn local instruments    "`varlist_iv'"
    }
end

*==============================================================================
* _apexogfit: workhorse.  At a given gamma-offset (already baked into `offvar'):
*   1. Fit `probit y X d chat0, offset(offvar) vce(robust)'
*   2. Pass results to Mata: compute Murphy-Topel-corrected V_beta, then
*      APE delta-method point estimate and SE.
*   3. Return ape, se, lo, hi via r().
*==============================================================================
program _apexogfit, rclass
    syntax,                                                                 ///
        touse(string)                                                       ///
        yvar(string) endoglist(string) chat0list(string) offvar(string)     ///
        aweightlist(string)                                                 ///
        exoglist(string) ivlist(string)                                     ///
        vcetype(string)                                                     ///
        level(real)                                                         ///
        [clustvar(string) aiex(string) aiez(string) aiexz(string)           ///
         aiecase(integer 0) debug]

    if "`debug'" != "" {
        qui sum `offvar' if `touse'
        di as text "[debug _apexogfit] offvar `offvar' mean=" %9.6f r(mean) " sd=" %9.6f r(sd)
        di as text "[debug _apexogfit] endoglist = `endoglist'   chat0list = `chat0list'"
    }

    * Build the vce string from vcetype + clustvar
    if "`vcetype'" == "cluster" local _vce = "vce(cluster `clustvar')"
    else                         local _vce = "vce(robust)"

    * Run probit with the offset.  Use `capture noisily quietly' pattern.
    capture {
        if "`debug'" != "" {
            probit `yvar' `exoglist' `endoglist' `chat0list' if `touse',    ///
                offset(`offvar') `_vce'
        }
        else {
            quietly probit `yvar' `exoglist' `endoglist' `chat0list' if `touse', ///
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
        exit
    }

    * Back-compat: also expose the FIRST endog's APE/SE as scalars for
    * existing callers that index r(ape), r(se), r(lo), r(hi).
    return scalar ape = _apexogAPE_vec[1, 1]
    return scalar se  = _apexogSE_vec[1, 1]
    return scalar lo  = _apexogLO_vec[1, 1]
    return scalar hi  = _apexogHI_vec[1, 1]
    return matrix ape_vec = _apexogAPE_vec
    return matrix se_vec  = _apexogSE_vec
    return matrix lo_vec  = _apexogLO_vec
    return matrix hi_vec  = _apexogHI_vec

    if "`debug'" != "" di as text "[debug _apexogfit] Mata returned: target APE=" %12.7f _apexogAPE_vec[1,1] " SE=" %12.7f _apexogSE_vec[1,1]
end

*==============================================================================
* _apexogerr: pretty error message
*==============================================================================
program _apexogerr
    di as error `"`0'"'
    exit 198
end

*==============================================================================
* Mata library
*==============================================================================
mata:

// ----------------------------------------------------------------------------
// _apexogsim_pct: per-endog percentile intervals for the LTZ ape simulation
// matrix.  _apexogAPESIM is N_iter x J (one column per endog).  Writes per-endog
// percentiles into Stata matrices _apexogLO_PCT (J x 1) and _apexogHI_PCT (J x 1).
// Also writes the FIRST endog's percentiles into locals lo_pct / hi_pct for
// back-compat with the single-endog display.
// ----------------------------------------------------------------------------
void _apexogsim_pct(real scalar level)
{
    real matrix      S, S_clean
    real colvector   col_clean, lo_vec, hi_vec
    real scalar      n_iter, J, j, n_valid, lo_pos, hi_pos, alpha
    S      = st_matrix("_apexogAPESIM")
    n_iter = rows(S)
    J      = cols(S)
    alpha  = 1 - level
    lo_vec = J(J, 1, .)
    hi_vec = J(J, 1, .)
    for (j = 1; j <= J; j++) {
        col_clean = select(S[., j], !missing(S[., j]))
        n_valid   = rows(col_clean)
        if (n_valid == 0) continue
        col_clean = sort(col_clean, 1)
        lo_pos    = max((1, round(alpha/2 * n_valid)))
        hi_pos    = min((n_valid, round((1 - alpha/2) * n_valid)))
        lo_vec[j] = col_clean[lo_pos]
        hi_vec[j] = col_clean[hi_pos]
    }
    st_matrix("_apexogLO_PCT", lo_vec)
    st_matrix("_apexogHI_PCT", hi_vec)
    // Back-compat: locals for the first endog
    st_local("lo_pct", strofreal(lo_vec[1], "%21.15g"))
    st_local("hi_pct", strofreal(hi_vec[1], "%21.15g"))
}

// ----------------------------------------------------------------------------
// _apexogfit_apese_mt: Murphy-Topel-corrected V_beta + APE delta-method SE.
//
// Unified for linear and probit first stages via a per-obs weight A_i:
//   linear:  A_i = 1                                  =>  w_dchat_i = -A_i = -1
//   probit:  A_i = chat0_i * (W_i'*pi + chat0_i)      =>  w_dchat_i = -A_i
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
//   Delta-method SE: sqrt( J_beta * V_beta_MT * J_beta' )
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
                       string scalar mname_ape,    string scalar mname_se,
                       string scalar mname_lo,     string scalar mname_hi,
                       real scalar   level)
{
    string rowvector evars, zvars, endogs, chats, Avars
    real matrix      X_exog, Z_iv, D, CHAT, AA, x_main, W
    real matrix      H_b, U_2, V_beta_MT, meat, s_2, s_total
    real matrix      M_pi_k, psi_pi_k, A1_k, J2pi_k
    real matrix      S_clust, cluster_ids, info
    real matrix      sort_idx, s_sorted, cluster_sorted
    real colvector   y, offset_vals, A_k, chat_k, eta, phi_eta, Phi_eta, Phim_eta
    real colvector   psi, psi_eta_psi, w_dchat_k, A2_k, ones_col, mean_term
    real colvector   eta_1, eta_0, Phi_1, Phi_0, phi_1, phi_0, diff_phi
    real rowvector   b_pr, J_beta
    real scalar      n_exog, n_iv, n_endog, n_obs, nx, nw, j, k, ci
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

    // x_main column order MUST match Stata's
    //   probit y EXOG d_1 ... d_J chat0_1 ... chat0_J
    // e(b) order: [X_1..X_p, d_1..d_J, chat0_1..chat0_J, _cons]
    if (n_exog > 0) x_main = (X_exog, D, CHAT, ones_col)
    else            x_main = (D, CHAT, ones_col)
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
    cons_idx        = nx
    d_target_idx    = d_start_idx                   // first endog = the target
    chat_target_idx = chat_start_idx                // its CF column

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

        A1_k = cross(x_main, psi_eta_psi :* w_dchat_k, W_k)   // nx x nw_k
        A2_k = cross(W_k,      psi          :* w_dchat_k)      // nw_k x 1
        J2pi_k = -b_chat_k :* A1_k
        J2pi_k[chat_start_idx + k - 1, .] = J2pi_k[chat_start_idx + k - 1, .] + A2_k'

        s_total = s_total + psi_pi_k * J2pi_k'       // N x nx (product)
    }

    // ---------- meat: cluster-robust or sandwich ----------
    if (clustvarname != "") {
        st_view(cluster_ids, ., clustvarname, touse)
        sort_idx       = order(cluster_ids, 1)
        cluster_sorted = cluster_ids[sort_idx, .]
        s_sorted       = s_total[sort_idx, .]
        info = panelsetup(cluster_sorted, 1)
        S_clust = J(rows(info), nx, 0)
        for (ci = 1; ci <= rows(info); ci++) {
            S_clust[ci, .] = colsum(s_sorted[|info[ci,1], 1 \ info[ci,2], nx|])
        }
        meat = cross(S_clust, S_clust)
    }
    else {
        meat = cross(s_total, s_total)
    }
    V_beta_MT = U_2 * meat * U_2

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
        // The XZ tempvar is always in evars (added by wrapper).
        real scalar X_xmcol, Z_xmcol, XZ_idx, k_search
        real colvector X_val, Z_val, XZ_val
        real rowvector J_beta_aie
        real scalar b_X_coef, b_Z_coef, b_XZ_coef, AIE_pt, se_aie, se_sq_aie
        X_xmcol = .; Z_xmcol = .; XZ_idx = .
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
        }
        if (X_xmcol == . | Z_xmcol == . | XZ_idx == .) {
            _error("_apexogfit_apese_mt: aie() variables not found in exog/endog list.")
        }
        XZ_val    = X_exog[., XZ_idx]
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

        // Delta-method SE (shared across cases)
        se_sq_aie = (J_beta_aie * V_beta_MT * J_beta_aie')[1, 1]
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
        return
    }

    // ====== APE branch (existing): per-endog loop ======
    ape_vec = J(n_endog, 1, .)
    se_vec  = J(n_endog, 1, .)
    lo_vec  = J(n_endog, 1, .)
    hi_vec  = J(n_endog, 1, .)

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

            // Jacobian:
            //   d/d(b_const)        = mean(phi_1 - phi_0)
            //   d/d(b_X_v)          = mean(diff_phi * X_v)
            //   d/d(b_d_t)          = mean(phi_1)                 (asymmetric)
            //   d/d(b_d_k, k!=t)    = mean(diff_phi * d_k)
            //   d/d(b_chat_k, all k)= mean(diff_phi * chat0_k)
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
            J_beta[cons_idx] = mean(diff_phi)
        }
        else {
            // Slope APE: mean(phi(eta)) * b_d_t
            mean_phi   = mean(phi_eta)
            ape_vec[t] = mean_phi * b_d_t
            mean_term  = -eta :* phi_eta

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
            J_beta[cons_idx] = b_d_t * mean(mean_term)
        }

        // Delta-method SE for target t
        se_sq_t = (J_beta * V_beta_MT * J_beta')[1,1]
        if (se_sq_t < 0) {
            printf("{err}_apexogfit_apese_mt: J*V*J' = %g < 0 for target %g.\n", se_sq_t, t)
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

end
