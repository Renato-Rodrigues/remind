*** |  (C) 2006-2024 Potsdam Institute for Climate Impact Research (PIK)
*** |  authors, and contributors see CITATION.cff file. This file is part
*** |  of REMIND and licensed under AGPL-3.0-or-later. Under Section 7 of
*** |  AGPL-3.0, you are granted additional permissions described in the
*** |  REMIND License Exception, version 1.0 (see LICENSE file).
*** |  Contact: remind@pik-potsdam.de
*** SOF ./modules/45_carbonprice/functionalForm/presolve.gms

*** ---------------------------------------------------------------------------
*** PFM political-feasibility coupling (cm_taxCO2_regiDiff = 11).
***
*** Runs INSIDE the Nash iteration loop, so the feasibility share phi can respond
*** to the energy system REMIND has just produced. This is why the coupling cannot
*** live in datainput.gms: $include is a COMPILE-TIME directive, so an .inc file is
*** baked in once and can never change between iterations. The interface therefore
*** mirrors EDGE-Transport's: shell out to R, then read the result back through a
*** gdx with Execute_Loadpoint.
***
*** The exchange each coupling iteration:
***   REMIND  --(fulldata.gdx: VRE share, electrification, fossil shares)-->  PFM
***   PFM     --(p45_regiDiff_phi.gdx: political feasibility share per region)--> REMIND
***
*** phi enters multiplicatively on the global anchor, so the budget iteration
*** (cm_iterative_target_adj = 5/7/9) keeps working: the anchor rescales until the
*** carbon budget is met and politics only redistributes WHERE abatement happens.
*** ---------------------------------------------------------------------------

file pfmcfg / "pfm-coupling-runtime.yml" /;

if(cm_taxCO2_regiDiff = 11,

*** --- the coupling call, skipped once phi has converged -----------------------
*** The loop is a fixed point: REMIND's energy system moves the ambition gaps, which
*** move phi, which moves the price, which moves the energy system. It has converged
*** when a further PFM call stops changing phi. After that the call is skipped for
*** the rest of the run - phi is FROZEN at its converged value, never reset to 1.
*** The R side reports "first call" as a huge delta, so the loop can never stop on
*** the first PFM iteration.
  if((pfmIter(iteration)) and (p45_pfmConverged = 0),

*** Track runtime, as EDGE-T does, so the coupling cost is visible in the log
    putclose runtime gyear(jnow):0:0 "-" gmonth(jnow):0:0 "-" gday(jnow):0:0 " " ghour(jnow):0:0 ":" gminute(jnow):0:0 ":" gsecond(jnow):0:0 ",iterativePFM," iteration.val:0;

*** Hand GAMS's own switches to R. The R process starts fresh with no arguments, so
*** without this the operator has to keep a duplicate copy of cm_pfmBindMode and
*** cm_pfmTheta in an .Rprofile - and a copy that disagrees produces a COMPLETE, WRONG
*** run. Writing them here makes the scenario config the single source of truth.
    put pfmcfg;
    put "# written by presolve.gms - do not edit, regenerated every coupling call" /;
    put "bindMode: ", cm_pfmBindMode:0:0 /;
    put "theta: ", cm_pfmTheta:0:6 /;
    put "convTol: ", cm_pfmConvTol:0:6 /;
    put "iteration: ", iteration.val:0:0 /;
    putclose pfmcfg;

    Execute "Rscript -e 'library(pfm); pfm::iterativePFM()'";

    putclose runtime gyear(jnow):0:0 "-" gmonth(jnow):0:0 "-" gday(jnow):0:0 " " ghour(jnow):0:0 ":" gminute(jnow):0:0 ":" gsecond(jnow):0:0 ",GAMS," iteration.val:0;

*** Load the updated feasibility shares. If the R side failed to produce the file
*** the previous iteration's phi is retained rather than silently reverting to 1 -
*** a failed coupling must not quietly turn into an uncoupled run.
    Execute_Loadpoint 'p45_regiDiff_phi' p45_regiDiff_phi_aux = p45_regiDiff_phi;
    p45_regiDiff_phi(regi)$(p45_regiDiff_phi_aux(regi) > 0) = p45_regiDiff_phi_aux(regi);
    p45_pfmCallCount = p45_pfmCallCount + 1;

*** Bind mode 2 also needs the ABSOLUTE politically feasible price. Same guard: a
*** zero means the R side did not supply it, so the previous value is kept rather
*** than capping the price at zero.
    if(cm_pfmBindMode = 2,
      Execute_Loadpoint 'p45_regiDiff_phi' p45_pfmPriceBound_aux = p45_pfmPriceBound;
      p45_pfmPriceBound(ttot,regi)$(p45_pfmPriceBound_aux(ttot,regi) > 0) = p45_pfmPriceBound_aux(ttot,regi);
    );

*** Convergence test. p45_pfmDelta is the largest change in ANY region's phi since
*** the previous PFM call - a max, not a mean, so one region still moving keeps the
*** loop open. A failed R call leaves the aux at its previous value; guarding on
*** "> 0" means a failure cannot be mistaken for convergence.
    Execute_Loadpoint 'p45_regiDiff_phi' p45_pfmDelta_aux = p45_pfmDelta;
    p45_pfmDelta = sum(regi, p45_pfmDelta_aux(regi)) / max(1, card(regi));
    if((p45_pfmDelta > 0) and (p45_pfmDelta <= cm_pfmConvTol),
      p45_pfmConverged = 1;
      display "PFM coupling CONVERGED - phi frozen for the remainder of the run";
    );
    display p45_regiDiff_phi, p45_pfmDelta, p45_pfmConverged, p45_pfmCallCount;
  );

*** --- apply phi, EVERY iteration ---------------------------------------------
*** Outside the pfmIter guard on purpose: the budget iteration keeps moving
*** p45_taxCO2eq_anchor, so the differentiated trajectory has to be rebuilt from the
*** current anchor even on iterations where phi was not recomputed. Leaving this
*** inside the guard would let pm_taxCO2eq go stale between coupling calls.
  p45_regiDiff_ratio(t,regi)$(t.val lt p45_regiDiff_startYr(regi)) = p45_regiDiff_phi(regi);
  p45_regiDiff_ratio(t,regi)$(t.val ge p45_regiDiff_startYr(regi)) =
    1 - (1 - p45_regiDiff_phi(regi))
        * rPower(1 - p45_regiDiff_lambda(regi), t.val - p45_regiDiff_startYr(regi));

  p45_taxCO2eq_regiDiff(t,regi) = p45_regiDiff_ratio(t,regi) * p45_taxCO2eq_anchor(t);

*** --- what phi actually binds (cm_pfmBindMode) --------------------------------
*** Mode 1 (RATIO): phi scales the anchor. The budget iteration can always meet the
*** budget by raising the anchor - and doing so raises the constrained region's price
*** too, so the political constraint weakens as ambition rises. Redistributive only.
  if(cm_pfmBindMode = 1,
    pm_taxCO2eq(t,regi)$(t.val ge cm_startyear) = p45_taxCO2eq_regiDiff(t,regi);
  );

*** Mode 2 (LEVEL): phi caps the ABSOLUTE price at phi * P_feasible. A region is held
*** to what its politics can deliver no matter how high the anchor goes, so the budget
*** may become unreachable. p45_pfmBinds records where the cap - rather than the
*** cost-optimal anchor - is setting the price; that flag is the raw material for both
*** headline B (how much budget politics costs) and headline C (who has to exceed
*** their limit, and by how much).
  if(cm_pfmBindMode = 2,
    pm_taxCO2eq(t,regi)$(t.val ge cm_startyear) =
      min(p45_taxCO2eq_anchor(t), p45_pfmPriceBound(t,regi));
    p45_pfmBinds(t,regi)$(t.val ge cm_startyear) =
      1$(p45_pfmPriceBound(t,regi) < p45_taxCO2eq_anchor(t));
    display p45_pfmBinds;
  );

*** --- infeasibility detection ------------------------------------------------
*** Runs every iteration. Nothing here changes the solution; it classifies the run so
*** a blown-up result cannot be quoted as a feasible one. p45_pfmInfesCount requires
*** cm_pfmInfesPatience CONSECUTIVE flagged iterations, so a single noisy Nash
*** iteration does not condemn the run - and one clean iteration resets it.
  p45_pfmMaxPrice = smax((t,regi)$(t.val ge cm_startyear), pm_taxCO2eq(t,regi));
  p45_pfmRescaleHist(iteration) = p45_factorRescale_taxCO2_Funneled(iteration);
  p45_pfmInfesCode = 0;

*** (1) Price explosion. The SILENT failure: the solve succeeds and the numbers look
*** like results. This is the mode this run family has already hit.
  if(p45_pfmMaxPrice > cm_pfmMaxPrice,
    p45_pfmInfesCode = 1;
  );

*** (2) Budget iteration diverging: the rescale factor should approach 1 as the anchor
*** settles. Persistently far from 1 means the budget cannot be met by rescaling, which
*** under bind mode 2 is the expected consequence of a binding political cap.
  if((p45_pfmInfesCode = 0) and (ord(iteration) > 5) and
     (abs(p45_factorRescale_taxCO2_Funneled(iteration) - 1) > 0.5),
    p45_pfmInfesCode = 2;
  );

*** (3) A missing or zero price bound under bind mode 2. The R side already refuses to
*** run mode 2 without a bound, but this is the second line of defence INSIDE GAMS: a
*** bound of zero would cap every price at zero and the solve would succeed, reporting
*** a "politically infeasible" world that is really a plumbing failure. Checked against
*** a nominal floor rather than a real price path, because P_ref is not available here.
  if((p45_pfmInfesCode = 0) and (cm_pfmBindMode = 2),
    if(smin((t,regi)$(t.val ge cm_startyear), p45_pfmPriceBound(t,regi)) <= 0,
      p45_pfmInfesCode = 3;
    );
  );

  if(p45_pfmInfesCode > 0,
    p45_pfmInfesCount = p45_pfmInfesCount + 1;
  else
    p45_pfmInfesCount = 0;
  );

  if(p45_pfmInfesCount >= cm_pfmInfesPatience,
    display "PFM COUPLING INFEASIBLE - see p45_pfmInfesCode (1 price explosion, 2 budget divergence, 3 bound below current policy)";
    display p45_pfmInfesCode, p45_pfmMaxPrice, p45_pfmInfesCount;
    execute_unload "pfm_infeasible.gdx", p45_pfmInfesCode, p45_pfmMaxPrice, p45_pfmInfesCount, pm_taxCO2eq, p45_pfmPriceBound, p45_regiDiff_phi;
  );
  display p45_pfmMaxPrice, p45_pfmInfesCode;

*** --- record this iteration ---------------------------------------------------
*** Written EVERY iteration, not only coupling ones, so the trace shows what happens
*** between PFM calls as well as at them - that is where the budget iteration moves
*** the anchor while phi is held fixed.
  p45_pfmPhi_iter(iteration,regi) = p45_regiDiff_phi(regi);
  p45_pfmDelta_iter(iteration) = p45_pfmDelta;
  p45_pfmMaxPrice_iter(iteration) = p45_pfmMaxPrice;
  p45_pfmInfes_iter(iteration) = p45_pfmInfesCode;
  p45_pfmConverged_iter(iteration) = p45_pfmConverged;
  p45_pfmAnchor_iter(iteration) = sum(t$(t.val eq 2050), p45_taxCO2eq_anchor(t));
  p45_pfmPriceMean_iter(iteration,regi) =
    sum(t$(t.val ge cm_startyear), pm_taxCO2eq(t,regi))
    / max(1, sum(t$(t.val ge cm_startyear), 1));
  if(cm_pfmBindMode = 2,
    p45_pfmBoundMean_iter(iteration,regi) =
      sum(t$(t.val ge cm_startyear), p45_pfmPriceBound(t,regi))
      / max(1, sum(t$(t.val ge cm_startyear), 1));
    p45_pfmBindShare_iter(iteration) =
      sum((t,regi)$(t.val ge cm_startyear), p45_pfmBinds(t,regi))
      / max(1, sum((t,regi)$(t.val ge cm_startyear), 1));
  );
);

*** EOF ./modules/45_carbonprice/functionalForm/presolve.gms
