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

if(cm_taxCO2_regiDiff = 11,

*** --- did the phi ratio survive the last postsolve? ---------------------------
*** FIRST thing in the file, deliberately: this reads p45_regiDiff_ratio exactly as the
*** previous iteration's postsolve left it, and exactly as core/presolve.gms has just
*** consumed it (core/loop.gms:54 runs before the module presolves, and turns
*** pm_taxCO2eq into pm_taxCO2eqSum - the parameter the equations see). Anywhere later
*** in this file the ratio has already been rebuilt from phi, which repairs the damage
*** and hides it. It must also come BEFORE the coupling call below, or the first call -
*** which makes phi non-uniform while the ratio is still legitimately uniform - would
*** trip it.
***
*** What it catches: postsolve.gms Step III.3 recomputes p45_regiDiff_ratio for every
*** cm_taxCO2_regiDiff except 0 and 3. Its cm_taxCO2_regiDiff = 11 branch was missing,
*** so mode 11 fell through to the else-branch and - because p45_regiDiff_endYr is 0 for
*** mode 11 - the "ratio = 1 from endYr" line set the ratio to 1 in EVERY year. phi was
*** computed, applied, and erased again before the solve ever saw it, while every
*** diagnostic in the gdx still reported a non-trivial phi. SSP2-PkBudg1000-PFMratio and
*** -PFMlevelC of the 2026-08-14/15 batch were lost to it and nothing flagged them.
  if((p45_pfmCallCount > 0) and (cm_pfmBindMode = 1) and
     (smax(regi, p45_regiDiff_phi(regi)) - smin(regi, p45_regiDiff_phi(regi)) > 1e-6),
    p45_pfmRatioSpread = smax((t,regi)$(t.val ge cm_startyear), p45_regiDiff_ratio(t,regi))
                       - smin((t,regi)$(t.val ge cm_startyear), p45_regiDiff_ratio(t,regi));
    if(p45_pfmRatioSpread < 1e-6,
      display p45_regiDiff_phi, p45_regiDiff_ratio, p45_pfmRatioSpread;
      abort "45_carbonprice: p45_regiDiff_ratio is uniform although phi is not - the PFM differentiation was erased between presolve and the solve. See postsolve.gms Step III.3.";
    );
  );

*** --- did the per-market markup survive the last postsolve? -------------------
*** The pm_taxemiMkt analogue of the p45_regiDiff_ratio check above, and here for the same
*** reason: this module writes pm_taxemiMkt in presolve, but 47_regipol's postsolve runs
*** AFTER this one and both zeroes it (cm_regiExoPrice, postsolve.gms:960/984) and rewrites
*** it (cm_emiMktTarget, postsolve.gms:384-436). datainput.gms aborts on the emiMktTarget
*** collision, and cm_regiExoPrice is documented as incompatible - but "documented as
*** incompatible" is exactly what defect 5 was, so this checks rather than assumes.
***
*** Reads pm_taxemiMkt as the previous postsolve left it, BEFORE the block further down
*** rebuilds it - anywhere later the damage is repaired and hidden. Fires only when this
*** module actually wrote a markup last time, so a legitimately zero markup (theta = 0,
*** or a region where the floor already is the higher sector) cannot trip it.
  if((cm_pfmSectorMarkup = 1) and (p45_pfmMarkupWritten > 1e-8),
    p45_pfmMarkupSeen = smax((t,regi,emiMkt)$(t.val ge cm_startyear), pm_taxemiMkt(t,regi,emiMkt));
    if(p45_pfmMarkupSeen < 1e-8,
      display p45_pfmMarkupWritten, p45_pfmMarkupSeen, pm_taxemiMkt;
      abort "45_carbonprice: the PFM per-market markup was written last presolve but is gone now - something erased pm_taxemiMkt between presolve and here. Check 47_regipol postsolve (cm_regiExoPrice, cm_emiMktTarget) - see ADR 0042 and COUPLING.md 11.4.";
    );
  );

*** --- the coupling call, skipped once phi has converged -----------------------
*** The loop is a fixed point: REMIND's energy system moves the ambition gaps, which
*** move phi, which moves the price, which moves the energy system. It has converged
*** when a further PFM call stops changing phi. After that the call is skipped for
*** the rest of the run - phi is FROZEN at its converged value, never reset to 1.
*** The R side reports "first call" as a huge delta, so the loop can never stop on
*** the first PFM iteration.
  if((pfmIter(iteration)) and (pm_pfmConverged = 0),

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
*** REMIND's solution before cm_startyear is FIXED to the reference run, so anything the
*** R side evaluates in that window cannot respond to the coupling. phi used to be pinned
*** at the first projection year (2025) while cm_startyear was 2030 - so phi was computed
*** entirely inside the frozen window and came back bit-identical on every call, delta
*** exactly 0, for the whole run. Handing the start year over lets the R side place its
*** tier year where the pathway can actually move.
    put "startYear: ", cm_startyear:0:0 /;
    putclose pfmcfg;

    Execute "Rscript -e 'library(pfm); pfm::iterativePFM()'";

    putclose runtime gyear(jnow):0:0 "-" gmonth(jnow):0:0 "-" gday(jnow):0:0 " " ghour(jnow):0:0 ":" gminute(jnow):0:0 ":" gsecond(jnow):0:0 ",GAMS," iteration.val:0;

*** Freshness FIRST, before anything is copied out of the gdx. The R side echoes the
*** iteration it was asked for; if that does not match, the file is a leftover from an
*** earlier call (or the call failed and wrote nothing) and no symbol in it may be
*** believed. Hoisted above the loads so every load below can be gated on it instead of
*** on a per-element "> 0" test - see the ETS block for why that distinction matters.
    Execute_Loadpoint 'p45_regiDiff_phi' p45_pfmIterSeen_aux = p45_pfmIterSeen;
    p45_pfmIterSeen = sum(regi, p45_pfmIterSeen_aux(regi)) / max(1, card(regi));
    p45_pfmFresh$(abs(p45_pfmIterSeen - iteration.val) <= 0.5) = 1;
    p45_pfmFresh$(abs(p45_pfmIterSeen - iteration.val) > 0.5) = 0;

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
*** The R side exports US$/tCO2; the anchor and pm_taxCO2eq are T$/GtC. Converting here once, so every later use is in model units.
      p45_pfmPriceBound(ttot,regi)$(p45_pfmPriceBound_aux(ttot,regi) > 0) =
        p45_pfmPriceBound_aux(ttot,regi) * sm_DptCO2_2_TDpGtC;
    );

*** Mode 3 carries its own price path rather than a share. Same "> 0" guard: a failed
*** R call leaves the previous path in place instead of zeroing the carbon price.
    if(cm_pfmBindMode = 3,
      Execute_Loadpoint 'p45_regiDiff_phi' p45_pfmMPPrice_aux = p45_pfmMPPrice;
*** Converted from US$/tCO2 to T$/GtC.
      p45_pfmMPPrice(ttot,regi)$(p45_pfmMPPrice_aux(ttot,regi) > 0) =
        p45_pfmMPPrice_aux(ttot,regi) * sm_DptCO2_2_TDpGtC;
    );

*** The PER-MARKET companions (ADR 0042). Loaded only when the markup is on, so a run
*** with cm_pfmSectorMarkup = 0 never touches them and stays bit-identical to the
*** pre-ADR behaviour.
***
*** Gated on p45_pfmFresh, NOT on a per-element "> 0". The economy-wide loads above use
*** "> 0" because their failure mode is "revert to phi = 1", i.e. silently UNCOUPLING a
*** coupled run, and no legitimate economy-wide phi is ever 0. Neither holds here:
*** a per-market phi of 0 is a legitimate value (theta -> 1 with that market's sector at
*** the bottom of the gap distribution), and "> 0" would silently discard it and fall
*** back to the floor - reading a maximally-constrained sector as an unconstrained one,
*** the wrong direction. The iteration stamp answers the question "> 0" was really
*** asking - is this gdx this call's? - without conflating it with the value being zero.
*** Same reasoning that retired the "> 0" test on p45_pfmDelta after the 2026-08-13
*** batch looped forever on a perfectly converged delta of exactly 0.
***
*** The aux is zeroed before each load so a gdx that is fresh but MISSING a companion
*** (an R-side export that failed for that symbol alone) cannot copy the previous
*** iteration's value in behind the stamp. Zero degrades to the floor, which is the
*** pre-ADR-0042 behaviour - the safe direction.
    if(cm_pfmSectorMarkup = 1,
      p45_pfmPhiMkt_aux(regi,emiMkt) = 0;
      p45_pfmLambdaMkt_aux(regi,emiMkt) = 0;
      Execute_Loadpoint 'p45_regiDiff_phi' p45_pfmPhiMkt_aux = p45_pfmPhiMkt;
      Execute_Loadpoint 'p45_regiDiff_phi' p45_pfmLambdaMkt_aux = p45_pfmLambdaMkt;
      if(p45_pfmFresh = 1,
        p45_pfmPhiMkt(regi,emiMkt) = p45_pfmPhiMkt_aux(regi,emiMkt);
*** Lambda keeps a "> 0" test ON TOP of the stamp, for a different reason: 0 is this
*** parameter's documented default ("the gap persists"), so an absent symbol and a
*** deliberate zero are indistinguishable, and falling back to the economy-wide rate is
*** the conservative reading of both.
        p45_pfmLambdaMkt(regi,emiMkt)$(p45_pfmLambdaMkt_aux(regi,emiMkt) > 0) =
          p45_pfmLambdaMkt_aux(regi,emiMkt);
      );
      if(cm_pfmBindMode = 2,
        p45_pfmPriceBoundMkt_aux(ttot,regi,emiMkt) = 0;
        Execute_Loadpoint 'p45_regiDiff_phi' p45_pfmPriceBoundMkt_aux = p45_pfmPriceBoundMkt;
        if(p45_pfmFresh = 1,
          p45_pfmPriceBoundMkt(ttot,regi,emiMkt) =
            p45_pfmPriceBoundMkt_aux(ttot,regi,emiMkt) * sm_DptCO2_2_TDpGtC;
        );
      );
      if(cm_pfmBindMode = 3,
        p45_pfmMPPriceMkt_aux(ttot,regi,emiMkt) = 0;
        Execute_Loadpoint 'p45_regiDiff_phi' p45_pfmMPPriceMkt_aux = p45_pfmMPPriceMkt;
        if(p45_pfmFresh = 1,
          p45_pfmMPPriceMkt(ttot,regi,emiMkt) =
            p45_pfmMPPriceMkt_aux(ttot,regi,emiMkt) * sm_DptCO2_2_TDpGtC;
        );
      );
      if(p45_pfmFresh = 0,
        display "45_carbonprice: STALE gdx - per-market companions not refreshed, markup keeps its previous values";
      );
      display p45_pfmPhiMkt, p45_pfmLambdaMkt;
    );

*** Convergence test. p45_pfmDelta is the largest change in ANY region's phi since
*** the previous PFM call - a max, not a mean, so one region still moving keeps the
*** loop open. A failed R call leaves the aux at its previous value; guarding on
*** "> 0" means a failure cannot be mistaken for convergence.
    Execute_Loadpoint 'p45_regiDiff_phi' p45_pfmDelta_aux = p45_pfmDelta;
    p45_pfmDelta = sum(regi, p45_pfmDelta_aux(regi)) / max(1, card(regi));
*** Freshness, not positivity. The old test required p45_pfmDelta > 0, meaning a
*** PERFECTLY converged loop - delta exactly 0 - was read as "not converged" and the
*** run kept calling R forever: the 2026-08-13 batch made 14 identical calls at ~2.6
*** minutes each and never stopped. The > 0 guard was there to stop a FAILED call
*** (which writes no gdx, so the previous values survive the load) from being mistaken
*** for convergence. That job now belongs to the iteration stamp the R side echoes
*** back: if it does not match the iteration we asked for, the gdx is a leftover and
*** nothing about it may be believed. The stamp is loaded once, at the top of this
*** block, and cached in p45_pfmFresh - the loads above are gated on the same test.
    if(p45_pfmFresh = 0,
      display "PFM coupling: the gdx is STALE - the R side did not answer this iteration. Keeping the previous phi.";
      display p45_pfmIterSeen, p45_pfmCallCount;
    else
      if((p45_pfmCallCount >= 2) and (p45_pfmDelta <= cm_pfmConvTol),
        pm_pfmConverged = 1;
        display "PFM coupling CONVERGED - phi frozen for the remainder of the run";
      );
    );
    display p45_regiDiff_phi, p45_pfmDelta, pm_pfmConverged, p45_pfmCallCount;
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
*** >>> Modes 2 and 3 below are MIRRORED in postsolve.gms Step IV.4. Change one, change
*** >>> the other. They cannot live here alone: core/loop.gms:54-55 runs core/presolve.gms
*** >>> - where pm_taxCO2eq becomes pm_taxCO2eqSum, the parameter the equations actually
*** >>> see - BEFORE the module presolves, so the solve consumes whatever pm_taxCO2eq was
*** >>> left holding at the END of the previous iteration, i.e. after postsolve. In a
*** >>> budget-iterating run postsolve Part IV rebuilds pm_taxCO2eq from the rescaled
*** >>> anchor and would otherwise throw the political layer away. Mode 1 needs no mirror:
*** >>> it IS ratio * anchor, so Step III.3 + Part IV reproduce it (and keep the
*** >>> path_gdx_ref interpolation, as every other cm_taxCO2_regiDiff mode does).

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
*** The first coupling call is c_pfmIter iterations in (default 15). Until then no
*** bound has been delivered and p45_pfmPriceBound is still at its GAMS default of
*** ZERO - and min(anchor, 0) is a ZERO CARBON PRICE in every region. Those early
*** iterations then solve a completely different model from the one the run is meant
*** to be, and under a held budget the Nash loop tears itself apart trying to close a
*** budget with no price: that is how SSP2-PkBudg1000-PFMlevelB died at iteration 5
*** with 42 infeasibilities, before the coupling had run even once. So cap only once
*** a real bound exists; before that the anchor stands unmodified.
  if(cm_pfmBindMode = 2,
    if(smax((t,regi)$(t.val ge cm_startyear), p45_pfmPriceBound(t,regi)) > 0,
      pm_taxCO2eq(t,regi)$(t.val ge cm_startyear) =
        min(p45_taxCO2eq_anchor(t), p45_pfmPriceBound(t,regi));
      p45_pfmBinds(t,regi)$(t.val ge cm_startyear) =
        1$(p45_pfmPriceBound(t,regi) < p45_taxCO2eq_anchor(t));
      display p45_pfmBinds;
    else
      pm_taxCO2eq(t,regi)$(t.val ge cm_startyear) = p45_taxCO2eq_regiDiff(t,regi);
      p45_pfmBinds(t,regi)$(t.val ge cm_startyear) = 0;
      display "45_carbonprice: bind mode 2 with no price bound yet - running on the anchor until the first PFM call";
    );
  );

*** Mode 3 (MILD PROGRESSION): the price is GENERATED by the political dynamics rather
*** than constraining the cost-optimal path. P(t+1) = P(t)(1 + lambda (S*-S)/S), seeded
*** with the observed current-policy price. There is no anchor and no budget in it, so
*** it cannot be infeasible - whatever emissions follow are the result.
*** The guard matters: if the R side never delivered a path, p45_pfmMPPrice is all
*** zero, and applying it would silently impose a ZERO carbon price everywhere - an
*** "uncoupled" run wearing a mild-progression label. Fail instead.
  if(cm_pfmBindMode = 3,
    if(smax((t,regi)$(t.val ge cm_startyear), p45_pfmMPPrice(t,regi)) <= 0,
      pm_pfmInfesCode = 3;
    else
      pm_taxCO2eq(t,regi)$(t.val ge cm_startyear) = p45_pfmMPPrice(t,regi);
      display p45_pfmMPPrice;
    );
  );

*** --- sector-differentiated delivery: the ETS markup (ADR 0042) ---------------
*** The three branches above set pm_taxCO2eq from the WORSE sector - the floor every
*** market pays. This adds back, PER MARKET, what that market's own sector could bear
*** beyond the floor.
***
*** Symmetric since 2026-08-17. The first version added the markup on ETS only and set
*** ES and other to zero, which capped the demand side at the Bulk price wherever BULK
*** was the worse sector - 14 of 48 countries on the deployed frontier, discarding up to
*** 0.30 of Diffuse phi. That reintroduced, on the other sector, exactly the information
*** loss ADR 0042 exists to remove. Now every market carries its own sector's price.
***
*** The invariant that matters: floor + markup(m) reproduces market m's own sector price
*** EXACTLY, so no sector is capped by the other. min() is arithmetic that keeps the
*** markup non-negative, not a modelling step that discards a sector.
***
*** Note it is NOT true that exactly one markup is positive. sectorRule = "min" takes the
*** worse SHARE and the slower SPEED, and those can come from different sectors, so the
*** floor is the most-constrained COMBINATION - belonging to neither sector and sometimes
*** strictly below both. Markups are still never negative. But pm_taxCO2eq is then a
*** price no market actually faces, and everything reading pm_taxCO2eqSum sees it: the
*** MAC curves, the land-use tax, the trade tariffs, the net-negative penalty. That is
*** the conservative direction, and it is deliberate, but it is a real distortion.
*** Pinned by test-exportFeasibilityBound.R "the min floor can sit BELOW both".
***
*** Why a markup and not a second price: pm_taxemiMkt is ADDITIVE. q21_taxrevGHG
*** charges pm_taxCO2eqSum on all CO2eq and q21_taxemiMkt then adds
*** pm_taxemiMkt(m) * vm_co2eqMkt(m) per market (21_tax/on/equations.gms), so the
*** effective price on market m is pm_taxCO2eq + pm_taxemiMkt(m). Keeping the floor in
*** pm_taxCO2eq means every OTHER consumer of pm_taxCO2eqSum keeps working untouched:
*** the MAC curves via p_priceCO2 (core/presolve.gms:246), the land-use CO2 tax
*** (q21_taxrevCO2luc), the trade tariffs (q21_tau_Import) and the net-negative penalty
*** (q21_taxrevNetNegEmi). Zeroing pm_taxCO2eq and carrying the whole price per market -
*** the way 47_regipol does - would silently zero ALL of those. The repair block that
*** rebuilds p_priceCO2 from pm_taxemiMkt (core/presolve.gms:255-270) is unreachable
*** here: it is gated on cm_emiMktTarget, which datainput.gms aborts on.
***
*** ETS ~ Bulk (electricity + industry), ES + other ~ Diffuse (buildings + transport);
*** the mapping is applied on the R side. It is good but not a bijection: sector2emiMkt
*** puts indst in BOTH ETS and ES, so industry's ES slice receives the Diffuse price.
*** Deliberate, and it errs toward LESS differentiation - i.e. toward the old min().
***
*** One-iteration lag, by design and pre-existing: core/presolve.gms runs before the
*** module presolves, so p_priceCO2forMAC sees the previous iteration's markup. The
*** SOLVE sees this one, because the equations read pm_taxemiMkt directly. The same
*** lag already applies to pm_taxCO2eqSum, and at a fixed point it vanishes.
  if(cm_pfmSectorMarkup = 1,
*** Each market's OWN closure rate, not p45_regiDiff_lambda. The economy-wide rate is the
*** one that survived sectorRule = "min", which takes the SLOWER of the two speeds -
*** Diffuse, 0.0770/yr against Bulk's 0.1023/yr (MODEL.md 4.3). Using it here would let
*** the faster sector converge on the anchor at the slower one's pace and understate the
*** markup by roughly a third. Modes 2 and 3 never had this problem: they receive
*** finished per-sector price paths from R. Mode 1 is the only branch that rebuilds the
*** path inside GAMS.
    if(cm_pfmBindMode = 1,
      p45_pfmPriceMkt(t,regi,emiMkt)$(t.val ge cm_startyear) =
        ( 1 - (1 - p45_pfmPhiMkt(regi,emiMkt))
              * rPower(1 - p45_pfmLambdaMkt(regi,emiMkt), t.val - p45_regiDiff_startYr(regi)) )
        * p45_taxCO2eq_anchor(t);
    elseif cm_pfmBindMode = 2,
      p45_pfmPriceMkt(t,regi,emiMkt)$(t.val ge cm_startyear) =
        min(p45_taxCO2eq_anchor(t), p45_pfmPriceBoundMkt(t,regi,emiMkt));
    elseif cm_pfmBindMode = 3,
      p45_pfmPriceMkt(t,regi,emiMkt)$(t.val ge cm_startyear) = p45_pfmMPPriceMkt(t,regi,emiMkt);
    );

*** Clamped at zero: pm_taxemiMkt is a markup, and a negative one would SUBSIDISE that
*** market's emissions relative to the rest of the economy. The market whose sector IS
*** the floor gets exactly zero, correctly. A missing or stale companion also lands here
*** as zero, which degrades to the old min() behaviour rather than to something invented.
    pm_taxemiMkt(t,regi,emiMkt)$(t.val ge cm_startyear) =
      max(p45_pfmPriceMkt(t,regi,emiMkt) - pm_taxCO2eq(t,regi), 0);

    p45_pfmMarkupShare_iter(iteration) =
      sum((t,regi,emiMkt)$(t.val ge cm_startyear), 1$(pm_taxemiMkt(t,regi,emiMkt) > 0))
      / max(1, sum((t,regi,emiMkt)$(t.val ge cm_startyear), 1));
    p45_pfmMarkupMean_iter(iteration) =
      sum((t,regi,emiMkt)$(t.val ge cm_startyear), pm_taxemiMkt(t,regi,emiMkt))
      / max(1, sum((t,regi,emiMkt)$(t.val ge cm_startyear), 1));

*** Remember what we wrote, so the next presolve can tell "erased" from "legitimately zero".
    p45_pfmMarkupWritten = smax((t,regi,emiMkt)$(t.val ge cm_startyear), pm_taxemiMkt(t,regi,emiMkt));

*** LIVENESS. The companions degrade silently by design: a missing or stale one lands as a
*** zero markup, i.e. the pre-ADR-0042 min() behaviour, which is the safe direction but is
*** indistinguishable from "the markup is switched on and working" in any output. So if the
*** per-market shares genuinely differ somewhere - the only case in which a markup is owed -
*** a markup of exactly zero everywhere means the companions never arrived.
*** Guarded on the shares, not on theta: at theta = 0 every share is 1, the spread is 0, and
*** a zero markup is correct.
    if(p45_pfmCallCount > 0,
      p45_pfmPhiMktSpread = sum(regi,
        smax(emiMkt, p45_pfmPhiMkt(regi,emiMkt)) - smin(emiMkt, p45_pfmPhiMkt(regi,emiMkt)));
      if((p45_pfmPhiMktSpread > 1e-6) and (p45_pfmMarkupWritten < 1e-8),
        display p45_pfmPhiMkt, p45_pfmPhiMktSpread, p45_pfmMarkupWritten;
        abort "45_carbonprice: the per-market shares differ but every markup is zero - the ADR 0042 companions did not reach GAMS, so this run is silently the old min() behaviour wearing a sector-differentiated label. Check that the R side wrote p45_pfmPhiMkt/p45_pfmPriceBoundMkt/p45_pfmMPPriceMkt and that the gdx is fresh (p45_pfmFresh).";
      );
    );

    display "45_carbonprice: per-market markups applied on top of the economy-wide floor";
    display p45_pfmMarkupShare_iter, p45_pfmMarkupWritten, pm_taxemiMkt;
  );

*** --- infeasibility detection ------------------------------------------------
*** Runs every iteration. Nothing here changes the solution; it classifies the run so
*** a blown-up result cannot be quoted as a feasible one. p45_pfmInfesCount requires
*** cm_pfmInfesPatience CONSECUTIVE flagged iterations, so a single noisy Nash
*** iteration does not condemn the run - and one clean iteration resets it.
  p45_pfmMaxPrice = smax((t,regi)$(t.val ge cm_startyear), pm_taxCO2eq(t,regi));
  p45_pfmRescaleHist(iteration) = p45_factorRescale_taxCO2_Funneled(iteration);
  pm_pfmInfesCode = 0;

*** (1) Price explosion. The SILENT failure: the solve succeeds and the numbers look
*** like results. This is the mode this run family has already hit.
  if(p45_pfmMaxPrice > cm_pfmMaxPrice,
    pm_pfmInfesCode = 1;
  );

*** (2) Budget iteration diverging: the rescale factor should approach 1 as the anchor
*** settles. Persistently far from 1 means the budget cannot be met by rescaling, which
*** under bind mode 2 is the expected consequence of a binding political cap.
  if((pm_pfmInfesCode = 0) and (ord(iteration) > 5) and
     (abs(p45_factorRescale_taxCO2_Funneled(iteration) - 1) > 0.5),
    pm_pfmInfesCode = 2;
  );

*** (3) A missing or zero price bound under bind mode 2. The R side already refuses to
*** run mode 2 without a bound, but this is the second line of defence INSIDE GAMS: a
*** bound of zero would cap every price at zero and the solve would succeed, reporting
*** a "politically infeasible" world that is really a plumbing failure. Checked against
*** a nominal floor rather than a real price path, because P_ref is not available here.
*** Only meaningful once a PFM call has actually happened: before the first one the
*** bound is legitimately zero (nothing has been delivered yet), and flagging that
*** would burn through cm_pfmInfesPatience on the uncoupled warm-up iterations and
*** condemn a run that has not yet done anything wrong.
  if((pm_pfmInfesCode = 0) and (cm_pfmBindMode = 2) and (p45_pfmCallCount > 0),
    if(smin((t,regi)$(t.val ge cm_startyear), p45_pfmPriceBound(t,regi)) <= 0,
      pm_pfmInfesCode = 3;
    );
  );

  if(pm_pfmInfesCode > 0,
    p45_pfmInfesCount = p45_pfmInfesCount + 1;
  else
    p45_pfmInfesCount = 0;
  );

  if(p45_pfmInfesCount >= cm_pfmInfesPatience,
    display "PFM COUPLING INFEASIBLE - see pm_pfmInfesCode (1 price explosion, 2 budget divergence, 3 bound below current policy)";
    display pm_pfmInfesCode, p45_pfmMaxPrice, p45_pfmInfesCount;
    execute_unload "pfm_infeasible.gdx", pm_pfmInfesCode, p45_pfmMaxPrice, p45_pfmInfesCount, pm_taxCO2eq, p45_pfmPriceBound, p45_regiDiff_phi;
  );
  display p45_pfmMaxPrice, pm_pfmInfesCode;

*** --- record this iteration ---------------------------------------------------
*** Written EVERY iteration, not only coupling ones, so the trace shows what happens
*** between PFM calls as well as at them - that is where the budget iteration moves
*** the anchor while phi is held fixed.
  p45_pfmPhi_iter(iteration,regi) = p45_regiDiff_phi(regi);
  p45_pfmDelta_iter(iteration) = p45_pfmDelta;
  p45_pfmMaxPrice_iter(iteration) = p45_pfmMaxPrice;
  p45_pfmInfes_iter(iteration) = pm_pfmInfesCode;
  p45_pfmConverged_iter(iteration) = pm_pfmConverged;
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
