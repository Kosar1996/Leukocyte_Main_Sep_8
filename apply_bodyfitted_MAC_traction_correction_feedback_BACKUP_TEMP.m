function [state, fluid, ok, stopReason, convergeInfo] = apply_bodyfitted_MAC_traction_correction_feedback( ...
    z, old, state, fluid, meshE, interfaceE, baseE, meshL, interfaceL, baseL, parL, par, opts)
%APPLY_BODYFITTED_MAC_TRACTION_CORRECTION_FEEDBACK
% Redesign of apply_bodyfitted_MAC_traction_correction.m to be a real
% feedback control loop,
%   "The solid fluid coupling should [be] strictly enforced in the ideal
%    case, and the correction should be chosen to satisfy this criteria
%    with a threshold of %mismatch or something. Dig into how the
%    correction works instead of just modifying its strength. Make sure
%    the correction considers both normal and tangential directions to
%    form a real feedback control loop."
%
% WHAT WAS WRONG WITH THE ORIGINAL: apply_bodyfitted_MAC_traction_
% correction.m runs a FIXED number of passes (maxBodyFittedTractionCorrections)
% at a FIXED relaxation factor, and never checks whether the interface
% traction mismatch has actually shrunk -- "convergence" was never
% measured, only assumed. Increasing the pass count and relaxation 10x
% (1 pass/0.05 -> 10 passes/0.5) helped one region (low-slope tangential
% mismatch 5.8->3.9 Pa) and did nothing for another (steep-slope
% tangential mismatch 23.7->23.8 Pa) -- exactly what you'd expect from an
% open-loop correction with no feedback signal telling it when to stop or
% whether it's actually working.
%
% WHAT THIS DOES INSTEAD: an actual closed loop --
%   1. MEASURE the current normal+tangential mismatch on both interfaces
%      (leukocyte|fluid and fluid|endothelium), using the same validated
%      methodology as check_interface_traction_mismatch_report.m /
%      compute_interface_traction_mismatch.m.
%   2. If the WORST of those four numbers (Leuko-en, Leuko-et, Endo-en,
%      Endo-et) is already below opts.mismatchThresholdPct, stop --
%      already converged, no correction needed.
%   3. Otherwise, do ONE partitioned correction pass (re-solve each solid
%      against the current fluid traction, under-relax, re-solve the
%      fluid on the corrected shape -- same mechanics as the original
%      function, now reusing the shared solve_finite_def_solid.m).
%   4. Go back to step 1 and re-measure. Repeat until converged or
%      opts.maxIterations passes have been used.
% The full mismatch-vs-iteration trajectory is returned in convergeInfo,
% so convergence behavior can be inspected/plotted directly instead of
% inferred from a single before/after comparison.
%
% Usage:
%   [state, fluid, ok, stopReason, convergeInfo] = ...
%       apply_bodyfitted_MAC_traction_correction_feedback( ...
%           z, old, state, fluid, meshE, interfaceE, baseE, ...
%           meshL, interfaceL, baseL, parL, par, opts)
%
% opts (all optional):
%   opts.mismatchThresholdPct  stop once the worst of the four %mismatch
%                              numbers (Leuko en/et, Endo en/et) is below
%                              this (default 10). NOT independently
%                              validated as "the right" tolerance for this
%                              physical problem -- a reasonable default,
%                              meant to be tuned.
%   opts.statForConvergence    'max' (default -- specified convergence
%                              must be judged by the largest-deviation point,
%                              not an average or a percentile; stated
%                              explicitly and repeated) or 'p90' (looser,
%                              robust to a few outlier sample points near
%                              mesh edges, opt-in only)
%   opts.maxIterations         safety cap on correction passes (default 10,
%                              matching the typical maxBodyFittedTraction
%                              Corrections usage in this codebase)
%   opts.relax                 fixed under-relaxation factor for the solid
%                              update (default par.bodyFittedTractionCorrectionRelax
%                              if present, else 1.0 -- same convention as
%                              the original function)
%   opts.mismatchNQuery        interface sample points used for the
%                              in-loop mismatch measurement (default 31 --
%                              coarser than check_interface_traction_
%                              mismatch_report.m's default 61, since this
%                              runs once per correction pass, not once per
%                              report)
%   opts.epsFrac                sampling offset fraction, same meaning/
%                              default (0.1) as compute_interface_
%                              traction_mismatch.m -- see check_interface_
%                              traction_epsFrac_sensitivity.m for why this
%                              matters: part of any given %mismatch
%                              reading is offset-dependent, not physical,
%                              so don't expect exact convergence to 0%.
%   opts.failMode              'warn' (default), 'error', or 'skip' -- what
%                              to do if opts.maxIterations is reached
%                              without converging. 'warn'/'skip' keep the
%                              best (lowest-mismatch) state found and
%                              continue; 'error' raises.
%   opts.verbose                true/false, print per-iteration mismatch
%                              (default true)
%
% Outputs:
%   state, fluid    same meaning as apply_bodyfitted_MAC_traction_correction.m
%   ok, stopReason  same meaning
%   convergeInfo    struct with:
%       .converged        true/false
%       .iterations       number of correction passes actually used
%       .history          Nx1 struct array, one row per measurement
%                          (iteration 0 = before any correction), each with
%                          .iter, .pctEnL, .pctEtL, .pctEnE, .pctEtE, .worstPct
%       .finalWorstPct    worst %mismatch at the returned state

if nargin < 13 || isempty(opts), opts = struct(); end
if ~isfield(opts,'mismatchThresholdPct') || isempty(opts.mismatchThresholdPct)
    opts.mismatchThresholdPct = 10;
end
if ~isfield(opts,'statForConvergence') || isempty(opts.statForConvergence)
    opts.statForConvergence = 'max';
end
if ~isfield(opts,'maxIterations') || isempty(opts.maxIterations)
    if isfield(par,'maxBodyFittedTractionCorrections') && isfinite(par.maxBodyFittedTractionCorrections) ...
            && par.maxBodyFittedTractionCorrections > 0
        opts.maxIterations = par.maxBodyFittedTractionCorrections;
    else
        opts.maxIterations = 10;
    end
end
if ~isfield(opts,'relax') || isempty(opts.relax)
    if isfield(par,'bodyFittedTractionCorrectionRelax') && isfinite(par.bodyFittedTractionCorrectionRelax)
        opts.relax = min(1.0, max(0.0, par.bodyFittedTractionCorrectionRelax));
    else
        opts.relax = 1.0;
    end
end
if ~isfield(opts,'mismatchNQuery') || isempty(opts.mismatchNQuery)
    opts.mismatchNQuery = 31;
end
if ~isfield(opts,'epsFrac') || isempty(opts.epsFrac)
    opts.epsFrac = 0.1;
end
if ~isfield(opts,'failMode') || isempty(opts.failMode)
    opts.failMode = 'warn';
end
if ~isfield(opts,'verbose') || isempty(opts.verbose)
    opts.verbose = true;
end
if ~isfield(opts,'useAitkenRelax') || isempty(opts.useAitkenRelax)
    % Fixed-fraction relax makes the correction chase a moving target: it
    % under-relax-steps the solid toward the CURRENT fluid.tractionE
    % snapshot, then the fluid gets re-solved on the new shape (lubrication
    % pressure is extremely sensitive to gap width, ~1/h^3), producing a new
    % traction the next pass is measured against -- so a step that reduced
    % mismatch against snapshot k can read as worse against snapshot k+1,
    % regardless of step size (confirmed empirically: relax=0.02..0.2 all
    % monotonically worsened from iteration 0). Aitken's Delta^2 relaxation
    % is the standard fix for exactly this loosely-coupled-partitioned-FSI
    % failure mode -- it derives the relaxation factor each pass from how
    % the correction's own residual actually changed, instead of a fixed
    % guess, and can shrink/grow/flip sign as needed. Default off --
    % opts.relax (fixed) remains the default behavior.
    opts.useAitkenRelax = false;
end
if ~isfield(opts,'aitkenRelaxBounds') || isempty(opts.aitkenRelaxBounds)
    opts.aitkenRelaxBounds = [-2, 2]; % clamp to avoid a wild single-pass jump
end
if ~isfield(opts,'debugRadialAxial') || isempty(opts.debugRadialAxial)
    % Diagnostic only -- prints, per correction pass, the norm of the
    % endothelium's RADIAL vs AXIAL displacement-update DOFs (uEcorr-uEold,
    % pre-relaxation), plus recovered sigma_rr/sigma_rz sampled at the
    % interface. Checks whether a frozen "EnE" (normal mismatch) is caused
    % by the correction's radial-direction update going to ~0 while the
    % axial update keeps moving, vs. a physical/measurement issue elsewhere.
    % Default off -- purely additive, no effect on returned state/fluid.
    opts.debugRadialAxial = false;
end

ok = true;
stopReason = '';

useRLoutInner = use_RLout_fluid_interface_for_solid_leukocyte(par);
hasL = ~useRLoutInner && ~isempty(meshL) && ~isempty(interfaceL) && ...
    isfield(state,'uL') && ~isempty(state.uL);

% Leukocyte-side parL for the mismatch measurement (Ge=GL, Ke=KL, etaE=etaL),
% same remap used everywhere else this comparison is made.
parLmismatch = par;
if isfield(par, 'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par, 'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par, 'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

mismatchOpts = struct('nQuery', opts.mismatchNQuery, 'epsFrac', opts.epsFrac, 'trimFrac', 0.05);

% par.solidAbsTol is tuned for the main monolithic solve's own force
% scale. Reused as-is here, it made solve_finite_def_solid's absolute
% residual check (resNorm < solidAbsTol) satisfied on the very first
% evaluation -- raw nodal forces at this micro-scale geometry are
% naturally ~1e-10 N, below the 1e-10 tolerance, regardless of how large
% the actual interface traction mismatch (in Pa, what this loop is
% supposed to be driving to zero) still is. That caused the correction to
% report "converged" (no update, uNew==uInitial) on iteration 1 every
% time, independent of opts.relax. Drop the inherited absolute tolerance
% here so solve_finite_def_solid falls back to its own tight built-in
% default (1e-13) instead, and actually takes Newton steps.
parCorrE = par;
if isfield(parCorrE, 'solidAbsTol'), parCorrE = rmfield(parCorrE, 'solidAbsTol'); end
parCorrL = parL;
if isfield(parCorrL, 'solidAbsTol'), parCorrL = rmfield(parCorrL, 'solidAbsTol'); end

history = struct('iter', {}, 'pctEnL', {}, 'pctEtL', {}, 'pctEnE', {}, 'pctEtE', {}, 'worstPct', {});
converged = false;
bestWorstPct = inf;
bestState = state;
bestFluid = fluid;

aitkenW = opts.relax;
aitkenPrevR = [];

ic = 0;
while true
    if ~isfield(fluid,'tractionE') || ~isfield(fluid,'tractionL')
        [fluid.tractionL, fluid.tractionE] = compute_bodyfitted_wall_traction(fluid.meshF, fluid, par);
    end

    [worstPct, row] = measure_mismatch(meshE, state.uE, old.uE, par, ...
        meshL, state.uL, old.uL, parLmismatch, par.dt, fluid, mismatchOpts, ...
        opts.statForConvergence, ic);
    history(end+1) = row; %#ok<AGROW>

    if opts.verbose
        fprintf(['   traction-correction feedback iter %d: pct(%s) EnL=%.2f%% EtL=%.2f%% ', ...
            'EnE=%.2f%% EtE=%.2f%%  worst=%.2f%% (target < %.2f%%)\n'], ...
            ic, opts.statForConvergence, row.pctEnL, row.pctEtL, row.pctEnE, row.pctEtE, ...
            worstPct, opts.mismatchThresholdPct);
    end

    if worstPct < bestWorstPct
        bestWorstPct = worstPct;
        bestState = state;
        bestFluid = fluid;
    end

    if worstPct < opts.mismatchThresholdPct
        converged = true;
        break;
    end

    if ic >= opts.maxIterations
        break;
    end

    ic = ic + 1;
    stateBeforeCorr = state;
    fluidBeforeCorr = fluid;
    try
        uEold = state.uE;
        uEcorr = solve_finite_def_solid(meshE, old.uE, fluid.tractionE, ...
            interfaceE, baseE, par.supportE, parCorrE, state.uE);

        if opts.debugRadialAxial
            dU = uEcorr - uEold;
            radialDU = dU(1:2:end);
            axialDU  = dU(2:2:end);
            stressDbg = recover_nodal_stress_axisym_viscoelastic(meshE, uEcorr, old.uE, par.dt, par);
            fprintf(['      [debug] pass %d: |dRadial|=%.4e |dAxial|=%.4e  ', ...
                'p90|sigma_rr|=%.4f p90|sigma_rz|=%.4f\n'], ...
                ic, norm(radialDU), norm(axialDU), ...
                prctile(abs(stressDbg.sigma_rr), 90), prctile(abs(stressDbg.sigma_rz), 90));
        end

        rE = uEcorr - uEold;

        if hasL
            uLold = state.uL;
            uLcorr = solve_finite_def_solid(meshL, old.uL, fluid.tractionL, ...
                interfaceL, baseL, par.supportL, parCorrL, state.uL);
            rL = uLcorr - uLold;
        else
            uLold = []; uLcorr = []; rL = [];
        end

        rCombined = [rE; rL];
        if opts.useAitkenRelax
            if isempty(aitkenPrevR) || numel(aitkenPrevR) ~= numel(rCombined)
                w = opts.relax; % first pass (or DOF count changed): fall back to the fixed initial guess
            else
                deltaR = rCombined - aitkenPrevR;
                denom = deltaR.' * deltaR;
                if denom > eps(class(denom)) * max(1, norm(rCombined))
                    w = -aitkenW * (aitkenPrevR.' * deltaR) / denom;
                    w = min(max(w, opts.aitkenRelaxBounds(1)), opts.aitkenRelaxBounds(2));
                else
                    w = aitkenW; % residual barely changed -- keep the previous factor
                end
            end
            aitkenW = w;
            aitkenPrevR = rCombined;
        else
            w = opts.relax;
        end

        state.uE = uEold + w * rE;
        [state.deltaE, state.UwE] = monolithic_interface_kinematics_value_only( ...
            meshE, state.uE, old.uE, interfaceE, z, par);

        if hasL
            state.uL = uLold + w * rL;
            [state.deltaL, state.UwL] = monolithic_interface_kinematics_value_only( ...
                meshL, state.uL, old.uL, interfaceL, z, par);
        end

        if opts.verbose && opts.useAitkenRelax
            fprintf('      [aitken] pass %d: w=%.4f\n', ic, w);
        end

        if useRLoutInner
            state.deltaL = par.RLout * ones(size(z));
            state.UwL = zeros(size(z));
        end

        state = attach_physical_solid_interface_fields( ...
            state, old, meshE, interfaceE, meshL, interfaceL, z, par);

        if any(state.deltaE(:) - state.deltaL(:) <= par.minGap)
            error('Body-fitted traction correction produced a gap below minGap.');
        end

        if ~isfield(par,'resolveFluidAfterTractionCorrection') || par.resolveFluidAfterTractionCorrection
            [fluid, ok, stopReason] = solve_selected_poststep_fluid(z, old, state, par);
            if ~ok
                state = bestState; fluid = bestFluid;
                convergeInfo = pack_converge_info(converged, ic, history, bestWorstPct);
                return;
            end
        end
    catch ME
        state = stateBeforeCorr;
        fluid = fluidBeforeCorr;
        stopReason = ME.message;

        failMode = lower(string(opts.failMode));
        if failMode == "warn" || failMode == "skip"
            ok = true;
            if failMode == "warn"
                warning(['apply_bodyfitted_MAC_traction_correction_feedback: pass %d failed, ', ...
                    'keeping best state found (worst pct mismatch=%.2f%%): %s'], ...
                    ic, bestWorstPct, ME.message);
            end
            state = bestState; fluid = bestFluid;
            convergeInfo = pack_converge_info(converged, ic, history, bestWorstPct);
            return;
        end
        ok = false;
        convergeInfo = pack_converge_info(converged, ic, history, bestWorstPct);
        return;
    end
end

if ~converged
    msg = sprintf(['apply_bodyfitted_MAC_traction_correction_feedback: did NOT converge within ', ...
        '%d passes. Best worst-case %%mismatch achieved = %.2f%% (target < %.2f%%). ', ...
        'Returning the best state found, not necessarily the last one.'], ...
        opts.maxIterations, bestWorstPct, opts.mismatchThresholdPct);
    failMode = lower(string(opts.failMode));
    if failMode == "error"
        error(msg); %#ok<SPERR>
    else
        warning(msg);
        stopReason = msg;
    end
end

state = bestState;
fluid = bestFluid;
convergeInfo = pack_converge_info(converged, ic, history, bestWorstPct);

end

function [worstPct, row] = measure_mismatch(meshE, uE, uEPrev, parE, meshL, uL, uLPrev, parL, dt, fluid, mismatchOpts, statName, iter)
[leuko, endo] = compute_interface_traction_mismatch( ...
    meshE, uE, uEPrev, parE, meshL, uL, uLPrev, parL, dt, fluid, mismatchOpts);

if strcmpi(statName, 'max')
    pctEnL = leuko.stats.pctMaxEn; pctEtL = leuko.stats.pctMaxEt;
    pctEnE = endo.stats.pctMaxEn;  pctEtE = endo.stats.pctMaxEt;
else
    pctEnL = leuko.stats.pctP90En; pctEtL = leuko.stats.pctP90Et;
    pctEnE = endo.stats.pctP90En;  pctEtE = endo.stats.pctP90Et;
end

worstPct = max([pctEnL, pctEtL, pctEnE, pctEtE]);

row = struct();
row.iter = iter;
row.pctEnL = pctEnL;
row.pctEtL = pctEtL;
row.pctEnE = pctEnE;
row.pctEtE = pctEtE;
row.worstPct = worstPct;
end

function convergeInfo = pack_converge_info(converged, iterations, history, finalWorstPct)
convergeInfo = struct();
convergeInfo.converged = converged;
convergeInfo.iterations = iterations;
convergeInfo.history = history;
convergeInfo.finalWorstPct = finalWorstPct;
end
