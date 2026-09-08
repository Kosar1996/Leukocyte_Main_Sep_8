%% TEST_FINESTEP_UNSTEADY_ROUTERFIX_PP10_NOSTOP_FROM_RUNB_STATE1
% Direct counterpart to test_finestep_UNSTEADY_ROUTERFIX_PPFIX_NOSTOP_
% from_runB_state1.m (Sep 3), which combined the rOuter fix with
% fluidPressurePenalty=100 and stopAtMinGap=false, and completed a full
% 22-step run cleanly, pressure decaying to numerical zero, gap pinned
% at the true 1nm floor from step 3 onward.
%
% This swaps PP=100 for PP=10 -- everything else identical (same rOuter
% fix, same stopAtMinGap=false, same starting state, same 20-step
% target) -- motivated by the Sep 3 condition-number sweep, which found
% PP=10 gives an 8x BETTER condition number than PP=100 at the exact
% near-floor gap state, in a single-solve snapshot test. This is the
% real test: does that conditioning advantage translate into equally
% good (or better) real time-marching behavior, isolating PP value as
% the only variable against the already-completed PP=100 result.
%
% Uses the ORIGINAL starting state (out_pure2dmac_dtlarge_7steps.mat,
% state1), same as every test this week. Writes to its own output/log
% files only -- does not touch the PP=100 NOSTOP run's files.

clc;

S = load('out_pure2dmac_dtlarge_7steps.mat');
out = S.out;

fprintf('Run B step 1: t=%.6e s, dt used=%.6e s\n', out.t(1), out.dtHist(1));
fprintf('Run B step 2 (actual, coarse dt, steady): t=%.6e s, max|p|=%.4f Pa\n', ...
    out.t(2), max(abs(out.fluidHist{2}.p)));

state1 = out.stateHist{1};
par = out.par;
par.dt = 3e-4;
par.useUnsteadyStokes = true;
par.rho = 1000;
nFineSteps = 20;
par.tEnd = state1.t + nFineSteps * par.dt;

fprintf('BEFORE FIX: par.global1DOuterRadius = %.4e m\n', par.global1DOuterRadius);
par.global1DOuterRadius = 15e-6;
fprintf('AFTER FIX:  par.global1DOuterRadius = %.4e m\n', par.global1DOuterRadius);

fprintf('Relaxation factor (unchanged): par.bodyFittedTractionCorrectionRelax = %.4f\n', ...
    par.bodyFittedTractionCorrectionRelax);
par.fluidPressurePenalty = 10;
fprintf('PROBE: par.fluidPressurePenalty = %g (vs 100 in the earlier comparable run)\n', ...
    par.fluidPressurePenalty);

fprintf('BEFORE FIX: par.stopAtMinGap = %d (gapStopFactor=%g, minGap=%.3e m)\n', ...
    par.stopAtMinGap, par.gapStopFactor, par.minGap);
par.stopAtMinGap = false;
fprintf('AFTER FIX:  par.stopAtMinGap = %d -- run will continue past the 2nm safety trigger\n', ...
    par.stopAtMinGap);

restart = struct();
restart.state = state1;
restart.par = par;
restart.meshE = out.meshE;
restart.interfaceE = out.interfaceE;
restart.meshL = out.meshL;
restart.interfaceL = out.interfaceL;
restart.z = out.z;
restart.PHist = out.PHist;
restart.urCHist = out.urCHist;
restart.uzCHist = out.uzCHist;
restart.speedCHist = out.speedCHist;
restart.RPHist = out.RPHist;
restart.ZPHist = out.ZPHist;

cfg = struct();
cfg.geometry.endotheliumPrestressFile = 'solid_endo_P300_wide.mat';
cfg.geometry.leukocytePrestressFile = 'solid_leu_P600.mat';

outFine = softlube_run_case_global_coupled(cfg, restart);

% Save FIRST, before any post-processing that could crash and lose the run.
save('out_finestep_UNSTEADY_ROUTERFIX_PP10_NOSTOP.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_ROUTERFIX_PP10_NOSTOP.mat\n');

fprintf('\n=== RESULT (rOuter fix + fluidPressurePenalty=10 + stopAtMinGap DISABLED) ===\n');
for k = 1:outFine.stopStep
    if k <= numel(outFine.fluidHist) && isstruct(outFine.fluidHist{k}) && isfield(outFine.fluidHist{k}, 'p')
        fprintf('Fine step %d: t=%.6e s, max|p|=%.4f Pa\n', ...
            k, outFine.t(k), max(abs(outFine.fluidHist{k}.p)));
    else
        fprintf('Fine step %d: t=%.6e s, INVALID fluidHist entry (class=%s)\n', ...
            k, outFine.t(k), class(outFine.fluidHist{k}));
    end
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
