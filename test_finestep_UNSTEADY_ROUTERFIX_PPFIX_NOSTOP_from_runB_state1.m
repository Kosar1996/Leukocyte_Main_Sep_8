%% TEST_FINESTEP_UNSTEADY_ROUTERFIX_PPFIX_NOSTOP_FROM_RUNB_STATE1
% Follow-up to test_finestep_UNSTEADY_ROUTERFIX_PPFIX_from_runB_state1.m
% (Sep 2), which combined the rOuter fix (par.global1DOuterRadius=15e-6)
% with the Issue 3 candidate fix (par.fluidPressurePenalty=100) and got
% past the exact point that permanently stuck the un-penalized version
% for 18+ hours. That run stopped after only 2 steps, but NOT because of
% any failure -- it hit the code's built-in stopAtMinGap safety trigger
% (par.gapStopFactor=2 x par.minGap=1e-9 m = 2 nm; actual hmin=1.29 nm),
% a numerical safety stop, not an independently-verified physical contact
% criterion. Sep 2's conclusion was too generous in calling that the
% "intended physical endpoint" -- this test corrects that by disabling
% the early stop (par.stopAtMinGap=false) to see what actually happens if
% the run is forced to keep going as the gap pushes toward and at the
% true 1 nm minGap floor: does the combined fix keep behaving well, or
% was trouble just deferred by the stop trigger firing early?
%
% Uses the ORIGINAL starting state (out_pure2dmac_dtlarge_7steps.mat,
% state1), same as every test this week. Writes to its own output/log
% files only -- does not touch the original ROUTERFIX_PPFIX run's files.

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
par.fluidPressurePenalty = 100;
fprintf('PROBE: par.fluidPressurePenalty = %g\n', par.fluidPressurePenalty);

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
save('out_finestep_UNSTEADY_ROUTERFIX_PPFIX_NOSTOP.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_ROUTERFIX_PPFIX_NOSTOP.mat\n');

fprintf('\n=== RESULT (rOuter fix + fluidPressurePenalty=100 + stopAtMinGap DISABLED) ===\n');
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
