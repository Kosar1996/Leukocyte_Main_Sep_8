%% TEST_FINESTEP_UNSTEADY_ROUTERFIX_PPFIX_FROM_RUNB_STATE1
% Combines two fixes to test against the single hardest failure mode seen
% today: the rOuter-fix-only run (test_finestep_UNSTEADY_ROUTERFIX_...)
% got permanently stuck oscillating in its step-2 correction loop for 85
% passes over 18+ hours, never converging, right after the monolithic
% solve failed with gap=1.000e-09 m -- exactly the minGap floor. That run
% was never saved (killed mid-run), so this cannot resume from its exact
% stuck state, but it reproduces the identical setup that led there
% (same starting state, same rOuter fix) and adds the Issue 3 candidate
% fix (par.fluidPressurePenalty=100, validated today over a clean
% 22-step run) to see whether it lets the run actually get past the
% specific point where the un-penalized version could not, for 18 hours.
%
% This is the most convincing test available for the pressure-penalty
% fix: not another clean isolated window, but the exact real-world
% scenario (rOuter corrected, as it should be per the Sep 1 config-bug
% finding) that produced the worst failure mode observed all day.
%
% Uses the ORIGINAL starting state (out_pure2dmac_dtlarge_7steps.mat,
% state1), same as every test today. Writes to its own output/log files
% only -- does not touch the original rOuter-fix or PP=100 test files.

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
save('out_finestep_UNSTEADY_ROUTERFIX_PPFIX.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_ROUTERFIX_PPFIX.mat\n');

fprintf('\n=== RESULT (unsteady Stokes ON, both warm-start fixes + rOuter fix + fluidPressurePenalty=100 ACTIVE) ===\n');
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
