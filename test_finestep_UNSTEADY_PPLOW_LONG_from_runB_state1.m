%% TEST_FINESTEP_UNSTEADY_PPLOW_LONG_FROM_RUNB_STATE1
% Extended-horizon follow-up to test_finestep_UNSTEADY_PPLOW_from_runB_state1.m
% (Issue 3, fluidPressurePenalty=100). That 6-step run completed cleanly
% -- physically plausible pressures throughout (104.6-7,024.5 Pa), only
% one recoverable solver failure (step 4, endothelium exitflag=0,
% survived via fallback), in ~68 minutes. Compare to the un-penalized
% baseline (rOuter-fix run), which got permanently stuck oscillating in
% its 2nd step for 85 passes over 18+ hours and never finished.
%
% That result is promising but was only a 6-step window. This run
% extends the horizon 5x (nFineSteps=20, ~26-30 real fine steps given the
% same tEnd-computation offset seen in the 6-step run) from the SAME
% starting state and SAME penalty (100, relax left at default 0.5), to
% see whether the improved behavior holds up over a longer horizon or
% whether the correction loop eventually falls back into the persistent-
% oscillation failure mode seen without the penalty.
%
% Uses the ORIGINAL (buggy rOuter=4e-6) starting state, unchanged, same
% as every other isolated test today. Writes to its own output/log files
% only -- does not touch or overwrite the original 6-step PPLOW result.

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

fprintf('Relaxation factor (unchanged): par.bodyFittedTractionCorrectionRelax = %.4f\n', ...
    par.bodyFittedTractionCorrectionRelax);
par.fluidPressurePenalty = 100;
fprintf('PROBE: par.fluidPressurePenalty = %g (LONG horizon, nFineSteps=%d)\n', ...
    par.fluidPressurePenalty, nFineSteps);

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
save('out_finestep_UNSTEADY_PPLOW_LONG.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_PPLOW_LONG.mat\n');

fprintf('\n=== RESULT (unsteady Stokes ON, both warm-start fixes + fluidPressurePenalty=100, LONG horizon) ===\n');
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
