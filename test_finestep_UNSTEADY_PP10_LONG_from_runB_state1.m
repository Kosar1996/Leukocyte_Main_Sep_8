%% TEST_FINESTEP_UNSTEADY_PP10_LONG_FROM_RUNB_STATE1
% Follow-up to the Sep 3 condition-number sweep, which found PP=100 is
% NOT the best-conditioned value in its own range: at the exact near-
% floor gap state (1.000e-09 m), PP=10 gave an 8x BETTER condition
% number than PP=100 (5.37e7 vs 4.17e8, both vs 3.65e9 for PP=0). That
% was a single-state, single-solve snapshot -- this is the real test:
% does PP=10 hold up as well as PP=100 did over a genuine time-marching
% run (22 steps, same starting state, same horizon as the PP=100 long
% run), or does its condition-number advantage not translate into
% comparable real-run stability? If it does hold up, PP=10 would be
% preferable to PP=100 -- comparable or better conditioning with a much
% smaller intervention in the underlying physics.
%
% Uses the ORIGINAL (buggy rOuter=4e-6) starting state, unchanged, same
% as the PP=100 long run, so results are directly comparable step-for-
% step. Writes to its own output/log files only -- does not touch the
% PP=100 long run's files.

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
par.fluidPressurePenalty = 10;
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
save('out_finestep_UNSTEADY_PP10_LONG.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_PP10_LONG.mat\n');

fprintf('\n=== RESULT (fluidPressurePenalty=10, LONG horizon) ===\n');
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
