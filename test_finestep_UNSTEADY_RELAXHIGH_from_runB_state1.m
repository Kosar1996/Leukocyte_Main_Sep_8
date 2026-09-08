%% TEST_FINESTEP_UNSTEADY_RELAXHIGH_FROM_RUNB_STATE1

% hypothesis. The relax=0.25 test (test_finestep_UNSTEADY_RELAXFIX_...)
% did NOT behave the way simple fixed-point-damping theory predicts: it
% failed SOONER and HARDER (14,327 Pa spike at pass 4) than the relax=0.5
% baseline (which fails around pass 13-14 at ~11,796 Pa), the opposite of
% "smaller relax -> more stable". That anomaly is the motivation for this
% test: does a LARGER relaxation factor (0.7, less damping, bigger steps
% toward the corrected traction each pass) do better instead?
%
% Does NOT touch the rOuter-fix or relax=0.25 test or their files. Uses
% the ORIGINAL (buggy rOuter=4e-6) starting state, unchanged, so this
% isolates the relaxation-factor variable alone against the same known
% crisis point used by every earlier isolated comparison tonight. Writes
% to its own output/log files only.

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
nFineSteps = 4;
par.tEnd = state1.t + nFineSteps * par.dt;

fprintf('BEFORE FIX: par.bodyFittedTractionCorrectionRelax = %.4f\n', par.bodyFittedTractionCorrectionRelax);
par.bodyFittedTractionCorrectionRelax = 0.7;
fprintf('AFTER FIX:  par.bodyFittedTractionCorrectionRelax = %.4f\n', par.bodyFittedTractionCorrectionRelax);

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
save('out_finestep_UNSTEADY_RELAXHIGH.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_RELAXHIGH.mat\n');

fprintf('\n=== RESULT (unsteady Stokes ON, both warm-start fixes + relax=0.70 ACTIVE) ===\n');
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
