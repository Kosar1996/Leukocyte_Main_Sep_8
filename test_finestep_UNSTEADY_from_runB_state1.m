%% TEST_FINESTEP_UNSTEADY_FROM_RUNB_STATE1
% Same decisive setup as test_finestep_from_runB_state1.m (resume from
% Run B's own accepted state at t=6e-4, take fine steps matching Run A's
% dt), but with the new unsteady-Stokes term switched ON. If the fix is
% doing what it's meant to, the pressure spike seen in the steady version
% (17,306 Pa) should be smoothed/damped here, not eliminated by changing
% the physics being tested -- same starting state, same fine dt, only
% difference is par.useUnsteadyStokes.

clc;

S = load('out_pure2dmac_dtlarge_7steps.mat');
out = S.out;

fprintf('Run B step 1: t=%.6e s, dt used=%.6e s\n', out.t(1), out.dtHist(1));
fprintf('Run B step 2 (actual, coarse dt, steady): t=%.6e s, max|p|=%.4f Pa\n', ...
    out.t(2), max(abs(out.fluidHist{2}.p)));
fprintf('Fine-step steady-Stokes result (already measured): 17306.22 Pa\n\n');

state1 = out.stateHist{1};
par = out.par;
par.dt = 3e-4;
par.useUnsteadyStokes = true;
par.rho = 1000;
nFineSteps = 2;
par.tEnd = state1.t + nFineSteps * par.dt;

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
save('out_finestep_UNSTEADY.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output (including z-profile pressure) to out_finestep_UNSTEADY.mat\n');

fprintf('\n=== DIAGNOSTICS ===\n');
fprintf('outFine.stopStep = %d\n', outFine.stopStep);
fprintf('numel(outFine.fluidHist) = %d\n', numel(outFine.fluidHist));
fprintf('numel(outFine.t) = %d\n', numel(outFine.t));
for k = 1:numel(outFine.fluidHist)
    fprintf('  fluidHist{%d}: class=%s', k, class(outFine.fluidHist{k}));
    if isstruct(outFine.fluidHist{k})
        fprintf(', has p field=%d', isfield(outFine.fluidHist{k}, 'p'));
    end
    fprintf('\n');
end

fprintf('\n=== RESULT (unsteady Stokes ON) ===\n');
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
