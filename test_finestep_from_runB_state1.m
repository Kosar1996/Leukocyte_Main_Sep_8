%% TEST_FINESTEP_FROM_RUNB_STATE1
% Decisive check: is the pressure spike Run A saw near t=1.2e-3 a real
% transient that ANY sufficiently fine dt would catch (a temporal-
% resolution/sampling issue), or is it specific to Run A's own small-dt
% trajectory (a genuine dt-dependent physical difference)?
%
% Method: resume from Run B's OWN actual accepted state at step 1
% (t=6e-4, dt=6e-4 nominal) -- not Run A's state -- and take fine steps
% (dt=3e-4, matching Run A) from there. If THIS produces a spike near
% t=1.2e-3 similar to what Run A saw (~14882 Pa), that confirms the
% spike is a real feature reachable from Run B's own trajectory, just
% missed by its coarse dt -- a resolution issue. If it stays smooth
% (matching Run B's own actual step 2 result, ~65.72 Pa), that means
% Run A's spike is tied to Run A's specific small-dt path, a more
% surprising, genuinely different dt-dependence finding.

clc;

S = load('out_pure2dmac_dtlarge_7steps.mat');
out = S.out;

fprintf('Run B step 1: t=%.6e s, dt used=%.6e s\n', out.t(1), out.dtHist(1));
fprintf('Run B step 2 (actual, coarse dt): t=%.6e s, max|p|=%.4f Pa\n\n', ...
    out.t(2), max(abs(out.fluidHist{2}.p)));

state1 = out.stateHist{1};
par = out.par;
par.dt = 3e-4;              % fine dt, matching Run A
nFineSteps = 2;             % reach t = 6e-4 + 2*3e-4 = 1.2e-3, matching Run A's spike time
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
save('out_finestep_STEADY.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output (including z-profile pressure) to out_finestep_STEADY.mat\n');

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

fprintf('\n=== RESULT ===\n');
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
