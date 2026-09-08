%% TEST_FINESTEP_UNSTEADY_WARMSTARTFIX3_FROM_RUNB_STATE1
% Identical setup to test_finestep_UNSTEADY_from_runB_state1.m (resume from
% Run B's own accepted state at t=6e-4, fine dt matching Run A, unsteady
% Stokes ON) -- run SEPARATELY, in parallel with that run, to verify the
% Aug 31 fix in solve_finite_def_solid.m: the top-of-Newton-iteration
% evaluation is now protected against element inversion (matching the
% line-search's own protection), and "Element inverted" now routes into
% the fsolve fallback instead of aborting the correction pass immediately.
% Writes to its OWN output/log files (out_finestep_UNSTEADY_WARMSTARTFIX3.mat)
% -- does not touch out_finestep_UNSTEADY.mat or anything the other run uses.

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
save('out_finestep_UNSTEADY_WARMSTARTFIX3.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_WARMSTARTFIX3.mat\n');

fprintf('\n=== RESULT (unsteady Stokes ON, both warm-start fixes + endothelium trust-region diagnostic ACTIVE) ===\n');
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
