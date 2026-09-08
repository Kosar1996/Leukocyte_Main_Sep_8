%% TEST_FINESTEP_UNSTEADY_ROUTERFIX_FROM_RUNB_STATE1
% Same setup as test_finestep_UNSTEADY_WARMSTARTFIX2_from_runB_state1.m
% (both solid-solver fixes active), PLUS testing the rOuter/global1DOuterRadius
% fix found Sep 1: the loaded state's par.global1DOuterRadius is 4e-6 (a
% confirmed config bug -- should be 15e-6, matching every comment in the
% source config files and the function's own built-in default in
% global_1d_outer_radius.m). This clamps the endothelium's assumed radius
% to never exceed 4 microns anywhere in global_1d_interface_kinematics.m,
% including inside the active gap region during a genuine pressure crisis,
% when the physically correct response is for the endothelium to bulge
% outward past that point. Testing directly: does correcting this to
% 15e-6 change the crisis behavior at the same step-5 spike point?
% Writes to its OWN output file -- does not touch any other run's files.

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
nFineSteps = 5;
par.tEnd = state1.t + nFineSteps * par.dt;

fprintf('BEFORE FIX: par.global1DOuterRadius = %.4e m\n', par.global1DOuterRadius);
par.global1DOuterRadius = 15e-6;
fprintf('AFTER FIX:  par.global1DOuterRadius = %.4e m\n', par.global1DOuterRadius);

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
save('out_finestep_UNSTEADY_ROUTERFIX.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_ROUTERFIX.mat\n');

fprintf('\n=== RESULT (unsteady Stokes ON, both warm-start fixes + rOuter fix ACTIVE) ===\n');
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
