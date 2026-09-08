%% TEST_FINESTEP_UNSTEADY_FSOLVEBUDGET_PPFIX_FROM_RUNB_STATE1
% Combines the two Issue 1 / Issue 3 fixes to test a physics-grounded
% hypothesis for why the fsolve-budget fix alone failed (Sep 2): the
% escalating traction near the crisis point (step 3, pass 14, ~11,796 Pa,
% then further escalating to 28,090 Pa) is consistent with the classic
% thin-film lubrication pressure singularity as gap h -> 0 (pressure
% scales roughly like 1/h^3 for a fixed flow constraint). That would mean
% the endothelium Newton solve's exitflag=-2 failure there (a GENUINE
% non-convergence, confirmed distinct from budget exhaustion) reflects
% the underlying linearized problem itself becoming ill-conditioned in
% that regime, not the solver merely needing more iterations. Since
% par.fluidPressurePenalty=100 has already been shown (Sep 2) to improve
% the fluid matrix's conditioning in exactly this near-floor regime, it
% may also be what actually lets the raised fsolve budget do its job --
% not by itself, but by first regularizing the underlying near-singular
% response the Newton solve is fighting.
%
% Uses the ORIGINAL (buggy rOuter=4e-6) starting state, unchanged, same
% as the original fsolve-budget test, so this isolates "budget alone" vs
% "budget + PP=100" against the exact same known crisis point (step 3,
% pass 14, ~11,796 Pa). Writes to its own output/log files only -- does
% not touch the original FSOLVEBUDGET test's files.

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

fprintf('BEFORE FIX: fsolve MaxIterations=200 (hardcoded), MaxFunctionEvaluations=2000 (hardcoded)\n');
par.solidFsolveMaxIterations = 2000;
par.solidFsolveMaxFunctionEvaluations = 20000;
fprintf('AFTER FIX:  par.solidFsolveMaxIterations=%d, par.solidFsolveMaxFunctionEvaluations=%d\n', ...
    par.solidFsolveMaxIterations, par.solidFsolveMaxFunctionEvaluations);

par.fluidPressurePenalty = 100;
fprintf('ADDED:      par.fluidPressurePenalty = %g\n', par.fluidPressurePenalty);

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
save('out_finestep_UNSTEADY_FSOLVEBUDGET_PPFIX.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_FSOLVEBUDGET_PPFIX.mat\n');

fprintf('\n=== RESULT (raised fsolve budget + fluidPressurePenalty=100 ACTIVE) ===\n');
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
