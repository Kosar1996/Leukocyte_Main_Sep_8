%% TEST_FINESTEP_UNSTEADY_FSOLVEBUDGET_FROM_RUNB_STATE1

% fail outright under large traction (~11,800 Pa) even with the trust
% region already fixed (that hypothesis was tested directly Sep 1 and
% disproven -- scaling the ceiling up 3.9x did not prevent the failure).
%
% New, different hypothesis: the failure showed fsolve exitflag=0, which
% per MATLAB's own documentation means the iteration/function-evaluation
% budget was exhausted (MaxIterations=200/MaxFunctionEvaluations=2000,
% both previously hardcoded), not that a genuine dead end was detected
% (that would be exitflag=-2/-3). Testing directly: does simply raising
% the budget (2000 iterations / 20000 function evaluations, both now
% overridable via par.solidFsolveMaxIterations/
% solidFsolveMaxFunctionEvaluations) let the endothelium's fsolve
% fallback actually finish converging at the same failure point?
%
% Uses the ORIGINAL (buggy rOuter=4e-6) starting state, unchanged, so
% this isolates the fsolve-budget variable alone against the exact same
% known failure point (step 3, pass 14, traction ~11,796 Pa) used for
% the trust-region test. Does not touch the rOuter-fix or relax-fix
% tests or their files. Writes to its own output/log files only.

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
save('out_finestep_UNSTEADY_FSOLVEBUDGET.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_FSOLVEBUDGET.mat\n');

fprintf('\n=== RESULT (unsteady Stokes ON, both warm-start fixes + raised fsolve budget ACTIVE) ===\n');
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
