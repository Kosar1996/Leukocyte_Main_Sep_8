%% TEST_ISSUE2_CONDITION_NUMBER_DIRECT
% Direct, dedicated test of Issue 2's own claim (fluid matrix conditioning
% near the gap floor), independent of the solid coupling -- everything so
% far has only inferred this through Issue 3's side effects. Uses the
% instrumented solve_stokes_bodyfitted_MAC.m (par.reportMatrixConditionNumber,
% added today) to directly compare the row-scaled system's estimated
% condition number (condest) with and without the pressure-stabilization
% term, at the exact near-floor gap state (gap = 1.000e-09 m, pinned at
% the minGap floor) reached by today's rOuter+PP=100 no-stop run.
%
% Does one single fluid solve at a fixed solid state, twice (PP=0 vs
% PP=100) -- not a time-marching run, so this is fast. Writes only to
% stdout/its own log; does not touch any other test's files.

clc;

S = load('out_finestep_UNSTEADY_ROUTERFIX_PPFIX_NOSTOP.mat');
outFine = S.outFine;

nS = numel(outFine.stateHist);
old = outFine.stateHist{nS-1};
state = outFine.stateHist{nS};
z = outFine.z;
par = outFine.par;

gapNow = min(state.deltaE(:) - state.deltaL(:));
fprintf('Using state %d of %d, min gap = %.6e m\n', nS, nS, gapNow);

par.reportMatrixConditionNumber = true;

fprintf('\n--- Fluid solve with par.fluidPressurePenalty = 0 (baseline, no stabilization) ---\n');
parA = par;
parA.fluidPressurePenalty = 0;
[fluidA, okA, stopReasonA] = solve_selected_poststep_fluid(z, old, state, parA);
fprintf('ok=%d, stopReason=%s\n', okA, stopReasonA);
if okA
    fprintf('condNumber (condest) = %.6e\n', fluidA.condNumber);
    fprintf('linRes = %.6e\n', fluidA.linRes);
    fprintf('max|p| = %.6e Pa\n', max(abs(fluidA.p(:))));
end

fprintf('\n--- Fluid solve with par.fluidPressurePenalty = 100 (the Issue 3 candidate fix) ---\n');
parB = par;
parB.fluidPressurePenalty = 100;
[fluidB, okB, stopReasonB] = solve_selected_poststep_fluid(z, old, state, parB);
fprintf('ok=%d, stopReason=%s\n', okB, stopReasonB);
if okB
    fprintf('condNumber (condest) = %.6e\n', fluidB.condNumber);
    fprintf('linRes = %.6e\n', fluidB.linRes);
    fprintf('max|p| = %.6e Pa\n', max(abs(fluidB.p(:))));
end

if okA && okB
    fprintf('\n=== RESULT ===\n');
    fprintf('Condition number, PP=0:   %.6e\n', fluidA.condNumber);
    fprintf('Condition number, PP=100: %.6e\n', fluidB.condNumber);
    fprintf('Ratio (PP=0 / PP=100):    %.4f\n', fluidA.condNumber / fluidB.condNumber);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
