%% TEST_ISSUE2_CONDITION_NUMBER_SWEEP
% Follow-up to test_issue2_condition_number_direct.m: that test showed
% PP=100 improves the condition number 8.76x at the near-floor gap state.
% This sweeps a range of magnitudes (0, 1, 10, 100, 1000, 10000, 100000)
% at the SAME fixed state to characterize the trend -- is 100 near a
% point of diminishing returns, or would a smaller value already capture
% most of the benefit (minimizing any risk of physics distortion, per
% the PP=1e5 pressure-suppression concern found Sep 2)?
%
% Same single-solve-per-value approach as the direct test -- fast, no
% time-marching. Writes only to stdout.

clc;

S = load('out_finestep_UNSTEADY_ROUTERFIX_PPFIX_NOSTOP.mat');
outFine = S.outFine;

nS = numel(outFine.stateHist);
old = outFine.stateHist{nS-1};
state = outFine.stateHist{nS};
z = outFine.z;
parBase = outFine.par;
parBase.reportMatrixConditionNumber = true;

gapNow = min(state.deltaE(:) - state.deltaL(:));
fprintf('Using state %d of %d, min gap = %.6e m\n\n', nS, nS, gapNow);

ppValues = [0, 1, 10, 100, 1000, 10000, 100000];
condVals = nan(size(ppValues));
pVals = nan(size(ppValues));

for k = 1:numel(ppValues)
    parK = parBase;
    parK.fluidPressurePenalty = ppValues(k);
    [fluidK, okK] = solve_selected_poststep_fluid(z, old, state, parK);
    if okK
        condVals(k) = fluidK.condNumber;
        pVals(k) = max(abs(fluidK.p(:)));
        fprintf('PP=%-10g condNumber=%.6e  max|p|=%.6e Pa\n', ppValues(k), condVals(k), pVals(k));
    else
        fprintf('PP=%-10g FAILED\n', ppValues(k));
    end
end

fprintf('\n=== SUMMARY ===\n');
for k = 1:numel(ppValues)
    if ~isnan(condVals(k))
        fprintf('PP=%-10g  cond=%.4e  ratio-vs-PP0=%.3fx  max|p|=%.4e Pa\n', ...
            ppValues(k), condVals(k), condVals(1)/condVals(k), pVals(k));
    end
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
