%% PLOT_STEP3_3PCT_10PASS_2D_WIDEDOMAIN
% Tests whether the EnE/EtE plateau (Item 3's symptom) is actually caused
% by the endothelium domain truncation (Item 6), by rerunning the exact
% same 3%/10-pass convergence test as plot_step3_3pct_10pass_2D.m, but on
% the wide-domain prestress (solid_endo_P300_wide.mat, z=[-2,6]um)
% instead of the original narrow one. If EnE converges much closer to 3%
% here, that confirms the two items share a root cause and Item 6's fix
% is also most of Item 3's fix.

clc; close all;
cd(fileparts(mfilename('fullpath')));

cfg = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', false);
cfg.geometry.endotheliumPrestressFile = fullfile(fileparts(mfilename('fullpath')), 'solid_endo_P300_wide.mat');

testRelax = 0.3;
testMaxIter = 10;
testThresholdPct = 3;
cfg.fluid.bodyFittedTractionCorrectionRelax = testRelax;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = testRelax;
cfg.fluid.maxBodyFittedTractionCorrections = testMaxIter;
cfg.parOverrides.maxBodyFittedTractionCorrections = testMaxIter;
cfg.fluid.tractionCorrectionMismatchThresholdPct = testThresholdPct;
cfg.parOverrides.tractionCorrectionMismatchThresholdPct = testThresholdPct;

fprintf('Running 2D, WIDE domain, threshold=%.1f%%, maxIterations=%d ...\n', testThresholdPct, testMaxIter);
out = softlube_run_case_global_coupled(cfg);

convergeInfo = out.tractionCorrectionHistory{1};
h = convergeInfo.history;
iters  = [h.iter];
pctEnL = [h.pctEnL];
pctEtL = [h.pctEtL];
pctEnE = [h.pctEnE];
pctEtE = [h.pctEtE];

fprintf('\nFull history (WIDE domain):\n');
fprintf('%6s %8s %8s %8s %8s %8s\n', 'iter', 'EnL', 'EtL', 'EnE', 'EtE', 'worst');
for i = 1:numel(iters)
    fprintf('%6d %8.2f %8.2f %8.2f %8.2f %8.2f\n', iters(i), pctEnL(i), pctEtL(i), pctEnE(i), pctEtE(i), h(i).worstPct);
end
fprintf('\nconverged=%d, finalWorstPct=%.2f%%\n', convergeInfo.converged, convergeInfo.finalWorstPct);
fprintf('SMOKE_TEST_STATUS: OK\n');
