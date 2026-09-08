clc; close all;
cd(fileparts(mfilename('fullpath')));

cfg = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', false);
cfg.geometry.endotheliumPrestressFile = fullfile(fileparts(mfilename('fullpath')), 'solid_endo_P300_wide.mat');

testRelax = 0.05;
testMaxIter = 30;
testThresholdPct = 3;
cfg.fluid.bodyFittedTractionCorrectionRelax = testRelax;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = testRelax;
cfg.fluid.maxBodyFittedTractionCorrections = testMaxIter;
cfg.parOverrides.maxBodyFittedTractionCorrections = testMaxIter;
cfg.fluid.tractionCorrectionMismatchThresholdPct = testThresholdPct;
cfg.parOverrides.tractionCorrectionMismatchThresholdPct = testThresholdPct;
cfg.fluid.useAitkenTractionCorrectionRelax = true;
cfg.parOverrides.useAitkenTractionCorrectionRelax = true;

fprintf('Running 2D, WIDE domain, Aitken relax (init=%.2f), threshold=%.1f%%, maxIterations=%d ...\n', testRelax, testThresholdPct, testMaxIter);
out = softlube_run_case_global_coupled(cfg);

convergeInfo = out.tractionCorrectionHistory{1};
h = convergeInfo.history;
iters  = [h.iter];
pctEnL = [h.pctEnL];
pctEtL = [h.pctEtL];
pctEnE = [h.pctEnE];
pctEtE = [h.pctEtE];

fprintf('\nFull history (WIDE domain, Aitken):\n');
fprintf('%6s %8s %8s %8s %8s %8s\n', 'iter', 'EnL', 'EtL', 'EnE', 'EtE', 'worst');
for i = 1:numel(iters)
    fprintf('%6d %8.2f %8.2f %8.2f %8.2f %8.2f\n', iters(i), pctEnL(i), pctEtL(i), pctEnE(i), pctEtE(i), h(i).worstPct);
end
fprintf('\nconverged=%d, finalWorstPct=%.2f%%\n', convergeInfo.converged, convergeInfo.finalWorstPct);
fprintf('SMOKE_TEST_STATUS: OK\n');
