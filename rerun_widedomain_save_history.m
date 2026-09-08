clc; close all;
cd(fileparts(mfilename('fullpath')));

cfg = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', false);
cfg.geometry.endotheliumPrestressFile = fullfile(fileparts(mfilename('fullpath')), 'solid_endo_P300_wide.mat');

testRelax = 0.05;
testMaxIter = 150;
testThresholdPct = 3;
cfg.fluid.bodyFittedTractionCorrectionRelax = testRelax;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = testRelax;
cfg.fluid.maxBodyFittedTractionCorrections = testMaxIter;
cfg.parOverrides.maxBodyFittedTractionCorrections = testMaxIter;
cfg.fluid.tractionCorrectionMismatchThresholdPct = testThresholdPct;
cfg.parOverrides.tractionCorrectionMismatchThresholdPct = testThresholdPct;

fprintf('Running 2D, WIDE domain, relax=%.2f, threshold=%.1f%%, maxIterations=%d ...\n', testRelax, testThresholdPct, testMaxIter);
tic;
out = softlube_run_case_global_coupled(cfg);
fprintf('Elapsed: %.1f s\n', toc);

convergeInfo = out.tractionCorrectionHistory{1};
h = convergeInfo.history;

fprintf('\nLast 10 passes:\n');
fprintf('%6s %8s %8s %8s %8s\n', 'iter', 'EnL', 'EtL', 'EnE', 'EtE');
n = numel(h);
for i = max(1,n-9):n
    fprintf('%6d %8.2f %8.2f %8.2f %8.2f\n', i, h(i).pctEnL, h(i).pctEtL, h(i).pctEnE, h(i).pctEtE);
end

fprintf('\nMax over ALL 150 passes: EnL=%.2f EtL=%.2f EnE=%.2f EtE=%.2f\n', ...
    max([h.pctEnL]), max([h.pctEtL]), max([h.pctEnE]), max([h.pctEtE]));
fprintf('Max over passes 20-150 (post-transient): EnL=%.2f EtL=%.2f EnE=%.2f EtE=%.2f\n', ...
    max([h(20:end).pctEnL]), max([h(20:end).pctEtL]), max([h(20:end).pctEnE]), max([h(20:end).pctEtE]));
fprintf('Final pass (150): EnL=%.2f EtL=%.2f EnE=%.2f EtE=%.2f\n', ...
    h(end).pctEnL, h(end).pctEtL, h(end).pctEnE, h(end).pctEtE);

save('widedomain_150pass_history.mat', 'h', 'testRelax', 'testMaxIter');
fprintf('\nSMOKE_TEST_STATUS: OK\n');
