%% PLOT_STEP3_WIDEDOMAIN_LONGRUN
% Same setup as plot_step3_3pct_10pass_2D_widedomain_lowrelax.m (wide
% domain, fixed relax=0.05, 3% target -- the only one of the three tested
% strategies that didn't blow up), but with a much bigger pass cap to see
% whether it eventually converges or is genuinely stuck. At pass 30 it had
% bottomed out at 23.17% (pass 14) and was drifting back up to 26.48% --
% this run checks whether that's a slow oscillation that comes back down,
% or a one-way drift.

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
iters  = [h.iter];
pctEnL = [h.pctEnL];
pctEtL = [h.pctEtL];
pctEnE = [h.pctEnE];
pctEtE = [h.pctEtE];
worstAll = [h.worstPct];

fprintf('\nFull history (WIDE domain, relax=%.2f, %d passes):\n', testRelax, testMaxIter);
fprintf('%6s %8s %8s %8s %8s %8s\n', 'iter', 'EnL', 'EtL', 'EnE', 'EtE', 'worst');
for i = 1:numel(iters)
    fprintf('%6d %8.2f %8.2f %8.2f %8.2f %8.2f\n', iters(i), pctEnL(i), pctEtL(i), pctEnE(i), pctEtE(i), worstAll(i));
end
[bestWorst, bestIdx] = min(worstAll);
fprintf('\nBest worst-case mismatch: %.2f%% at pass %d\n', bestWorst, iters(bestIdx));
fprintf('converged=%d, finalWorstPct=%.2f%%\n', convergeInfo.converged, convergeInfo.finalWorstPct);

fig = figure('Position',[100 100 900 600],'Color','w');
plot(iters, pctEnL, 'o-', 'DisplayName','EnL'); hold on;
plot(iters, pctEtL, 's-', 'DisplayName','EtL');
plot(iters, pctEnE, '^-', 'DisplayName','EnE');
plot(iters, pctEtE, 'd-', 'DisplayName','EtE');
yline(3, '--r', '3% target');
grid on; legend('Location','best'); xlabel('correction pass'); ylabel('%% mismatch');
title(sprintf('Wide domain, relax=%.2f, %d passes', testRelax, testMaxIter));
outPng = fullfile(fileparts(mfilename('fullpath')), 'step3_widedomain_longrun.png');
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('Saved: %s\n', outPng);
fprintf('SMOKE_TEST_STATUS: OK\n');
