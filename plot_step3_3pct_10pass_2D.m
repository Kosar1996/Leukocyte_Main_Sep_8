%% PLOT_STEP3_3PCT_10PASS_2D

% 3% mismatch threshold (instead of the 10% default) and a 10-pass cap,
% on the FULL 2D solver.

clc; close all;
cd(fileparts(mfilename('fullpath')));

cfg = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', false);

testRelax = 0.3;
testMaxIter = 10;
testThresholdPct = 3;
cfg.fluid.bodyFittedTractionCorrectionRelax = testRelax;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = testRelax;
cfg.fluid.maxBodyFittedTractionCorrections = testMaxIter;
cfg.parOverrides.maxBodyFittedTractionCorrections = testMaxIter;
cfg.fluid.tractionCorrectionMismatchThresholdPct = testThresholdPct;
cfg.parOverrides.tractionCorrectionMismatchThresholdPct = testThresholdPct;

fprintf('Running 2D, threshold=%.1f%%, maxIterations=%d ...\n', testThresholdPct, testMaxIter);
out = softlube_run_case_global_coupled(cfg);

convergeInfo = out.tractionCorrectionHistory{1};
if isempty(convergeInfo)
    error('No traction-correction history recorded for step 1 -- check useFeedbackTractionCorrection is on.');
end

h = convergeInfo.history;
iters  = [h.iter];
pctEnL = [h.pctEnL];
pctEtL = [h.pctEtL];
pctEnE = [h.pctEnE];
pctEtE = [h.pctEtE];
target = testThresholdPct;

fig = figure('Position', [100 100 900 700], 'Color', 'w');
tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, sprintf('Traction-correction convergence -- 2D (relax=%.2f, %d passes, target=%.1f%%)', ...
    testRelax, testMaxIter, testThresholdPct));

panelData = {pctEnL, 'Leukocyte, normal (EnL)'; ...
             pctEtL, 'Leukocyte, tangential (EtL)'; ...
             pctEnE, 'Endothelium, normal (EnE)'; ...
             pctEtE, 'Endothelium, tangential (EtE)'};

for i = 1:4
    ax = nexttile(tl);
    plot(ax, iters, panelData{i,1}, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.75 0.35 0.10]);
    hold(ax, 'on');
    yline(ax, target, '--', sprintf('%g%% target', target), 'Color', [0.75 0.10 0.10], ...
        'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
    hold(ax, 'off');
    grid(ax, 'on');
    xlabel(ax, 'correction pass');
    ylabel(ax, '%% mismatch');
    title(ax, panelData{i,2});
    xlim(ax, [0 max(iters)]);
end

outPng = fullfile(fileparts(mfilename('fullpath')), 'step3_3pct_10pass_2D.png');
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('Saved: %s\n', outPng);
fprintf('converged=%d, finalWorstPct=%.2f%%\n', convergeInfo.converged, convergeInfo.finalWorstPct);
