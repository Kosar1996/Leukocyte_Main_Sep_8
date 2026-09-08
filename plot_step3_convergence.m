%% PLOT_STEP3_CONVERGENCE
% Runs the traction-correction feedback loop (relax=0.30, maxIterations=30)
% and plots the %mismatch trajectory for all four interface quantities
% (leukocyte normal/tangential, endothelium normal/tangential) vs.
% correction-pass index, against the mismatchThresholdPct target line.
%

% reaches a genuine, stable fixed point well above the 10% target, rather
% than slowly converging toward it.
%
% Requires softlube_run_case_global_coupled.m's out.tractionCorrectionHistory
% (added so the convergeInfo the feedback function already computes isn't
% discarded).

clc; close all;
cd(fileparts(mfilename('fullpath')));

cfg = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', true);

testRelax = 0.3;
testMaxIter = 30;
cfg.fluid.bodyFittedTractionCorrectionRelax = testRelax;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = testRelax;
cfg.fluid.maxBodyFittedTractionCorrections = testMaxIter;
cfg.parOverrides.maxBodyFittedTractionCorrections = testMaxIter;

fprintf('Running relax=%.2f, maxIterations=%d ...\n', testRelax, testMaxIter);
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
target = 10; % opts.mismatchThresholdPct default

fig = figure('Position', [100 100 900 700], 'Color', 'w');
tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, sprintf('Traction-correction convergence (relax=%.2f, %d passes)', testRelax, testMaxIter));

panelData = {pctEnL, 'Leukocyte, normal (EnL)'; ...
             pctEtL, 'Leukocyte, tangential (EtL)'; ...
             pctEnE, 'Endothelium, normal (EnE)'; ...
             pctEtE, 'Endothelium, tangential (EtE)'};

for i = 1:4
    ax = nexttile(tl);
    plot(ax, iters, panelData{i,1}, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.10 0.30 0.75]);
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

outPng = fullfile(fileparts(mfilename('fullpath')), 'step3_convergence.png');
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('Saved: %s\n', outPng);
fprintf('converged=%d, finalWorstPct=%.2f%%\n', convergeInfo.converged, convergeInfo.finalWorstPct);
