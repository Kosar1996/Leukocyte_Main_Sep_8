clc; close all;
cd(fileparts(mfilename('fullpath')));

cfg = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', false);
cfg.geometry.endotheliumPrestressFile = fullfile(fileparts(mfilename('fullpath')), 'solid_endo_P300_wide.mat');

testRelax = 0.05;
testMaxIter = 40;
testThresholdPct = 3;
cfg.fluid.bodyFittedTractionCorrectionRelax = testRelax;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = testRelax;
cfg.fluid.maxBodyFittedTractionCorrections = testMaxIter;
cfg.parOverrides.maxBodyFittedTractionCorrections = testMaxIter;
cfg.fluid.tractionCorrectionMismatchThresholdPct = testThresholdPct;
cfg.parOverrides.tractionCorrectionMismatchThresholdPct = testThresholdPct;

fprintf('Running 2D, WIDE domain, NEW local-fluid-relative %%mismatch definition, relax=%.2f, threshold=%.1f%%, maxIterations=%d ...\n', testRelax, testThresholdPct, testMaxIter);
tic;
out = softlube_run_case_global_coupled(cfg);
fprintf('Elapsed: %.1f s\n', toc);

convergeInfo = out.tractionCorrectionHistory{1};
h = convergeInfo.history;
n = numel(h);

fprintf('\nAll %d passes:\n', n);
fprintf('%6s %8s %8s %8s %8s %12s %12s %12s %12s\n', 'iter', 'EnL', 'EtL', 'EnE', 'EtE', 'rawEtL[Pa]', 'ttFluidL[Pa]', 'rawEtE[Pa]', 'ttFluidE[Pa]');
for i = 1:n
    fprintf('%6d %8.2f %8.2f %8.2f %8.2f %12.4f %12.4f %12.4f %12.4f\n', i, h(i).pctEnL, h(i).pctEtL, h(i).pctEnE, h(i).pctEtE, ...
        h(i).rawEtL_Pa, h(i).ttFluidL_atMax, h(i).rawEtE_Pa, h(i).ttFluidE_atMax);
end

fprintf('\nFinal pass (%d): EnL=%.2f EtL=%.2f EnE=%.2f EtE=%.2f\n', n, h(end).pctEnL, h(end).pctEtL, h(end).pctEnE, h(end).pctEtE);

fig = figure('Position',[100 100 900 500]);
iters = [h.iter];
plot(iters, [h.pctEnL], '-o', iters, [h.pctEtL], '-s', iters, [h.pctEnE], '-^', iters, [h.pctEtE], '-d');
xlabel('correction pass'); ylabel('%% mismatch');
title(sprintf('Wide domain, relax=%.2f, %d passes, NEW local-fluid-relative definition', testRelax, testMaxIter));
legend('EnL','EtL','EnE','EtE','Location','best');
grid on;
yline(testThresholdPct, 'r--', sprintf('%g%% target', testThresholdPct));
saveas(fig, 'widedomain_newdef_40pass_plot.png');

save('widedomain_newdef_40pass_history.mat', 'h', 'testRelax', 'testMaxIter');
fprintf('\nSMOKE_TEST_STATUS: OK\n');
