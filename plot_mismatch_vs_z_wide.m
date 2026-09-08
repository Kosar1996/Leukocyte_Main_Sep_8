%% PLOT_MISMATCH_VS_Z_WIDE
% Same exercise as plot_mismatch_vs_z_narrow.m but for the WIDE domain
% (z=[-2,6]), after the 150-pass correction, so the plot spans both the
% thin-gap edges AND the wide-gap bulk across the whole domain.

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

k = out.stopStep;
st = out.stateHist{k};
fl = out.fluidHist{k};
par = out.par;
dtStep = out.dtHist(k);

parLmismatch = par;
if isfield(par, 'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par, 'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par, 'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

mismatchOpts = struct('nQuery', 91, 'epsFrac', 0.1, 'trimFrac', 0.05);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);

pctEnL = 100*abs(leuko.en)/leuko.stats.refScale;
pctEtL = 100*abs(leuko.et)/leuko.stats.refScale;
pctEnE = 100*abs(endo.en)/endo.stats.refScale;
pctEtE = 100*abs(endo.et)/endo.stats.refScale;

fig = figure('Position',[100 100 1000 700]);

subplot(2,1,1);
plot(leuko.z*1e6, pctEnL, 'b-o', 'MarkerSize',3, 'DisplayName','Normal (En)'); hold on;
plot(leuko.z*1e6, pctEtL, 'r-s', 'MarkerSize',3, 'DisplayName','Tangential (Et)');
xlabel('z [\mum]'); ylabel('% mismatch');
title('Leukocyte-fluid interface: %mismatch vs z (WIDE domain, z=[-2,6], after 150-pass correction)');
legend('Location','best'); grid on;

subplot(2,1,2);
plot(endo.z*1e6, pctEnE, 'b-o', 'MarkerSize',3, 'DisplayName','Normal (En)'); hold on;
plot(endo.z*1e6, pctEtE, 'r-s', 'MarkerSize',3, 'DisplayName','Tangential (Et)');
xlabel('z [\mum]'); ylabel('% mismatch');
title('Endothelium-fluid interface: %mismatch vs z (WIDE domain, z=[-2,6], after 150-pass correction)');
legend('Location','best'); grid on;

saveas(fig, 'mismatch_vs_z_wide_domain.png');
fprintf('Saved mismatch_vs_z_wide_domain.png\n');
fprintf('Max En leuko=%.2f%% at z=%.4f, Max En endo=%.2f%% at z=%.4f\n', ...
    max(pctEnL), leuko.z(find(pctEnL==max(pctEnL),1))*1e6, max(pctEnE), endo.z(find(pctEnE==max(pctEnE),1))*1e6);
fprintf('SMOKE_TEST_STATUS: OK\n');
