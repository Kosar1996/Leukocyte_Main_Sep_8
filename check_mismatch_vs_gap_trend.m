%% CHECK_MISMATCH_VS_GAP_TREND
% Checks whether the interface %mismatch (normal + tangential, both
% interfaces) trends alongside the gap-narrowing behavior found in
% check_shape_degradation_over_time.m -- i.e. does the mismatch get worse
% over the same steps where the gap starts shrinking (steps ~17-40)?
% Pure post-processing of the already-saved 40-step run, no rerun.

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

nSteps = out.stopStep;
par = out.par;

parLmismatch = par;
if isfield(par, 'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par, 'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par, 'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

mismatchOpts = struct('nQuery', 41, 'epsFrac', 0.1, 'trimFrac', 0.05);

gapMin = nan(nSteps,1);
pctMaxEnL = nan(nSteps,1); pctMaxEtL = nan(nSteps,1);
pctMaxEnE = nan(nSteps,1); pctMaxEtE = nan(nSteps,1);

fprintf('Computing interface mismatch at each of %d steps (this takes a bit)...\n', nSteps);
for k = 1:nSteps
    d = out.diagHist{k};
    if ~isempty(d) && isfield(d,'gapMin'), gapMin(k) = d.gapMin; end

    st = out.stateHist{k};
    fl = out.fluidHist{k};
    dtStep = out.dtHist(k);

    try
        [leuko, endo] = compute_interface_traction_mismatch( ...
            out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);
        pctMaxEnL(k) = leuko.stats.pctMaxEn;
        pctMaxEtL(k) = leuko.stats.pctMaxEt;
        pctMaxEnE(k) = endo.stats.pctMaxEn;
        pctMaxEtE(k) = endo.stats.pctMaxEt;
    catch ME
        fprintf('  step %d: mismatch computation failed (%s)\n', k, ME.message);
    end

    if mod(k,10)==0, fprintf('  ...step %d done\n', k); end
end

fprintf('\n%6s | %10s | %8s %8s %8s %8s\n', 'step', 'gapMin[um]', 'EnL%', 'EtL%', 'EnE%', 'EtE%');
for k = 1:nSteps
    fprintf('%6d | %10.4f | %8.2f %8.2f %8.2f %8.2f\n', ...
        k, gapMin(k)*1e6, pctMaxEnL(k), pctMaxEtL(k), pctMaxEnE(k), pctMaxEtE(k));
end

fig = figure('Position',[100 100 900 600],'Color','w');
tl = tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');

ax1 = nexttile(tl);
plot(ax1, 1:nSteps, gapMin*1e6, 'k-o','LineWidth',1.5);
xlabel(ax1,'step'); ylabel(ax1,'min gap [\mum]'); grid(ax1,'on');
title(ax1,'Minimum gap (from earlier check)');

ax2 = nexttile(tl);
plot(ax2, 1:nSteps, pctMaxEnL, '-o','DisplayName','EnL'); hold(ax2,'on');
plot(ax2, 1:nSteps, pctMaxEtL, '-s','DisplayName','EtL');
plot(ax2, 1:nSteps, pctMaxEnE, '-^','DisplayName','EnE');
plot(ax2, 1:nSteps, pctMaxEtE, '-v','DisplayName','EtE');
xlabel(ax2,'step'); ylabel(ax2,'%% mismatch (max stat)'); grid(ax2,'on');
legend(ax2,'Location','best'); title(ax2,'Interface mismatch, all four quantities');

sgtitle(fig,'Does interface mismatch track the gap-narrowing trend?');
outPng = fullfile(fileparts(mfilename('fullpath')), 'mismatch_vs_gap_trend.png');
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('\nSaved: %s\nSMOKE_TEST_STATUS: OK\n', outPng);
