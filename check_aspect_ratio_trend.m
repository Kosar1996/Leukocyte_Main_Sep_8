%% CHECK_ASPECT_RATIO_TREND
% Tracks the fluid gap-cell aspect ratio (dz / radial cell size) across
% all 40 steps, to see when it crosses into a range likely to cause real
% numerical noise in cross-derivative (sigma_rz) stress recovery, and
% whether that timing lines up with when the noisy sigma_rz behavior
% was actually observed (steps 17+, contact zone locked in).

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

zGrid = out.z(:);
nSteps = out.stopStep;
gapMin = nan(nSteps,1);
aspectRatio = nan(nSteps,1);
dzVals = nan(nSteps,1);

for k = 1:nSteps
    fl = out.fluidHist{k};
    meshF = fl.meshF;
    deltaE_k = out.deltaEHist(:,k);
    deltaL_k = out.deltaLHist(:,k);
    gap = deltaE_k - deltaL_k;
    [minGap, ~] = min(gap);
    gapMin(k) = minGap;
    radialCell = minGap / (meshF.Nr - 1);
    aspectRatio(k) = meshF.dz / radialCell;
    dzVals(k) = meshF.dz;
end

fprintf('%6s %12s %12s %14s\n', 'step', 'gapMin[nm]', 'dz[nm]', 'aspectRatio');
for k = 1:nSteps
    fprintf('%6d %12.4f %12.2f %14.1f\n', k, gapMin(k)*1e9, dzVals(k)*1e9, aspectRatio(k));
end

fig = figure('Position',[100 100 900 600],'Color','w');
yyaxis left
semilogy(1:nSteps, aspectRatio, 'o-');
ylabel('aspect ratio (dz / radial cell)');
yline(10, 'r--', '10:1 (common quality threshold)');
yyaxis right
plot(1:nSteps, gapMin*1e9, 's-');
ylabel('min gap [nm]');
xlabel('step');
grid on;
title('Fluid gap-cell aspect ratio vs. gap closure');
outPng = fullfile(fileparts(mfilename('fullpath')), 'aspect_ratio_trend.png');
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('\nSaved: %s\n', outPng);

idxCross = find(aspectRatio >= 10, 1);
if ~isempty(idxCross)
    fprintf('\nAspect ratio first crosses 10:1 at step %d (gap=%.2f nm)\n', idxCross, gapMin(idxCross)*1e9);
end
fprintf('SMOKE_TEST_STATUS: OK\n');
