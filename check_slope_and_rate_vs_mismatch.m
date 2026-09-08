%% CHECK_SLOPE_AND_RATE_VS_MISMATCH
% Checks the two remaining candidate explanations for the mismatch blowup
% at small gap (steps 37-40):
%   1. Rate of geometric change -- is the gap closing faster (relatively)
%      near the end, giving the correction less effective time to track?
%   2. Local interface slope at the contact zone -- does it steepen as
%      the gap narrows, making the normal/tangential rotation (which
%      depends on slope) more sensitive?

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

nSteps = out.stopStep;
zGrid = out.z(:);
gapMin = nan(nSteps,1);
slopeAtContact = nan(nSteps,1);
relRate = nan(nSteps,1);

dz_fd = 1e-8;

for k = 1:nSteps
    d = out.diagHist{k};
    if ~isempty(d) && isfield(d,'gapMin'), gapMin(k) = d.gapMin; end

    deltaE_k = out.deltaEHist(:,k);
    deltaL_k = out.deltaLHist(:,k);
    gap = deltaE_k - deltaL_k;
    [~, iMin] = min(gap);
    zContact = zGrid(iMin);

    deltaEofz = @(z) interp1(zGrid, deltaE_k, z, 'linear', 'extrap');
    slopeE = (deltaEofz(zContact+dz_fd) - deltaEofz(zContact-dz_fd)) / (2*dz_fd);
    slopeAtContact(k) = slopeE;
end

for k = 2:nSteps
    if isfinite(gapMin(k)) && isfinite(gapMin(k-1)) && gapMin(k-1) ~= 0
        relRate(k) = (gapMin(k)-gapMin(k-1)) / gapMin(k-1);
    end
end

fprintf('%6s | %10s | %10s | %10s\n', 'step', 'gap[um]', '|slope|', 'relRate%');
for k = 1:nSteps
    fprintf('%6d | %10.4f | %10.4f | %10.2f\n', k, gapMin(k)*1e6, abs(slopeAtContact(k)), relRate(k)*100);
end

fig = figure('Position',[100 100 900 600],'Color','w');
tl = tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');

ax1 = nexttile(tl);
plot(ax1, 1:nSteps, abs(slopeAtContact), 'o-'); grid(ax1,'on');
xlabel(ax1,'step'); ylabel(ax1,'|slope| at contact');
title(ax1,'Local interface slope at the min-gap location');

ax2 = nexttile(tl);
plot(ax2, 2:nSteps, relRate(2:end)*100, 'o-'); grid(ax2,'on');
xlabel(ax2,'step'); ylabel(ax2,'relative gap change [%/step]');
title(ax2,'Relative rate of gap closure');

sgtitle(fig,'Candidate explanations: slope and closure rate');
outPng = fullfile(fileparts(mfilename('fullpath')), 'slope_and_rate_check.png');
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('\nSaved: %s\nSMOKE_TEST_STATUS: OK\n', outPng);
