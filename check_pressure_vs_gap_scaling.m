%% CHECK_PRESSURE_VS_GAP_SCALING
% Tests whether the mismatch blowup at small gap (steps 37-40) is
% consistent with classical lubrication scaling: fluid pressure grows
% steeply (roughly ~1/h^n) as the gap h shrinks. If pressure is genuinely
% blowing up faster than the correction (fixed relax, fixed pass count)
% can track, that's a physical explanation for why mismatch grows near
% contact -- not evidence of a separate bug.

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

par = out.par;
stepsToCheck = 1:out.stopStep;

gapMin = nan(numel(stepsToCheck),1);
pMaxAtContact = nan(numel(stepsToCheck),1);

for si = 1:numel(stepsToCheck)
    k = stepsToCheck(si);
    d = out.diagHist{k};
    if ~isempty(d) && isfield(d,'gapMin'), gapMin(si) = d.gapMin; end

    fl = out.fluidHist{k};
    zGrid = out.z(:);
    deltaE_k = out.deltaEHist(:,k);
    deltaL_k = out.deltaLHist(:,k);
    gap = deltaE_k - deltaL_k;
    [~, iMin] = min(gap);

    meshF = fl.meshF;
    zc = meshF.zc(:);
    P = fl.P;
    P_wall_raw = P(end,:).';
    zContact = zGrid(iMin);
    pMaxAtContact(si) = interp1(zc, P_wall_raw, zContact, 'linear', 'extrap');
end

fprintf('%6s | %10s | %12s\n', 'step', 'gap[um]', 'P_contact[Pa]');
for si = 1:numel(stepsToCheck)
    fprintf('%6d | %10.4f | %12.2f\n', stepsToCheck(si), gapMin(si)*1e6, pMaxAtContact(si));
end

% Fit log(P) vs log(gap) over the later steps where gap is monotonically shrinking (steps 17-40)
idxFit = find(stepsToCheck >= 17);
logGap = log(gapMin(idxFit));
logP = log(abs(pMaxAtContact(idxFit)));
validFit = isfinite(logGap) & isfinite(logP);
pfit = polyfit(logGap(validFit), logP(validFit), 1);
fprintf('\nPower-law fit over steps 17-40 (gap shrinking phase): P ~ gap^%.2f\n', pfit(1));
fprintf('(Classical lubrication squeeze-film scaling would predict something like gap^-1 to gap^-3\n');
fprintf('depending on geometry/flow regime -- a strongly negative exponent here supports the\n');
fprintf('"pressure genuinely steepens near contact" explanation.)\n');

fig = figure('Position',[100 100 800 500],'Color','w');
loglog(gapMin(idxFit)*1e6, abs(pMaxAtContact(idxFit)), 'o-');
xlabel('min gap [\mum]'); ylabel('|P| at contact [Pa]'); grid on;
title(sprintf('Pressure vs gap (log-log), steps 17-40 -- fitted slope = %.2f', pfit(1)));
outPng = fullfile(fileparts(mfilename('fullpath')), 'pressure_vs_gap_scaling.png');
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('\nSaved: %s\nSMOKE_TEST_STATUS: OK\n', outPng);
