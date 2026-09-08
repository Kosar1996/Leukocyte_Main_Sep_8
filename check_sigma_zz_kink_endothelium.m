clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

k = out.stopStep;
st = out.stateHist{k};
par = out.par;
dtStep = out.dtHist(k);

zGrid = out.z(:);
deltaE_k = out.deltaEHist(:,k);
deltaL_k = out.deltaLHist(:,k);
gap = deltaE_k - deltaL_k;
[~, iMin] = min(gap);
zContact = zGrid(iMin);

stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, par);
rE = out.meshE.nodes(:,1) + st.uE(1:2:end);
zE = out.meshE.nodes(:,2) + st.uE(2:2:end);
FzzE = scatteredInterpolant(rE, zE, stressE.sigma_zz, 'linear', 'nearest');

rWallE = interp1(zGrid, deltaE_k, zContact, 'linear', 'extrap');
rMaxE = max(rE(abs(zE - zContact) < 0.3e-6));
rQuery = linspace(rWallE*1.001, rMaxE*0.999, 100)';
szz = arrayfun(@(rr) FzzE(rr, zContact), rQuery);

fprintf('%8s | %10s\n', 'r[um]', 'sigma_zz');
for i = 1:5:numel(rQuery)
    fprintf('%8.3f | %10.4f\n', rQuery(i)*1e6, szz(i));
end
fprintf('...\n');
for i = max(1,numel(rQuery)-15):numel(rQuery)
    fprintf('%8.3f | %10.4f\n', rQuery(i)*1e6, szz(i));
end

fig = figure('Position',[100 100 800 500],'Color','w');
plot(rQuery*1e6, szz, 'o-', 'Color',[0.7 0.2 0.5]);
grid on; xlabel('r [\mum]'); ylabel('\sigma_{zz} [Pa]');
title(sprintf('Endothelium \\sigma_{zz} vs r at z=%.2f\\mum (contact zone)', zContact*1e6));
yline(0,'k--');
outPng = fullfile(fileparts(mfilename('fullpath')), 'sigma_zz_kink_check.png');
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('\nSaved: %s\n', outPng);
fprintf('SMOKE_TEST_STATUS: OK\n');
