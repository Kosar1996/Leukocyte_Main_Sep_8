%% PLOT_EXPECTED_TRENDS_ALL_THREE

% the leukocyte, as well as for the fluid as well." Produces matching
% sigma-vs-r trend plots for all three domains, 2D solver, last saved step
% (out_2D_t10_for_review.mat), consistent style so they can go side by
% side for Wednesday.
%
% Endothelium: sigma_rr, sigma_tt vs r at the load-adjacent z (contact
%   zone). Expect sigma_rr+sigma_tt ~ constant (Lame), sigma_rr -> ~0-ish
%   trend toward r=REout (roller: radially fixed, not stress-free, so
%   sigma_rr does not have to hit exactly 0 there -- shown for reference).
%
% Leukocyte: sigma_rr, sigma_tt vs r at the midline z. Expect sigma_rr =
%   sigma_tt at r=0 (regularity), both roughly uniform near center,
%   diverging toward the loaded outer surface.
%
% Fluid: pressure and sigma_rz vs z, spanning from inside the gap out to
%   the domain edges. Expect decay toward 0 outside the pore (no applied
%   far-field pressure).

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

k = out.stopStep;
st = out.stateHist{k};
fl = out.fluidHist{k};
par = out.par;
dtStep = out.dtHist(k);

zGrid = out.z(:);
deltaE_k = out.deltaEHist(:,k);
deltaL_k = out.deltaLHist(:,k);
gap = deltaE_k - deltaL_k;
[~, iMin] = min(gap);
zContact = zGrid(iMin);

%% Endothelium: sigma_rr, sigma_tt vs r at the contact-adjacent z
stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, par);
rE = out.meshE.nodes(:,1) + st.uE(1:2:end);
zE = out.meshE.nodes(:,2) + st.uE(2:2:end);
FrrE = scatteredInterpolant(rE, zE, stressE.sigma_rr, 'linear', 'nearest');
FttE = scatteredInterpolant(rE, zE, stressE.sigma_tt, 'linear', 'nearest');

rWallE = interp1(zGrid, deltaE_k, zContact, 'linear', 'extrap');
rMaxE = max(rE(abs(zE - zContact) < 0.3e-6));
rQueryE = linspace(rWallE*1.001, rMaxE*0.999, 60)';
srrE = arrayfun(@(rr) FrrE(rr, zContact), rQueryE);
sttE = arrayfun(@(rr) FttE(rr, zContact), rQueryE);

%% Leukocyte: sigma_rr, sigma_tt vs r at midline
parL = par;
if isfield(par, 'GL'),   parL.Ge   = par.GL;   end
if isfield(par, 'KL'),   parL.Ke   = par.KL;   end
if isfield(par, 'etaL'), parL.etaE = par.etaL; end
stressL = recover_nodal_stress_axisym_viscoelastic(out.meshL, st.uL, st.uLPrev, dtStep, parL);
rL = out.meshL.nodes(:,1) + st.uL(1:2:end);
zL = out.meshL.nodes(:,2) + st.uL(2:2:end);
FrrL = scatteredInterpolant(rL, zL, stressL.sigma_rr, 'linear', 'nearest');
FttL = scatteredInterpolant(rL, zL, stressL.sigma_tt, 'linear', 'nearest');

zMidL = 0.5*(min(zL) + max(zL));
rMaxL = max(rL(abs(zL - zMidL) < 0.3e-6));
rQueryL = linspace(0, rMaxL*0.98, 60)';
srrL = arrayfun(@(rr) FrrL(rr, zMidL), rQueryL);
sttL = arrayfun(@(rr) FttL(rr, zMidL), rQueryL);

%% Fluid: pressure and sigma_rz vs z, gap centerline out to domain edges
meshF2 = add_fluid_nodes(fl.meshF);
[pCellF, sigmaCellF, centerF] = recover_fluid_nodes_pressure_stress_Q4(meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell);
FpFl = scatteredInterpolant(centerF(:,1), centerF(:,2), pCellF, 'linear', 'nearest');
FrzFl = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,4), 'linear', 'nearest');

zQueryF = linspace(min(zGrid), max(zGrid), 120)';
rMidGap = 0.5*(interp1(zGrid, deltaE_k, zQueryF, 'linear','extrap') + ...
               interp1(zGrid, deltaL_k, zQueryF, 'linear','extrap'));
pAlongGap = arrayfun(@(i) FpFl(rMidGap(i), zQueryF(i)), (1:numel(zQueryF))');
rzAlongGap = arrayfun(@(i) FrzFl(rMidGap(i), zQueryF(i)), (1:numel(zQueryF))');

%% Plot
fig = figure('Position',[100 100 1200 800],'Color','w');
tl = tiledlayout(fig,3,2,'Padding','compact','TileSpacing','compact');
title(tl, 'Expected trends -- endothelium, leukocyte, fluid (2D solver, step 40)');

ax1 = nexttile(tl);
plot(ax1, rQueryE*1e6, srrE, 'o-', 'Color',[0.85 0.2 0.2], 'DisplayName','\sigma_{rr}'); hold(ax1,'on');
plot(ax1, rQueryE*1e6, sttE, 's-', 'Color',[0.2 0.4 0.85], 'DisplayName','\sigma_{\theta\theta}');
plot(ax1, rQueryE*1e6, srrE+sttE, '--', 'Color',[0.3 0.3 0.3], 'DisplayName','\sigma_{rr}+\sigma_{\theta\theta}');
grid(ax1,'on'); legend(ax1,'Location','best'); xlabel(ax1,'r [\mum]'); ylabel(ax1,'Pa');
title(ax1, sprintf('Endothelium at z=%.2f\\mum (contact zone)', zContact*1e6));

ax2 = nexttile(tl);
plot(ax2, rQueryL*1e6, srrL, 'o-', 'Color',[0.85 0.2 0.2], 'DisplayName','\sigma_{rr}'); hold(ax2,'on');
plot(ax2, rQueryL*1e6, sttL, 's-', 'Color',[0.2 0.4 0.85], 'DisplayName','\sigma_{\theta\theta}');
grid(ax2,'on'); legend(ax2,'Location','best'); xlabel(ax2,'r [\mum]'); ylabel(ax2,'Pa');
title(ax2, sprintf('Leukocyte at z=%.2f\\mum (midline) -- expect \\sigma_{rr}=\\sigma_{\\theta\\theta} at r=0', zMidL*1e6));

ax3 = nexttile(tl,[1 2]);
yyaxis(ax3,'left');
plot(ax3, zQueryF*1e6, pAlongGap, 'o-', 'Color',[0.15 0.6 0.3]);
ylabel(ax3, 'P [Pa]');
yyaxis(ax3,'right');
plot(ax3, zQueryF*1e6, rzAlongGap, '^-', 'Color',[0.6 0.3 0.7]);
ylabel(ax3, '\sigma_{rz} [Pa]');
grid(ax3,'on'); xlabel(ax3,'z [\mum]');
title(ax3, 'Fluid, along gap centerline -- expect decay to 0 outside the pore');
xline(ax3, min(zL)*1e6, ':', 'leuko z_{min}', 'Color',[0.5 0.5 0.5]);
xline(ax3, max(zL)*1e6, ':', 'leuko z_{max}', 'Color',[0.5 0.5 0.5]);

ax4 = nexttile(tl);
plot(ax4, rQueryE*1e6, srrE, 'o-', 'Color',[0.85 0.2 0.2]);
grid(ax4,'on'); xlabel(ax4,'r [\mum]'); ylabel(ax4,'Pa');
title(ax4, sprintf('Endothelium \\sigma_{rr} only -- REout=%.1f\\mum', par.REout*1e6));

ax5 = nexttile(tl);
diffPct = 100*abs(srrL-sttL)./max(abs(sttL),1e-9);
plot(ax5, rQueryL*1e6, diffPct, 'o-', 'Color',[0.5 0.3 0.1]);
grid(ax5,'on'); xlabel(ax5,'r [\mum]'); ylabel(ax5,'%% diff');
title(ax5, '|\sigma_{rr}-\sigma_{\theta\theta}|/\sigma_{\theta\theta}, leukocyte');

outPng = fullfile(fileparts(mfilename('fullpath')), 'expected_trends_all_three.png');
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('Saved: %s\n', outPng);
fprintf('At leukocyte r=0: sigma_rr=%.4f, sigma_tt=%.4f (diff=%.3f%%)\n', srrL(1), sttL(1), 100*abs(srrL(1)-sttL(1))/max(abs(sttL(1)),1e-9));
fprintf('Fluid P at z far outside pore (z=%.2fum): %.4f Pa; at z=%.2fum: %.4f Pa\n', ...
    zQueryF(1)*1e6, pAlongGap(1), zQueryF(end)*1e6, pAlongGap(end));
fprintf('SMOKE_TEST_STATUS: OK\n');
