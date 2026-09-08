%% CHECK_LEUKOCYTE_FLUID_EXPECTED_TRENDS

% explicit request for the "equivalent" exercise on all three.
%
% Leukocyte: solid cylinder (no hole) under external pressure. Enforcing
% regularity at r=0 (no B/r^2 singularity term) gives the classical
% result sigma_rr = sigma_theta_theta, UNIFORM across the WHOLE
% cross-section, not just at r=0. Testing that here.
%
% Fluid: thin-film lubrication theory assumes pressure doesn't vary
% meaningfully ACROSS the gap (r-direction) at a given z -- dP/dr ~ 0
% within the gap. Testing that directly.

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

%% Leukocyte: sigma_rr vs sigma_theta_theta across the whole radius
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
rMaxL = max(rL(abs(zL - zMidL) < 0.3e-6)); % leukocyte's own outer radius near midline
rQuery = linspace(0, rMaxL*0.95, 12)';

fprintf('=== LEUKOCYTE: sigma_rr vs sigma_theta_theta across the radius (z=midline=%.3f um) ===\n', zMidL*1e6);
fprintf('%8s | %10s %10s | %8s\n', 'r[um]', 'sigma_rr', 'sigma_tt', 'diff%%');
for i = 1:numel(rQuery)
    srr = FrrL(rQuery(i), zMidL);
    stt = FttL(rQuery(i), zMidL);
    fprintf('%8.3f | %10.4f %10.4f | %8.2f\n', rQuery(i)*1e6, srr, stt, 100*abs(srr-stt)/max(abs(stt),1e-9));
end

%% Fluid: pressure across the gap (r-direction) at fixed z, should be ~flat
meshF2 = add_fluid_nodes(fl.meshF);
[pCellF, ~, centerF] = recover_fluid_nodes_pressure_stress_Q4(meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell);
FpFl = scatteredInterpolant(centerF(:,1), centerF(:,2), pCellF, 'linear', 'nearest');

zGrid = out.z(:);
deltaE_k = out.deltaEHist(:,k);
deltaL_k = out.deltaLHist(:,k);
zTest = zMidL; % check at the load center, where pressure is largest and most meaningful
rEwall = interp1(zGrid, deltaE_k, zTest, 'linear', 'extrap');
rLwall = interp1(zGrid, deltaL_k, zTest, 'linear', 'extrap');
rAcrossGap = linspace(rLwall*1.02, rEwall*0.98, 10)';

fprintf('\n=== FLUID: pressure across the gap at z=%.3f um (rLwall=%.3f, rEwall=%.3f um) ===\n', ...
    zTest*1e6, rLwall*1e6, rEwall*1e6);
fprintf('%8s | %10s\n', 'r[um]', 'P[Pa]');
Pvals = nan(size(rAcrossGap));
for i = 1:numel(rAcrossGap)
    Pvals(i) = FpFl(rAcrossGap(i), zTest);
    fprintf('%8.4f | %10.4f\n', rAcrossGap(i)*1e6, Pvals(i));
end
fprintf('Spread across gap: %.4f Pa (%.2f%% of mean)\n', ...
    max(Pvals)-min(Pvals), 100*(max(Pvals)-min(Pvals))/max(abs(mean(Pvals)),1e-9));

fprintf('\nSMOKE_TEST_STATUS: OK\n');
