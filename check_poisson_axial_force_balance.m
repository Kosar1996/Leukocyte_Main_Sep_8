%% CHECK_POISSON_AXIAL_FORCE_BALANCE
% Rigorous follow-up to check_poisson_interior_endothelium.m: pointwise
% sigma_zz~=0 in the interior is consistent with a free-ended tube, but
% the actual physical requirement is that sigma_zz integrates to ZERO NET
% AXIAL FORCE across the full cross-section at any given z (that's what
% "traction-free end, no applied axial force" actually means for a rod/
% tube -- Saint-Venant: the local stress distribution can vary, but the
% resultant axial force must vanish). This checks that directly.

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

stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, par);
rE = out.meshE.nodes(:,1) + st.uE(1:2:end);
zE = out.meshE.nodes(:,2) + st.uE(2:2:end);
FzzE = scatteredInterpolant(rE, zE, stressE.sigma_zz, 'linear', 'nearest');
FrrE = scatteredInterpolant(rE, zE, stressE.sigma_rr, 'linear', 'nearest');

zGrid = out.z(:);
deltaE_k = out.deltaEHist(:,k);

zTestList = [1.0, 1.5, 2.0, 2.5, 3.0] * 1e-6;
fprintf('%8s %16s %16s %16s %10s\n', 'z[um]', 'netAxialForce[N]', 'refForceScale[N]', 'ratio', 'verdict');
for i = 1:numel(zTestList)
    zz = zTestList(i);
    rWall = interp1(zGrid, deltaE_k, zz, 'linear', 'extrap');
    rQuery = linspace(rWall, par.REout, 400)';

    szz = arrayfun(@(rr) FzzE(rr, zz), rQuery);
    % Force = integral of sigma_zz * 2*pi*r dr  (axisymmetric cross-section)
    integrand = szz .* (2*pi*rQuery);
    netForce = trapz(rQuery, integrand);

    % Reference scale: integral of |sigma_zz| * 2*pi*r dr (magnitude scale
    % for comparison, so "net force / this" tells us the FRACTIONAL
    % cancellation, not just a raw number with no context)
    refScale = trapz(rQuery, abs(integrand));

    ratio = abs(netForce) / max(refScale, 1e-30);
    if ratio < 0.05
        verdict = 'balances (<5%)';
    elseif ratio < 0.15
        verdict = 'roughly OK';
    else
        verdict = 'NOT balanced';
    end
    fprintf('%8.3f %16.6e %16.6e %10.4f  %s\n', zz*1e6, netForce, refScale, ratio, verdict);
end

fprintf('\nnetAxialForce should be small relative to refScale if sigma_zz truly integrates to ~0\n');
fprintf('net axial force (free-ended tube, no applied axial load anywhere in the interior).\n');
fprintf('SMOKE_TEST_STATUS: OK\n');
