%% CHECK_POISSON_AXIAL_FORCE_BALANCE_WIDEDOMAIN
% Same net-axial-force check as check_poisson_axial_force_balance.m, same
% z test points, but on the wide-domain (z=[-2,6]um) run from earlier
% tonight -- tests directly whether extending the domain reduces the
% residual force-balance imbalance found on the narrow domain, which
% would confirm (not just plausibly connect) the short-domain hypothesis.

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/code/leukocyte-main/out_2D_wide_domain_40step.mat';
S = load(file);
out = S.out;

k = out.stopStep;
st = out.stateHist{k};
par = out.par;
dtStep = out.dtHist(k);

fprintf('Wide-domain run: par.zMin=%.2f, par.zMax=%.2f um (meshE actual range: [%.4f, %.4f] um)\n', ...
    par.zMin*1e6, par.zMax*1e6, min(out.meshE.nodes(:,2))*1e6, max(out.meshE.nodes(:,2))*1e6);

stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, par);
rE = out.meshE.nodes(:,1) + st.uE(1:2:end);
zE = out.meshE.nodes(:,2) + st.uE(2:2:end);
FzzE = scatteredInterpolant(rE, zE, stressE.sigma_zz, 'linear', 'nearest');

zGrid = out.z(:);
deltaE_k = out.deltaEHist(:,k);

% Same absolute z test points as the narrow-domain check, for direct comparison
zTestList = [1.0, 1.5, 2.0, 2.5, 3.0] * 1e-6;
fprintf('\n%8s %16s %16s %16s %10s\n', 'z[um]', 'netAxialForce[N]', 'refForceScale[N]', 'ratio', 'verdict');
for i = 1:numel(zTestList)
    zz = zTestList(i);
    rWall = interp1(zGrid, deltaE_k, zz, 'linear', 'extrap');
    rQuery = linspace(rWall, par.REout, 400)';

    szz = arrayfun(@(rr) FzzE(rr, zz), rQuery);
    integrand = szz .* (2*pi*rQuery);
    netForce = trapz(rQuery, integrand);
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

fprintf('\nCompare directly against the narrow-domain result at the same z values:\n');
fprintf('  narrow: z=1.0->17.0%%, z=1.5->7.7%%, z=2.0->5.4%%, z=2.5->8.0%%, z=3.0->17.3%%\n');
fprintf('SMOKE_TEST_STATUS: OK\n');
