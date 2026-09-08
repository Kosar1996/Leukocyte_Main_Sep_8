%% CHECK_POISSON_INTERIOR_ENDOTHELIUM

% a point INSIDE the endothelium (not at the interface, not at the
% leukocyte's r=0 axis -- both already checked). Picks a genuine interior
% point: mid-wall radius, z away from both the domain edges (z=0,4) and
% the contact-loading region (z~3.6-3.9), so this reflects the material's
% own bulk elastic response, not a boundary or load-concentration effect.

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

fprintf('nuE = %.4f\n', par.nuE);

stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, par);
rE = out.meshE.nodes(:,1) + st.uE(1:2:end);
zE = out.meshE.nodes(:,2) + st.uE(2:2:end);
FrrE = scatteredInterpolant(rE, zE, stressE.sigma_rr, 'linear', 'nearest');
FttE = scatteredInterpolant(rE, zE, stressE.sigma_tt, 'linear', 'nearest');
FzzE = scatteredInterpolant(rE, zE, stressE.sigma_zz, 'linear', 'nearest');

% Interior point: z=2um (mid-domain, away from edges z=0/4 and the
% contact-loading region z~3.6-3.9), r = midway between inner wall and
% REout at that z.
zTest = 2.0e-6;
zGrid = out.z(:);
deltaE_k = out.deltaEHist(:,k);
rWallE = interp1(zGrid, deltaE_k, zTest, 'linear', 'extrap');
rInterior = 0.5*(rWallE + par.REout);

fprintf('\nInterior test point: r=%.4f um, z=%.4f um (rWall=%.4f, REout=%.4f)\n', ...
    rInterior*1e6, zTest*1e6, rWallE*1e6, par.REout*1e6);

srr = FrrE(rInterior, zTest);
stt = FttE(rInterior, zTest);
szz = FzzE(rInterior, zTest);

fprintf('sigma_rr = %.4f Pa\n', srr);
fprintf('sigma_tt = %.4f Pa\n', stt);
fprintf('sigma_zz = %.4f Pa\n', szz);

szz_predicted = par.nuE * (srr + stt);
fprintf('\nPlane-strain prediction: sigma_zz ~= nu*(sigma_rr+sigma_tt) = %.4f * %.4f = %.4f Pa\n', ...
    par.nuE, srr+stt, szz_predicted);
fprintf('Actual sigma_zz = %.4f Pa\n', szz);
fprintf('Difference = %.4f Pa (%.2f%% of predicted)\n', szz-szz_predicted, ...
    100*abs(szz-szz_predicted)/max(abs(szz_predicted),1e-9));

% Also scan a range of interior points (several z, several radial depths)
% to check whether this holds generally in the bulk, not just at one point
fprintf('\n=== Scan: several interior points (mid-wall radius, various z away from load/edges) ===\n');
fprintf('%8s %10s %10s %10s %10s %10s %12s\n', 'z[um]', 'r[um]', 'srr', 'stt', 'szz_actual', 'szz_pred', 'diff%%');
zScan = [1.0, 1.5, 2.0, 2.5, 3.0]*1e-6;
for i = 1:numel(zScan)
    zz = zScan(i);
    rWall = interp1(zGrid, deltaE_k, zz, 'linear', 'extrap');
    rMid = 0.5*(rWall + par.REout);
    srr_i = FrrE(rMid, zz);
    stt_i = FttE(rMid, zz);
    szz_i = FzzE(rMid, zz);
    szz_pred_i = par.nuE * (srr_i + stt_i);
    diffPct = 100*abs(szz_i-szz_pred_i)/max(abs(szz_pred_i),1e-9);
    fprintf('%8.3f %10.4f %10.4f %10.4f %10.4f %10.4f %12.2f\n', ...
        zz*1e6, rMid*1e6, srr_i, stt_i, szz_i, szz_pred_i, diffPct);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
