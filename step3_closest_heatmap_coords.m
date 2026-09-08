%% STEP3_CLOSEST_HEATMAP_COORDS

% heatmap matrix the closest heatmap coordinates... the closest heatmap
% coordinate may be on the solid side or the fluid side. Record the
% coordinates as well as all stress and pressure components. If on the
% heatmap it shows as a fluid point but its coordinate in comparison to
% the interface coordinates suggest it should be a solid point. Alert.
% Using the spatial derivatives identified in 2 to calculate the
% tangential and normal stresses and compare with 2, see how off they are"
%
% For each of the 4 step-2 points (L1,L2 on leuko-fluid interface,
% E1,E2 on endo-fluid interface), finds the nearest ACTUAL computational
% grid point on both the fluid side (MAC cell centers) and the solid side
% (FE nodes), reports which is truly closest overall, flags any
% solid/fluid classification confusion, and recomputes tangential/normal
% stress at that nearest real grid point (same slope-based normal vector
% as step 2) to see how much it differs from step 2's interpolated value.

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

parLmismatch = par;
if isfield(par, 'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par, 'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par, 'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

mu = par.mu;
meshF = add_fluid_nodes(fl.meshF);
[pCellF, sigmaCellFluid, centerFluid] = recover_fluid_nodes_pressure_stress_Q4( ...
    meshF, fl.ur2D, fl.uz2D, mu, fl.pCell);

stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, par);
stressL = recover_nodal_stress_axisym_viscoelastic(out.meshL, st.uL, st.uLPrev, dtStep, parLmismatch);

rE = out.meshE.nodes(:,1) + st.uE(1:2:end);
zE = out.meshE.nodes(:,2) + st.uE(2:2:end);
rL = out.meshL.nodes(:,1) + st.uL(1:2:end);
zL = out.meshL.nodes(:,2) + st.uL(2:2:end);

mesh = fl.meshF;
deltaLofz = @(z) interp1(mesh.zc, mesh.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh.zc, mesh.deltaE_c, z, 'linear', 'extrap');
dz_fd = 1e-8;
slopeLofz = @(z) (deltaLofz(z+dz_fd) - deltaLofz(z-dz_fd)) / (2*dz_fd);
slopeEofz = @(z) (deltaEofz(z+dz_fd) - deltaEofz(z-dz_fd)) / (2*dz_fd);

zPoints = [3.45e-6, 3.50e-6];

function report_point(label, rWall, z, isLeuko, rE, zE, rL, zL, centerFluid, ...
        sigmaCellFluid, pCellF, stressE, stressL, nvec)
    % nearest FLUID grid point
    df = hypot(centerFluid(:,1)-rWall, centerFluid(:,2)-z);
    [dfMin, idxF] = min(df);
    % nearest SOLID grid point (on the relevant solid for this interface)
    if isLeuko
        rS = rL; zS = zL; stressS = stressL;
    else
        rS = rE; zS = zE; stressS = stressE;
    end
    ds = hypot(rS-rWall, zS-z);
    [dsMin, idxS] = min(ds);

    fprintf('\n--- %s: wall/interface point (r=%.4f, z=%.4f) um ---\n', label, rWall*1e6, z*1e6);
    fprintf('  Nearest FLUID grid point: (r=%.4f, z=%.4f) um, distance=%.2f nm\n', ...
        centerFluid(idxF,1)*1e6, centerFluid(idxF,2)*1e6, dfMin*1e9);
    fprintf('    sigma_rr=%.4f, sigma_zz=%.4f, sigma_rz=%.4f Pa, P=%.4f Pa\n', ...
        sigmaCellFluid(idxF,1), sigmaCellFluid(idxF,3), sigmaCellFluid(idxF,4), pCellF(idxF));
    fprintf('  Nearest SOLID grid point: (r=%.4f, z=%.4f) um, distance=%.2f nm\n', ...
        rS(idxS)*1e6, zS(idxS)*1e6, dsMin*1e9);
    fprintf('    sigma_rr=%.4f, sigma_zz=%.4f, sigma_rz=%.4f Pa\n', ...
        stressS.sigma_rr(idxS), stressS.sigma_zz(idxS), stressS.sigma_rz(idxS));

    % Classification check: which is actually closer overall?
    if dfMin < dsMin
        fprintf('  --> Overall closest grid point is on the FLUID side (as expected for an interface point).\n');
    else
        fprintf('  --> ALERT: overall closest grid point is on the SOLID side, even though we are checking the fluid interface.\n');
    end

    % Recompute tangential/normal stress AT the nearest fluid grid node
    % and AT the nearest solid grid node, using the SAME normal vector as step 2
    that = [-nvec(2); nvec(1)];
    sigF_near = [sigmaCellFluid(idxF,1), sigmaCellFluid(idxF,4); sigmaCellFluid(idxF,4), sigmaCellFluid(idxF,3)];
    tF = sigF_near * nvec(:); tnF_near = tF.'*nvec(:); ttF_near = tF.'*that;
    sigS_near = [stressS.sigma_rr(idxS), stressS.sigma_rz(idxS); stressS.sigma_rz(idxS), stressS.sigma_zz(idxS)];
    tS = sigS_near * nvec(:); tnS_near = tS.'*nvec(:); ttS_near = tS.'*that;
    fprintf('  At nearest grid points -- sigma_normal: solid=%.4f, fluid=%.4f Pa (diff=%.4f)\n', tnS_near, tnF_near, tnS_near-tnF_near);
    fprintf('  At nearest grid points -- sigma_tangent: solid=%.4f, fluid=%.4f Pa (diff=%.4f)\n', ttS_near, ttF_near, ttS_near-ttF_near);
end

fprintf('================ LEUKOCYTE-FLUID INTERFACE ================\n');
for i = 1:numel(zPoints)
    z = zPoints(i);
    rWall = deltaLofz(z);
    slopeL = slopeLofz(z);
    nvec = [1, -slopeL] / norm([1, -slopeL]);
    report_point(sprintf('L%d',i), rWall, z, true, rE, zE, rL, zL, centerFluid, sigmaCellFluid, pCellF, stressE, stressL, nvec);
end

fprintf('\n================ ENDOTHELIUM-FLUID INTERFACE ================\n');
for i = 1:numel(zPoints)
    z = zPoints(i);
    rWall = deltaEofz(z);
    slopeE = slopeEofz(z);
    nvec = [1, -slopeE] / norm([1, -slopeE]);
    report_point(sprintf('E%d',i), rWall, z, false, rE, zE, rL, zL, centerFluid, sigmaCellFluid, pCellF, stressE, stressL, nvec);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
