%% STEP2_FULL_COMPONENT_TRACE

% pressure components. Using the spatial derivatives identified in 2 to
% calculate the tangential and normal stresses and compare with 2, see
% how off they are."
%
% step2_four_point_trace.m already computes the raw sigma_rr, sigma_zz,
% sigma_rz tensors at each of the 4 points for both solid and fluid (to
% build the normal/tangential projection) but never PRINTS them -- only
% the final projected tn/tt values were reported. This script is the same
% exact computation, with the raw components and the normal-vector
% components explicitly recorded alongside the projected result, so the
% full derivation is visible, not just the end result.

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

mismatchOpts = struct('nQuery', 61, 'epsFrac', 0.1, 'trimFrac', 0.05);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);
refScaleL = leuko.stats.refScale;
refScaleE = endo.stats.refScale;

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

interpF = @(x,y,v) scatteredInterpolant(x(:), y(:), v(:), 'linear', 'none');
FfluidRR = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,1));
FfluidTT = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,2));
FfluidZZ = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,3));
FfluidRZ = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,4));
FfluidP  = interpF(centerFluid(:,1), centerFluid(:,2), pCellF);

FLrr = interpF(rL, zL, stressL.sigma_rr);
FLtt = interpF(rL, zL, stressL.sigma_tt);
FLzz = interpF(rL, zL, stressL.sigma_zz);
FLrz = interpF(rL, zL, stressL.sigma_rz);
FErr = interpF(rE, zE, stressE.sigma_rr);
FEtt = interpF(rE, zE, stressE.sigma_tt);
FEzz = interpF(rE, zE, stressE.sigma_zz);
FErz = interpF(rE, zE, stressE.sigma_rz);

mesh = fl.meshF;
deltaLofz = @(z) interp1(mesh.zc, mesh.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh.zc, mesh.deltaE_c, z, 'linear', 'extrap');
dz_fd = 1e-8;
slopeLofz = @(z) (deltaLofz(z+dz_fd) - deltaLofz(z-dz_fd)) / (2*dz_fd);
slopeEofz = @(z) (deltaEofz(z+dz_fd) - deltaEofz(z-dz_fd)) / (2*dz_fd);

epsFrac = mismatchOpts.epsFrac;
zPoints = [3.45e-6, 3.50e-6];

fprintf('================ LEUKOCYTE-FLUID INTERFACE -- FULL COMPONENT TRACE ================\n');
for i = 1:numel(zPoints)
    z = zPoints(i);
    rLwall = deltaLofz(z);
    gapHere = deltaEofz(z) - rLwall;
    eps_ = epsFrac * gapHere;
    slopeL = slopeLofz(z);
    nL = [1, -slopeL] / norm([1, -slopeL]);

    rq  = rLwall + eps_;
    rq2 = rLwall - eps_;

    fRR = FfluidRR(rq,z); fTT = FfluidTT(rq,z); fZZ = FfluidZZ(rq,z); fRZ = FfluidRZ(rq,z); fP = FfluidP(rq,z);
    sRR = FLrr(rq2,z); sTT = FLtt(rq2,z); sZZ = FLzz(rq2,z); sRZ = FLrz(rq2,z);

    sigF = [fRR, fRZ; fRZ, fZZ];
    sigS = [sRR, sRZ; sRZ, sZZ];
    [tnF, ttF] = traction_components_local(sigF, nL);
    [tnS, ttS] = traction_components_local(sigS, nL);
    pctEn = 100*abs(tnS-tnF)/refScaleL;
    pctEt = 100*abs(ttS-ttF)/refScaleL;

    fprintf('\n--- Point L%d: z = %.4f um, wall r = %.4f um, slope = %.4f, normal n = (%.4f, %.4f) ---\n', ...
        i, z*1e6, rLwall*1e6, slopeL, nL(1), nL(2));
    fprintf('  Sample coords: fluid (r=%.4f, z=%.4f) um,  solid (r=%.4f, z=%.4f) um\n', rq*1e6, z*1e6, rq2*1e6, z*1e6);
    fprintf('  FLUID raw components:  sigma_rr=%.4f  sigma_tt=%.4f  sigma_zz=%.4f  sigma_rz=%.4f  P=%.4f  [Pa]\n', fRR, fTT, fZZ, fRZ, fP);
    fprintf('  SOLID raw components:  sigma_rr=%.4f  sigma_tt=%.4f  sigma_zz=%.4f  sigma_rz=%.4f  [Pa]\n', sRR, sTT, sZZ, sRZ);
    fprintf('  Projected normal:   solid=%.4f  fluid=%.4f  (diff=%.4f, %.2f%% of refScale)\n', tnS, tnF, tnS-tnF, pctEn);
    fprintf('  Projected tangent:  solid=%.4f  fluid=%.4f  (diff=%.4f, %.2f%% of refScale)\n', ttS, ttF, ttS-ttF, pctEt);
end

fprintf('\n================ ENDOTHELIUM-FLUID INTERFACE -- FULL COMPONENT TRACE ================\n');
for i = 1:numel(zPoints)
    z = zPoints(i);
    rEwall = deltaEofz(z);
    gapHere = rEwall - deltaLofz(z);
    eps_ = epsFrac * gapHere;
    slopeE = slopeEofz(z);
    nE = [1, -slopeE] / norm([1, -slopeE]);

    rq  = rEwall - eps_;
    rq2 = rEwall + eps_;

    fRR = FfluidRR(rq,z); fTT = FfluidTT(rq,z); fZZ = FfluidZZ(rq,z); fRZ = FfluidRZ(rq,z); fP = FfluidP(rq,z);
    sRR = FErr(rq2,z); sTT = FEtt(rq2,z); sZZ = FEzz(rq2,z); sRZ = FErz(rq2,z);

    sigF = [fRR, fRZ; fRZ, fZZ];
    sigS = [sRR, sRZ; sRZ, sZZ];
    [tnF, ttF] = traction_components_local(sigF, nE);
    [tnS, ttS] = traction_components_local(sigS, nE);
    pctEn = 100*abs(tnS-tnF)/refScaleE;
    pctEt = 100*abs(ttS-ttF)/refScaleE;

    fprintf('\n--- Point E%d: z = %.4f um, wall r = %.4f um, slope = %.4f, normal n = (%.4f, %.4f) ---\n', ...
        i, z*1e6, rEwall*1e6, slopeE, nE(1), nE(2));
    fprintf('  Sample coords: fluid (r=%.4f, z=%.4f) um,  solid (r=%.4f, z=%.4f) um\n', rq*1e6, z*1e6, rq2*1e6, z*1e6);
    fprintf('  FLUID raw components:  sigma_rr=%.4f  sigma_tt=%.4f  sigma_zz=%.4f  sigma_rz=%.4f  P=%.4f  [Pa]\n', fRR, fTT, fZZ, fRZ, fP);
    fprintf('  SOLID raw components:  sigma_rr=%.4f  sigma_tt=%.4f  sigma_zz=%.4f  sigma_rz=%.4f  [Pa]\n', sRR, sTT, sZZ, sRZ);
    fprintf('  Projected normal:   solid=%.4f  fluid=%.4f  (diff=%.4f, %.2f%% of refScale)\n', tnS, tnF, tnS-tnF, pctEn);
    fprintf('  Projected tangent:  solid=%.4f  fluid=%.4f  (diff=%.4f, %.2f%% of refScale)\n', ttS, ttF, ttS-ttF, pctEt);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');

function [tn, tt] = traction_components_local(sig, n)
n = n(:);
that = [-n(2); n(1)];
t = sig * n;
tn = t.' * n;
tt = t.' * that;
end
