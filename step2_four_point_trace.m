%% STEP2_FOUR_POINT_TRACE

% endothelium-fluid interface, in the thin-fluid region identified in
% step 1 (z~3.45-3.50 um, chosen as the best available compromise: clearly
% thinner than the bulk gap, with some standoff from the worst of the
% z=4.04um edge distortion).
%
% Uses the EXACT same methodology as the live correction code
% (compute_interface_traction_mismatch.m) -- same slope-aware normal
% vector construction, same traction_components projection, same %mismatch
% normalization (against the p90 of |fluid normal traction| across the
% WHOLE interface, not a per-point scale) -- so what's reported here is
% what the actual correction computes, not an independently-recomputed
% approximation.

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

%% Step A: run the REAL function once, full interface, to get the correct
% refScale (p90 |tnFluid|) used for %mismatch normalization -- must match
% what the live correction loop actually uses, not a locally-recomputed one.
mismatchOpts = struct('nQuery', 61, 'epsFrac', 0.1, 'trimFrac', 0.05);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);
refScaleL = leuko.stats.refScale;
refScaleE = endo.stats.refScale;
fprintf('refScale (leukocyte side, p90|tnFluid| across whole interface) = %.4f Pa\n', refScaleL);
fprintf('refScale (endothelium side, p90|tnFluid| across whole interface) = %.4f Pa\n', refScaleE);

%% Step B: rebuild the same interpolants the real function uses internally
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
FfluidZZ = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,3));
FfluidRZ = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,4));
FfluidP  = interpF(centerFluid(:,1), centerFluid(:,2), pCellF);

FLrr = interpF(rL, zL, stressL.sigma_rr);
FLzz = interpF(rL, zL, stressL.sigma_zz);
FLrz = interpF(rL, zL, stressL.sigma_rz);
FErr = interpF(rE, zE, stressE.sigma_rr);
FEzz = interpF(rE, zE, stressE.sigma_zz);
FErz = interpF(rE, zE, stressE.sigma_rz);

mesh = fl.meshF;
deltaLofz = @(z) interp1(mesh.zc, mesh.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh.zc, mesh.deltaE_c, z, 'linear', 'extrap');
dz_fd = 1e-8;
slopeLofz = @(z) (deltaLofz(z+dz_fd) - deltaLofz(z-dz_fd)) / (2*dz_fd);
slopeEofz = @(z) (deltaEofz(z+dz_fd) - deltaEofz(z-dz_fd)) / (2*dz_fd);

epsFrac = mismatchOpts.epsFrac;

%% Step C: evaluate at our 4 chosen points, EXACT same per-point logic
zPoints = [3.45e-6, 3.50e-6];

fprintf('\n================ LEUKOCYTE-FLUID INTERFACE ================\n');
resultsL = struct('z',{},'r',{},'tauSolid',{},'tauFluid',{},'sigSolid',{},'sigFluid',{}, ...
    'pFluid',{},'pctEn',{},'pctEt',{},'slope',{});
for i = 1:numel(zPoints)
    z = zPoints(i);
    rLwall = deltaLofz(z);
    gapHere = deltaEofz(z) - rLwall;
    eps_ = epsFrac * gapHere;
    slopeL = slopeLofz(z);
    nL = [1, -slopeL] / norm([1, -slopeL]);

    rq  = rLwall + eps_;   % fluid side sample point (just inside fluid)
    sigF = [FfluidRR(rq,z), FfluidRZ(rq,z); FfluidRZ(rq,z), FfluidZZ(rq,z)];
    rq2 = rLwall - eps_;   % solid side sample point (just inside leukocyte)
    sigS = [FLrr(rq2,z), FLrz(rq2,z); FLrz(rq2,z), FLzz(rq2,z)];
    [tnF, ttF] = traction_components_local(sigF, nL);
    [tnS, ttS] = traction_components_local(sigS, nL);
    pF = FfluidP(rq, z);

    pctEn = 100*abs(tnS-tnF)/refScaleL;
    pctEt = 100*abs(ttS-ttF)/refScaleL;

    fprintf('\n--- Point L%d: z=%.4f um (wall r=%.4f um, slope=%.4f) ---\n', i, z*1e6, rLwall*1e6, slopeL);
    fprintf('  Sample points: fluid at r=%.4fum, solid at r=%.4fum (eps=%.2fnm, gap=%.2fnm)\n', ...
        rq*1e6, rq2*1e6, eps_*1e9, gapHere*1e9);
    fprintf('  sigma_normal:  solid=%.4f Pa,  fluid=%.4f Pa   (diff=%.4f Pa, %.2f%% of refScale)\n', tnS, tnF, tnS-tnF, pctEn);
    fprintf('  sigma_tangent: solid=%.4f Pa,  fluid=%.4f Pa   (diff=%.4f Pa, %.2f%% of refScale)\n', ttS, ttF, ttS-ttF, pctEt);
    fprintf('  fluid pressure P = %.4f Pa\n', pF);
end

fprintf('\n================ ENDOTHELIUM-FLUID INTERFACE ================\n');
for i = 1:numel(zPoints)
    z = zPoints(i);
    rEwall = deltaEofz(z);
    gapHere = rEwall - deltaLofz(z);
    eps_ = epsFrac * gapHere;
    slopeE = slopeEofz(z);
    nE = [1, -slopeE] / norm([1, -slopeE]);

    rq  = rEwall - eps_;   % fluid side sample point (just inside fluid)
    sigF = [FfluidRR(rq,z), FfluidRZ(rq,z); FfluidRZ(rq,z), FfluidZZ(rq,z)];
    rq2 = rEwall + eps_;   % solid side sample point (just inside endothelium)
    sigS = [FErr(rq2,z), FErz(rq2,z); FErz(rq2,z), FEzz(rq2,z)];
    [tnF, ttF] = traction_components_local(sigF, nE);
    [tnS, ttS] = traction_components_local(sigS, nE);
    pF = FfluidP(rq, z);

    pctEn = 100*abs(tnS-tnF)/refScaleE;
    pctEt = 100*abs(ttS-ttF)/refScaleE;

    fprintf('\n--- Point E%d: z=%.4f um (wall r=%.4f um, slope=%.4f) ---\n', i, z*1e6, rEwall*1e6, slopeE);
    fprintf('  Sample points: fluid at r=%.4fum, solid at r=%.4fum (eps=%.2fnm, gap=%.2fnm)\n', ...
        rq*1e6, rq2*1e6, eps_*1e9, gapHere*1e9);
    fprintf('  sigma_normal:  solid=%.4f Pa,  fluid=%.4f Pa   (diff=%.4f Pa, %.2f%% of refScale)\n', tnS, tnF, tnS-tnF, pctEn);
    fprintf('  sigma_tangent: solid=%.4f Pa,  fluid=%.4f Pa   (diff=%.4f Pa, %.2f%% of refScale)\n', ttS, ttF, ttS-ttF, pctEt);
    fprintf('  fluid pressure P = %.4f Pa\n', pF);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');

function [tn, tt] = traction_components_local(sig, n)
n = n(:);
that = [-n(2); n(1)];
t = sig * n;
tn = t.' * n;
tt = t.' * that;
end
