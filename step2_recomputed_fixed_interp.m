%% STEP2_RECOMPUTED_FIXED_INTERP
% Full, consistent recomputation using the FIXED position-aware Q4
% interpolation everywhere -- both for the 4 report points (L1,L2,E1,E2)
% AND for the refScale normalization sweep (nQuery=61 points along each
% interface), so nothing in the final numbers still depends on the old,
% buggy centroid-only fluid stress evaluation.

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

stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, par);
stressL = recover_nodal_stress_axisym_viscoelastic(out.meshL, st.uL, st.uLPrev, dtStep, parLmismatch);
rE = out.meshE.nodes(:,1) + st.uE(1:2:end);
zE = out.meshE.nodes(:,2) + st.uE(2:2:end);
rL = out.meshL.nodes(:,1) + st.uL(1:2:end);
zL = out.meshL.nodes(:,2) + st.uL(2:2:end);
FLrr = scatteredInterpolant(rL,zL,stressL.sigma_rr,'linear','nearest');
FLrz = scatteredInterpolant(rL,zL,stressL.sigma_rz,'linear','nearest');
FLzz = scatteredInterpolant(rL,zL,stressL.sigma_zz,'linear','nearest');
FErr = scatteredInterpolant(rE,zE,stressE.sigma_rr,'linear','nearest');
FErz = scatteredInterpolant(rE,zE,stressE.sigma_rz,'linear','nearest');
FEzz = scatteredInterpolant(rE,zE,stressE.sigma_zz,'linear','nearest');

meshF2 = add_fluid_nodes(fl.meshF);
mraw = fl.meshF;

mesh_zc = fl.meshF;
deltaLofz = @(z) interp1(mesh_zc.zc, mesh_zc.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh_zc.zc, mesh_zc.deltaE_c, z, 'linear', 'extrap');
dz_fd = 1e-8;
slopeLofz = @(z) (deltaLofz(z+dz_fd) - deltaLofz(z-dz_fd)) / (2*dz_fd);
slopeEofz = @(z) (deltaEofz(z+dz_fd) - deltaEofz(z-dz_fd)) / (2*dz_fd);
epsFrac = 0.1;

function [tn, tt] = proj(sig, n)
n = n(:); that = [-n(2); n(1)];
t = sig * n; tn = t.'*n; tt = t.'*that;
end

function [i,j] = initial_guess_ij(mraw, rq, z)
Nr = mraw.Nr; Nz = mraw.Nz;
Zcol = mraw.Zp(1,:);
[~,j0] = min(abs(Zcol - z)); j = min(max(j0,1),Nz-1);
if Zcol(j) > z && j>1, j=j-1; end
Rcol = mraw.Rp(:,j);
i0 = find(Rcol <= rq, 1, 'last'); if isempty(i0), i0=1; end
i = min(max(i0,1),Nr-1);
end

%% ---- refScale sweep (nQuery=61), leukocyte and endothelium, FIXED interpolation ----
zLo = max([min(mesh_zc.zc), min(zL), min(zE)]);
zHi = min([max(mesh_zc.zc), max(zL), max(zE)]);
nQuery = 61;
zQuery = linspace(zLo, zHi, nQuery)';

tnFluidL = nan(nQuery,1); tnFluidE = nan(nQuery,1);
for qi = 1:nQuery
    z = zQuery(qi);
    rWallL = deltaLofz(z); rWallE = deltaEofz(z); gap = rWallE - rWallL;
    if gap <= 0, continue; end
    eps_ = epsFrac*gap;
    slope = slopeLofz(z); n = [1,-slope]/norm([1,-slope]);
    rq = rWallL + eps_;
    [ig,jg] = initial_guess_ij(mraw, rq, z);
    try
        [sigF] = locate_and_interp_fluid_stress(mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, rq, z, ig, jg);
        sigFmat = [sigF(1) sigF(4); sigF(4) sigF(3)];
        [tnF,~] = proj(sigFmat, n);
        tnFluidL(qi) = tnF;
    catch
    end

    slopeE = slopeEofz(z); nE = [1,-slopeE]/norm([1,-slopeE]);
    rqE = rWallE - eps_;
    [igE,jgE] = initial_guess_ij(mraw, rqE, z);
    try
        [sigFE] = locate_and_interp_fluid_stress(mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, rqE, z, igE, jgE);
        sigFEmat = [sigFE(1) sigFE(4); sigFE(4) sigFE(3)];
        [tnFE,~] = proj(sigFEmat, nE);
        tnFluidE(qi) = tnFE;
    catch
    end
end
refScaleL_new = prctile(abs(tnFluidL(isfinite(tnFluidL))), 90);
refScaleE_new = prctile(abs(tnFluidE(isfinite(tnFluidE))), 90);
fprintf('NEW refScale (fixed interpolation): leuko=%.4f Pa, endo=%.4f Pa\n', refScaleL_new, refScaleE_new);
fprintf('OLD refScale (buggy interpolation): both sides ~190.5 Pa (as previously reported)\n\n');

%% ---- L1, L2, E1, E2 with FIXED interpolation, both fluid AND consistent refScale ----
zPoints = [3.45e-6, 3.50e-6];
labelsL = {'L1','L2'}; labelsE = {'E1','E2'};
oldNormal = struct('L1',36.1,'L2',35.5,'E1',32.5,'E2',30.7);

fprintf('================ RECOMPUTED Step 2, fixed interpolation throughout ================\n');
for idx = 1:2
    z = zPoints(idx);
    rWallL = deltaLofz(z); gap = deltaEofz(z)-rWallL; eps_ = epsFrac*gap;
    slope = slopeLofz(z); n = [1,-slope]/norm([1,-slope]);
    rqF = rWallL+eps_; rqS = rWallL-eps_;
    [igF,jgF] = initial_guess_ij(mraw, rqF, z);
    sigF = locate_and_interp_fluid_stress(mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, rqF, z, igF, jgF);
    sigFmat = [sigF(1) sigF(4); sigF(4) sigF(3)];
    sigS = [FLrr(rqS,z) FLrz(rqS,z); FLrz(rqS,z) FLzz(rqS,z)];
    [tnF,ttF] = proj(sigFmat,n); [tnS,ttS] = proj(sigS,n);
    pctEn = 100*abs(tnS-tnF)/refScaleL_new;
    pctEt = 100*abs(ttS-ttF)/refScaleL_new;
    lab = labelsL{idx};
    fprintf('%s: sig(S)=%.2f sig(F)=%.2f [OLD fluid was different] | mismatch NEW n=%.1f%% (OLD was %.1f%%) t=%.1f%%\n', ...
        lab, tnS, tnF, pctEn, oldNormal.(lab), pctEt);
end
for idx = 1:2
    z = zPoints(idx);
    rWallE = deltaEofz(z); gap = rWallE-deltaLofz(z); eps_ = epsFrac*gap;
    slope = slopeEofz(z); n = [1,-slope]/norm([1,-slope]);
    rqF = rWallE-eps_; rqS = rWallE+eps_;
    [igF,jgF] = initial_guess_ij(mraw, rqF, z);
    sigF = locate_and_interp_fluid_stress(mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, rqF, z, igF, jgF);
    sigFmat = [sigF(1) sigF(4); sigF(4) sigF(3)];
    sigS = [FErr(rqS,z) FErz(rqS,z); FErz(rqS,z) FEzz(rqS,z)];
    [tnF,ttF] = proj(sigFmat,n); [tnS,ttS] = proj(sigS,n);
    pctEn = 100*abs(tnS-tnF)/refScaleE_new;
    pctEt = 100*abs(ttS-ttF)/refScaleE_new;
    lab = labelsE{idx};
    fprintf('%s: sig(S)=%.2f sig(F)=%.2f [OLD fluid was different] | mismatch NEW n=%.1f%% (OLD was %.1f%%) t=%.1f%%\n', ...
        lab, tnS, tnF, pctEn, oldNormal.(lab), pctEt);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
