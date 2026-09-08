%% REGENERATE_SLIDE4_UNIFIED
% The actual "regenerated slide 4" deliverable: L1, L2, E1, E2 (same
% z=3.45/3.50 points from Step 1/2) evaluated on the NEW unified-domain
% (wide, now-default) 40-step run, using the now-patched (position-aware)
% compute_interface_traction_mismatch.m -- both fixes active.

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/code/leukocyte-main/out_2D_unified_domain_40step.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

k = out.stopStep;
st = out.stateHist{k};
fl = out.fluidHist{k};
par = out.par;
dtStep = out.dtHist(k);
fprintf('Unified-domain run: stopStep=%d\n', k);

parLmismatch = par;
if isfield(par, 'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par, 'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par, 'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

mismatchOpts = struct('nQuery', 61, 'epsFrac', 0.1, 'trimFrac', 0.05);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);
refScaleL = leuko.stats.refScale;
refScaleE = endo.stats.refScale;
fprintf('refScale: leuko=%.4f endo=%.4f Pa\n', refScaleL, refScaleE);

fprintf('\nOverall worst-case (this final state, whole interface):\n');
fprintf('  Leuko: maxEn=%.2f%% maxEt=%.2f%%\n', leuko.stats.pctMaxEn, leuko.stats.pctMaxEt);
fprintf('  Endo:  maxEn=%.2f%% maxEt=%.2f%%\n', endo.stats.pctMaxEn, endo.stats.pctMaxEt);

% Now the specific L1/L2/E1/E2 points -- same z=3.45/3.50 as Step 1/2 --
% evaluated within this NEW unified-domain geometry.
mesh_zc = fl.meshF;
if ~all(isfield(mesh_zc, {'zc','deltaL_c','deltaE_c'}))
    error('missing zc/deltaL_c/deltaE_c');
end
deltaLofz = @(z) interp1(mesh_zc.zc, mesh_zc.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh_zc.zc, mesh_zc.deltaE_c, z, 'linear', 'extrap');
dz_fd = 1e-8;
slopeLofz = @(z) (deltaLofz(z+dz_fd) - deltaLofz(z-dz_fd)) / (2*dz_fd);
slopeEofz = @(z) (deltaEofz(z+dz_fd) - deltaEofz(z-dz_fd)) / (2*dz_fd);
epsFrac = 0.1;

meshF2 = add_fluid_nodes(fl.meshF);
mraw = fl.meshF;

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

function [i,j] = ig_ij(mraw, rq, z)
Nr = mraw.Nr; Nz = mraw.Nz;
Zcol = mraw.Zp(1,:);
[~,j0] = min(abs(Zcol - z)); j = min(max(j0,1),Nz-1);
if Zcol(j) > z && j>1, j=j-1; end
Rcol = mraw.Rp(:,j);
i0 = find(Rcol <= rq, 1, 'last'); if isempty(i0), i0=1; end
i = min(max(i0,1),Nr-1);
end

function [tn, tt] = proj(sig, n)
n = n(:); that = [-n(2); n(1)];
t = sig * n; tn = t.'*n; tt = t.'*that;
end

zPoints = [3.45e-6, 3.50e-6];
labelsL = {'L1','L2'}; labelsE = {'E1','E2'};

fprintf('\n================ REGENERATED SLIDE 4: L1, L2, E1, E2 on unified domain ================\n');
for idx = 1:2
    z = zPoints(idx);
    rWallL = deltaLofz(z); gap = deltaEofz(z)-rWallL; eps_ = epsFrac*gap;
    slope = slopeLofz(z); n = [1,-slope]/norm([1,-slope]);
    rqF = rWallL+eps_; rqS = rWallL-eps_;
    [ig,jg] = ig_ij(mraw, rqF, z);
    sigFvec = locate_and_interp_fluid_stress(mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, rqF, z, ig, jg);
    sigF = [sigFvec(1), sigFvec(4); sigFvec(4), sigFvec(3)];
    sigS = [FLrr(rqS,z), FLrz(rqS,z); FLrz(rqS,z), FLzz(rqS,z)];
    [tnF,ttF] = proj(sigF,n); [tnS,ttS] = proj(sigS,n);
    pctEn = 100*abs(tnS-tnF)/refScaleL; pctEt = 100*abs(ttS-ttF)/refScaleL;
    fprintf('%s: gap=%.1fnm | sig(S)=%.2f sig(F)=%.2f | tau(S)=%.2f tau(F)=%.2f | mismatch n=%.1f%% t=%.1f%%\n', ...
        labelsL{idx}, gap*1e9, tnS, tnF, ttS, ttF, pctEn, pctEt);
end
for idx = 1:2
    z = zPoints(idx);
    rWallE = deltaEofz(z); gap = rWallE-deltaLofz(z); eps_ = epsFrac*gap;
    slope = slopeEofz(z); n = [1,-slope]/norm([1,-slope]);
    rqF = rWallE-eps_; rqS = rWallE+eps_;
    [ig,jg] = ig_ij(mraw, rqF, z);
    sigFvec = locate_and_interp_fluid_stress(mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, rqF, z, ig, jg);
    sigF = [sigFvec(1), sigFvec(4); sigFvec(4), sigFvec(3)];
    sigS = [FErr(rqS,z), FErz(rqS,z); FErz(rqS,z), FEzz(rqS,z)];
    [tnF,ttF] = proj(sigF,n); [tnS,ttS] = proj(sigS,n);
    pctEn = 100*abs(tnS-tnF)/refScaleE; pctEt = 100*abs(ttS-ttF)/refScaleE;
    fprintf('%s: gap=%.1fnm | sig(S)=%.2f sig(F)=%.2f | tau(S)=%.2f tau(F)=%.2f | mismatch n=%.1f%% t=%.1f%%\n', ...
        labelsE{idx}, gap*1e9, tnS, tnF, ttS, ttF, pctEn, pctEt);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
