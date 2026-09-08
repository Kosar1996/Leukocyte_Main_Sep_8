clc; close all;
file = 'out_2D_unified_domain_40step.mat';
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

fprintf('================ HIGH-PRECISION TABLE (4 decimals) ================\n');
fprintf('%4s %12s %12s %12s %12s %12s %14s\n', 'Pt','tang(solid)','tang(fluid)','norm(solid)','norm(fluid)','P[Pa]','MismN/MismT');

for idx = 1:2
    z = zPoints(idx);
    rWallL = deltaLofz(z); gap = deltaEofz(z)-rWallL; eps_ = epsFrac*gap;
    slope = slopeLofz(z); n = [1,-slope]/norm([1,-slope]);
    rqF = rWallL+eps_; rqS = rWallL-eps_;
    [ig,jg] = ig_ij(mraw, rqF, z);
    [sigFvec, pF] = locate_and_interp_fluid_stress(mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, rqF, z, ig, jg);
    sigF = [sigFvec(1), sigFvec(4); sigFvec(4), sigFvec(3)];
    sigS = [FLrr(rqS,z), FLrz(rqS,z); FLrz(rqS,z), FLzz(rqS,z)];
    [tnF,ttF] = proj(sigF,n); [tnS,ttS] = proj(sigS,n);
    pctN = 100*abs(tnS-tnF)/abs(tnF);
    pctT = 100*abs(ttS-ttF)/abs(ttF);
    fprintf('%4s %12.4f %12.4f %12.4f %12.4f %12.4f   N=%.2f%% T=%.2f%%\n', ...
        labelsL{idx}, ttS, ttF, tnS, tnF, pF, pctN, pctT);
    fprintf('     hand-check with 4-decimal values: |%.4f-%.4f|/|%.4f| = %.4f%%\n', ttS, ttF, ttF, 100*abs(ttS-ttF)/abs(ttF));
end
for idx = 1:2
    z = zPoints(idx);
    rWallE = deltaEofz(z); gap = rWallE-deltaLofz(z); eps_ = epsFrac*gap;
    slope = slopeEofz(z); n = [1,-slope]/norm([1,-slope]);
    rqF = rWallE-eps_; rqS = rWallE+eps_;
    [ig,jg] = ig_ij(mraw, rqF, z);
    [sigFvec, pF] = locate_and_interp_fluid_stress(mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, rqF, z, ig, jg);
    sigF = [sigFvec(1), sigFvec(4); sigFvec(4), sigFvec(3)];
    sigS = [FErr(rqS,z), FErz(rqS,z); FErz(rqS,z), FEzz(rqS,z)];
    [tnF,ttF] = proj(sigF,n); [tnS,ttS] = proj(sigS,n);
    pctN = 100*abs(tnS-tnF)/abs(tnF);
    pctT = 100*abs(ttS-ttF)/abs(ttF);
    fprintf('%4s %12.4f %12.4f %12.4f %12.4f %12.4f   N=%.2f%% T=%.2f%%\n', ...
        labelsE{idx}, ttS, ttF, tnS, tnF, pF, pctN, pctT);
    fprintf('     hand-check with 4-decimal values: |%.4f-%.4f|/|%.4f| = %.4f%%\n', ttS, ttF, ttF, 100*abs(ttS-ttF)/abs(ttF));
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
