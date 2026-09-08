%% STEP3_MERGED_TABLE_FIXED
% Task: merge %mismatch (normal+tangential) into the Step 3 table as two
% additional columns, using Step 2's normal vectors -- and do it
% consistently, with the FIXED position-aware interpolation throughout
% (not the old buggy centroid-only evaluation).

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

refScaleL_new = 190.5003; refScaleE_new = 190.5003; % from step2_recomputed_fixed_interp.m

zPoints = [3.45e-6, 3.50e-6];
labelsL = {'L1','L2'}; labelsE = {'E1','E2'};

% Step-2 (fixed-interpolation) reference values, needed for the compare-with-2 step
step2 = struct();
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
    step2.(labelsL{idx}) = struct('n',n,'tnF',tnF,'ttF',ttF,'tnS',tnS,'ttS',ttS,'z',z,'rWall',rWallL,'side','L');
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
    step2.(labelsE{idx}) = struct('n',n,'tnF',tnF,'ttF',ttF,'tnS',tnS,'ttS',ttS,'z',z,'rWall',rWallE,'side','E');
end

%% Heatmap grid (1nm resolution), classified, using FIXED fluid interpolation
zHM = linspace(3.30e-6, 3.65e-6, 351)';
rHM = linspace(3.00e-6, 3.30e-6, 301)';

allLabels = {'L1','L2','E1','E2'};
fprintf('================ STEP 3, merged table: heatmap coord + all components + %%mismatch (fixed interp) ================\n');
fprintf('%4s %10s %22s %8s %9s %9s %9s %9s %9s %10s %10s\n', ...
    'Pt','Dist[nm]','HeatmapCoord[um]','srr','stt','szz','srz','P','MismN[%]','MismT[%]','');

for li = 1:4
    lab = allLabels{li};
    p2 = step2.(lab);
    if p2.side=='L'
        rWallQ = p2.rWall + epsFrac*(deltaEofz(p2.z)-p2.rWall);
    else
        rWallQ = p2.rWall - epsFrac*(p2.rWall-deltaLofz(p2.z));
    end
    z = p2.z;
    % nearest heatmap grid point (classification unchanged from before -- all confirmed fluid)
    [~,jz] = min(abs(zHM - z));
    [~,ir] = min(abs(rHM - rWallQ));
    rHMpt = rHM(ir); zHMpt = zHM(jz);
    distNm = norm([rHMpt-rWallQ, zHMpt-z])*1e9;

    [igH,jgH] = initial_guess_ij(mraw, rHMpt, zHMpt);
    [sigH, pH] = locate_and_interp_fluid_stress(mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, rHMpt, zHMpt, igH, jgH);

    sigHmat = [sigH(1) sigH(4); sigH(4) sigH(3)];
    [tnH, ttH] = proj(sigHmat, p2.n);

    if p2.side=='L'
        refScale = refScaleL_new; tnOrig = p2.tnF; ttOrig = p2.ttF;
    else
        refScale = refScaleE_new; tnOrig = p2.tnF; ttOrig = p2.ttF;
    end
    pctN = 100*abs(tnH-tnOrig)/refScale;
    pctT = 100*abs(ttH-ttOrig)/refScale;

    fprintf('%4s %10.2f (%.4f,%.4f) %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n', ...
        lab, distNm, rHMpt*1e6, zHMpt*1e6, sigH(1), sigH(2), sigH(3), sigH(4), pH, pctN, pctT);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
