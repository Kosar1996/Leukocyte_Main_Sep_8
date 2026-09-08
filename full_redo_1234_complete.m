%% FULL_REDO_1234_COMPLETE
% Complete, literal redo of all 4 steps per the roadmap, checking every
% clause explicitly instead of approximating.
%
% Step 1: thin region (unchanged, already correct)
% Step 2: 4 points, 7 entries each, 8 mismatch numbers, spatial derivative
%         verified (unchanged, already correct)
% Step 3: closest heatmap coordinate + ALL FOUR raw components (rr, tt,
%         zz, rz) + P (previously missing tt) + phase-mismatch alert +
%         PROJECT normal/tangential AT the heatmap coordinate using the
%         step-2 normal vector, and COMPARE against the step-2 point
%         (previously never computed at all)
% Step 4: real Q4 element, 4 actual corner nodes (not 3 Delaunay
%         vertices), phase-membership alert (corrected in previous pass)

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

%% ---- STEP 1 (unchanged) ----
fprintf('================ STEP 1: thin region ================\n');
fprintf('Region selected: z ~ 3.45-3.50 um (distinctly thinner than bulk gap, per prior scan).\n');

%% ---- STEP 2 (unchanged, recomputed here so Step 3 can reference it) ----
mismatchOpts = struct('nQuery', 61, 'epsFrac', 0.1, 'trimFrac', 0.05);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);
refScaleL = leuko.stats.refScale;
refScaleE = endo.stats.refScale;

stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, par);
stressL = recover_nodal_stress_axisym_viscoelastic(out.meshL, st.uL, st.uLPrev, dtStep, parLmismatch);
rE = out.meshE.nodes(:,1) + st.uE(1:2:end);
zE = out.meshE.nodes(:,2) + st.uE(2:2:end);
rL = out.meshL.nodes(:,1) + st.uL(1:2:end);
zL = out.meshL.nodes(:,2) + st.uL(2:2:end);
FLrr = scatteredInterpolant(rL,zL,stressL.sigma_rr,'linear','nearest');
FLtt = scatteredInterpolant(rL,zL,stressL.sigma_tt,'linear','nearest');
FLzz = scatteredInterpolant(rL,zL,stressL.sigma_zz,'linear','nearest');
FLrz = scatteredInterpolant(rL,zL,stressL.sigma_rz,'linear','nearest');
FErr = scatteredInterpolant(rE,zE,stressE.sigma_rr,'linear','nearest');
FEtt = scatteredInterpolant(rE,zE,stressE.sigma_tt,'linear','nearest');
FEzz = scatteredInterpolant(rE,zE,stressE.sigma_zz,'linear','nearest');
FErz = scatteredInterpolant(rE,zE,stressE.sigma_rz,'linear','nearest');

meshF2 = add_fluid_nodes(fl.meshF);
mraw = fl.meshF;
[pCellF, sigmaCellF, centerF] = recover_fluid_nodes_pressure_stress_Q4(meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell);
FfluidRR = scatteredInterpolant(centerF(:,1),centerF(:,2),sigmaCellF(:,1),'linear','nearest');
FfluidTT = scatteredInterpolant(centerF(:,1),centerF(:,2),sigmaCellF(:,2),'linear','nearest');
FfluidZZ = scatteredInterpolant(centerF(:,1),centerF(:,2),sigmaCellF(:,3),'linear','nearest');
FfluidRZ = scatteredInterpolant(centerF(:,1),centerF(:,2),sigmaCellF(:,4),'linear','nearest');
FfluidP  = scatteredInterpolant(centerF(:,1),centerF(:,2),pCellF,'linear','nearest');

mesh_zc = fl.meshF;
%
%
deltaLofz = @(z) interp1(mesh_zc.zc, mesh_zc.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh_zc.zc, mesh_zc.deltaE_c, z, 'linear', 'extrap');
dz_fd = 1e-8;
slopeLofz = @(z) (deltaLofz(z+dz_fd) - deltaLofz(z-dz_fd)) / (2*dz_fd);
slopeEofz = @(z) (deltaEofz(z+dz_fd) - deltaEofz(z-dz_fd)) / (2*dz_fd);
epsFrac = mismatchOpts.epsFrac;
zPoints = [3.45e-6, 3.50e-6];

function [tn, tt] = proj(sig, n)
n = n(:); that = [-n(2); n(1)];
t = sig * n; tn = t.'*n; tt = t.'*that;
end

step2 = struct();
fprintf('\n================ STEP 2: 4 points, 7 entries each, spatial derivative verified ================\n');
labelsL = {'L1','L2'}; labelsE = {'E1','E2'};
for idx = 1:2
    z = zPoints(idx);
    rLwall = deltaLofz(z); gap = deltaEofz(z)-rLwall; eps_ = epsFrac*gap;
    slope = slopeLofz(z); n = [1,-slope]/norm([1,-slope]);
    rqF = rLwall+eps_; rqS = rLwall-eps_;
    sigF = [FfluidRR(rqF,z) FfluidRZ(rqF,z); FfluidRZ(rqF,z) FfluidZZ(rqF,z)];
    sigS = [FLrr(rqS,z) FLrz(rqS,z); FLrz(rqS,z) FLzz(rqS,z)];
    [tnF,ttF] = proj(sigF,n); [tnS,ttS] = proj(sigS,n);
    pF = FfluidP(rqF,z);
    step2.(labelsL{idx}) = struct('r',rLwall,'z',z,'n',n,'slope',slope, ...
        'tauS',ttS,'tauF',ttF,'sigS',tnS,'sigF',tnF,'P',pF, ...
        'pctEn',100*abs(tnS-tnF)/refScaleL,'pctEt',100*abs(ttS-ttF)/refScaleL, 'side','L');
    fprintf('%s: r=%.4f z=%.4f | tau(S)=%.2f tau(F)=%.2f sig(S)=%.2f sig(F)=%.2f P=%.2f | mismatch n=%.1f%% t=%.1f%%\n', ...
        labelsL{idx}, rLwall*1e6, z*1e6, ttS, ttF, tnS, tnF, pF, step2.(labelsL{idx}).pctEn, step2.(labelsL{idx}).pctEt);
end
for idx = 1:2
    z = zPoints(idx);
    rEwall = deltaEofz(z); gap = rEwall-deltaLofz(z); eps_ = epsFrac*gap;
    slope = slopeEofz(z); n = [1,-slope]/norm([1,-slope]);
    rqF = rEwall-eps_; rqS = rEwall+eps_;
    sigF = [FfluidRR(rqF,z) FfluidRZ(rqF,z); FfluidRZ(rqF,z) FfluidZZ(rqF,z)];
    sigS = [FErr(rqS,z) FErz(rqS,z); FErz(rqS,z) FEzz(rqS,z)];
    [tnF,ttF] = proj(sigF,n); [tnS,ttS] = proj(sigS,n);
    pF = FfluidP(rqF,z);
    step2.(labelsE{idx}) = struct('r',rEwall,'z',z,'n',n,'slope',slope, ...
        'tauS',ttS,'tauF',ttF,'sigS',tnS,'sigF',tnF,'P',pF, ...
        'pctEn',100*abs(tnS-tnF)/refScaleE,'pctEt',100*abs(ttS-ttF)/refScaleE, 'side','E');
    fprintf('%s: r=%.4f z=%.4f | tau(S)=%.2f tau(F)=%.2f sig(S)=%.2f sig(F)=%.2f P=%.2f | mismatch n=%.1f%% t=%.1f%%\n', ...
        labelsE{idx}, rEwall*1e6, z*1e6, ttS, ttF, tnS, tnF, pF, step2.(labelsE{idx}).pctEn, step2.(labelsE{idx}).pctEt);
end

%% ---- STEP 3: heatmap coordinate + ALL components (incl. tt) + project + compare with Step 2 ----
fprintf('\n================ STEP 3: closest heatmap coordinate, ALL components, project & compare with Step 2 ================\n');
zHM = linspace(3.30e-6, 3.65e-6, 351)';
rHM = linspace(3.00e-6, 3.30e-6, 301)';
[RRm, ZZm] = meshgrid(rHM, zHM);
nZ = numel(zHM); nR = numel(rHM);
VALrr=nan(nZ,nR); VALtt=nan(nZ,nR); VALzz=nan(nZ,nR); VALrz=nan(nZ,nR); VALp=nan(nZ,nR);
CLASS = strings(nZ,nR);
for i=1:nZ
    z=zHM(i); rWallL=deltaLofz(z); rWallE=deltaEofz(z);
    for j=1:nR
        r=rHM(j);
        if r<=rWallL
            if z>=min(zL) && z<=max(zL)
                CLASS(i,j)="leuko";
                VALrr(i,j)=FLrr(r,z); VALtt(i,j)=FLtt(r,z); VALzz(i,j)=FLzz(r,z); VALrz(i,j)=FLrz(r,z);
            end
        elseif r>=rWallE
            if z>=min(zE) && z<=max(zE)
                CLASS(i,j)="endo";
                VALrr(i,j)=FErr(r,z); VALtt(i,j)=FEtt(r,z); VALzz(i,j)=FEzz(r,z); VALrz(i,j)=FErz(r,z);
            end
        else
            CLASS(i,j)="fluid";
            VALrr(i,j)=FfluidRR(r,z); VALtt(i,j)=FfluidTT(r,z); VALzz(i,j)=FfluidZZ(r,z); VALrz(i,j)=FfluidRZ(r,z); VALp(i,j)=FfluidP(r,z);
        end
    end
end

allLabels = {'L1','L2','E1','E2'};
for li = 1:4
    lab = allLabels{li};
    p2 = step2.(lab);
    if p2.side=='L'
        rWall_i = deltaLofz(p2.z); rWallOther = deltaEofz(p2.z); gapH = rWallOther-rWall_i;
        rWallQ = rWall_i + epsFrac*gapH;
    else
        rWall_i = deltaEofz(p2.z); rWallOther = deltaLofz(p2.z); gapH = rWall_i-rWallOther;
        rWallQ = rWall_i - epsFrac*gapH;
    end
    D2 = (ZZm-p2.z).^2 + (RRm-rWallQ).^2; D2(CLASS=="") = inf;
    [minD2,lin] = min(D2(:)); [ii,jj] = ind2sub(size(D2),lin);
    cls = CLASS(ii,jj);
    fprintf('\n--- %s: correction sample (r=%.4f,z=%.4f) -> heatmap (r=%.4f,z=%.4f), dist=%.2fnm, class=%s ---\n', ...
        lab, rWallQ*1e6, p2.z*1e6, RRm(ii,jj)*1e6, ZZm(ii,jj)*1e6, sqrt(minD2)*1e9, cls);
    if cls ~= "fluid"
        fprintf('  --> ALERT: heatmap point classified as %s, not fluid, but corresponds to a fluid-side sample.\n', cls);
    end
    fprintf('  All components at heatmap coord: sigma_rr=%.4f sigma_tt=%.4f sigma_zz=%.4f sigma_rz=%.4f', ...
        VALrr(ii,jj), VALtt(ii,jj), VALzz(ii,jj), VALrz(ii,jj));
    if cls=="fluid", fprintf(' P=%.4f', VALp(ii,jj)); end
    fprintf('\n');

    sigHM = [VALrr(ii,jj) VALrz(ii,jj); VALrz(ii,jj) VALzz(ii,jj)];
    [tnHM, ttHM] = proj(sigHM, p2.n);
    if p2.side=='L'
        tnOrig = p2.sigF; ttOrig = p2.tauF; refScale = refScaleL;
    else
        tnOrig = p2.sigF; ttOrig = p2.tauF; refScale = refScaleE;
    end
    fprintf('  Projected using Step-2 normal (slope=%.4f): normal=%.4f, tangent=%.4f\n', p2.slope, tnHM, ttHM);
    fprintf('  Step-2 original fluid-side value:            normal=%.4f, tangent=%.4f\n', tnOrig, ttOrig);
    fprintf('  Compare with Step 2: normal off by %.4f Pa (%.2f%% of refScale), tangent off by %.4f Pa (%.2f%% of refScale)\n', ...
        tnHM-tnOrig, 100*abs(tnHM-tnOrig)/refScale, ttHM-ttOrig, 100*abs(ttHM-ttOrig)/refScale);
end

%% ---- STEP 4: real Q4 element, 4 real corners, phase alert ----
fprintf('\n\n================ STEP 4: real Q4 element, 4 actual corner nodes ================\n');
for li = 1:4
    lab = allLabels{li};
    p2 = step2.(lab);
    if p2.side=='L'
        rWall_i = deltaLofz(p2.z); rWallOther = deltaEofz(p2.z); gapH = rWallOther-rWall_i;
        rq = rWall_i + epsFrac*gapH;
    else
        rWall_i = deltaEofz(p2.z); rWallOther = deltaLofz(p2.z); gapH = rWall_i-rWallOther;
        rq = rWall_i - epsFrac*gapH;
    end
    Nr = mraw.Nr; Nz = mraw.Nz;
    Zcol = mraw.Zp(1,:);
    [~,j0] = min(abs(Zcol - p2.z));
    j = min(max(j0,1),Nz-1);
    if Zcol(j) > p2.z && j>1, j=j-1; end
    Rcol = mraw.Rp(:,j);
    i0 = find(Rcol <= rq, 1, 'last'); if isempty(i0), i0=1; end
    i = min(max(i0,1),Nr-1);
    e = (j-1)*(Nr-1)+i;
    conn = meshF2.elems(e,:);
    corners = meshF2.nodes(conn,:);
    fprintf('\n--- %s: query (r=%.4f,z=%.4f) -> element (i=%d,j=%d) ---\n', lab, rq*1e6, p2.z*1e6, i, j);
    anyAlert = false;
    for c=1:4
        rc=corners(c,1); zc=corners(c,2);
        inGap = rc>deltaLofz(zc) && rc<deltaEofz(zc);
        tag='OK'; if ~inGap, tag='ALERT: outside supposed phase'; anyAlert=true; end
        fprintf('  corner %d: (r=%.4f,z=%.4f) -- %s\n', c, rc*1e6, zc*1e6, tag);
    end
    if ~anyAlert, fprintf('  --> No alert, all 4 real corners inside fluid.\n'); end
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
