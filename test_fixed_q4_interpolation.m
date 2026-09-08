%% TEST_FIXED_Q4_INTERPOLATION
% Verify the position-aware fix against the known-broken centroid-only
% version, for all 4 points (L1,L2,E1,E2), and confirm the returned
% "rz_actual" matches the true query point (not the element centroid).

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

meshF2 = add_fluid_nodes(fl.meshF);
mraw = fl.meshF;

mesh_zc = fl.meshF;
deltaLofz = @(z) interp1(mesh_zc.zc, mesh_zc.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh_zc.zc, mesh_zc.deltaE_c, z, 'linear', 'extrap');
epsFrac = 0.1;
zPoints = [3.45e-6, 3.50e-6];

points = struct('label',{},'rq',{},'zq',{},'i',{},'j',{});
% L1, L2
for idx = 1:2
    z = zPoints(idx);
    rWallL = deltaLofz(z); gap = deltaEofz(z)-rWallL; eps_ = epsFrac*gap;
    rq = rWallL + eps_;
    Nr = mraw.Nr; Nz = mraw.Nz;
    Zcol = mraw.Zp(1,:);
    [~,j0] = min(abs(Zcol - z)); j = min(max(j0,1),Nz-1);
    if Zcol(j) > z && j>1, j=j-1; end
    Rcol = mraw.Rp(:,j);
    i0 = find(Rcol <= rq, 1, 'last'); if isempty(i0), i0=1; end
    i = min(max(i0,1),Nr-1);
    lab = sprintf('L%d', idx);
    points(end+1) = struct('label',lab,'rq',rq,'zq',z,'i',i,'j',j); %#ok<SAGROW>
end
% E1, E2
for idx = 1:2
    z = zPoints(idx);
    rWallE = deltaEofz(z); gap = rWallE-deltaLofz(z); eps_ = epsFrac*gap;
    rq = rWallE - eps_;
    Nr = mraw.Nr; Nz = mraw.Nz;
    Zcol = mraw.Zp(1,:);
    [~,j0] = min(abs(Zcol - z)); j = min(max(j0,1),Nz-1);
    if Zcol(j) > z && j>1, j=j-1; end
    Rcol = mraw.Rp(:,j);
    i0 = find(Rcol <= rq, 1, 'last'); if isempty(i0), i0=1; end
    i = min(max(i0,1),Nr-1);
    lab = sprintf('E%d', idx);
    points(end+1) = struct('label',lab,'rq',rq,'zq',z,'i',i,'j',j); %#ok<SAGROW>
end

fprintf('================ OLD (broken, centroid-only) vs NEW (position-aware) ================\n');
[pCellF_old, sigmaCellF_old, centerF_old] = recover_fluid_nodes_pressure_stress_Q4(meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell);

for pi = 1:numel(points)
    pt = points(pi);
    e = (pt.j-1)*(mraw.Nr-1) + pt.i;
    oldSigma = sigmaCellF_old(e,:);
    oldCenter = centerF_old(e,:);

    [newSigma, newP, rzActual, xieta] = recover_fluid_stress_at_point( ...
        mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, pt.i, pt.j, pt.rq, pt.zq);

    fprintf('\n--- %s: query (r=%.4f, z=%.4f) um, element (i=%d,j=%d) ---\n', pt.label, pt.rq*1e6, pt.zq*1e6, pt.i, pt.j);
    fprintf('  OLD (centroid, always at r=%.4f,z=%.4f): rr=%.4f tt=%.4f zz=%.4f rz=%.4f\n', ...
        oldCenter(1)*1e6, oldCenter(2)*1e6, oldSigma(1), oldSigma(2), oldSigma(3), oldSigma(4));
    fprintf('  NEW (xi=%.3f,eta=%.3f, resolves to r=%.4f,z=%.4f): rr=%.4f tt=%.4f zz=%.4f rz=%.4f\n', ...
        xieta(1), xieta(2), rzActual(1)*1e6, rzActual(2)*1e6, newSigma(1), newSigma(2), newSigma(3), newSigma(4));
    fprintf('  Query point match check: target=(%.4f,%.4f), solved=(%.4f,%.4f), error=%.2e um\n', ...
        pt.rq*1e6, pt.zq*1e6, rzActual(1)*1e6, rzActual(2)*1e6, norm(rzActual-[pt.rq,pt.zq])*1e6);
    fprintf('  Change from old to new: d_rr=%.4f d_tt=%.4f d_zz=%.4f d_rz=%.4f Pa\n', ...
        newSigma(1)-oldSigma(1), newSigma(2)-oldSigma(2), newSigma(3)-oldSigma(3), newSigma(4)-oldSigma(4));
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
