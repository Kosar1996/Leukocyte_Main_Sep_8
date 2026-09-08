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
    points(end+1) = struct('label',sprintf('L%d',idx),'rq',rq,'zq',z,'i',i,'j',j); %#ok<SAGROW>
end
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
    points(end+1) = struct('label',sprintf('E%d',idx),'rq',rq,'zq',z,'i',i,'j',j); %#ok<SAGROW>
end

[pCellF_old, sigmaCellF_old, centerF_old] = recover_fluid_nodes_pressure_stress_Q4(meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell);

fprintf('================ ROBUST position-aware interpolation, all 4 points ================\n');
for pi = 1:numel(points)
    pt = points(pi);
    e_guess = (pt.j-1)*(mraw.Nr-1) + pt.i;
    oldSigma = sigmaCellF_old(e_guess,:);

    [newSigma, newP, rzActual, xieta, ij] = locate_and_interp_fluid_stress( ...
        mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, pt.rq, pt.zq, pt.i, pt.j);

    fprintf('\n--- %s: query (r=%.4f, z=%.4f) um ---\n', pt.label, pt.rq*1e6, pt.zq*1e6);
    fprintf('  Initial guess element (i=%d,j=%d) -> final element (i=%d,j=%d)\n', pt.i, pt.j, ij(1), ij(2));
    fprintf('  Local (xi,eta) = (%.4f, %.4f)\n', xieta(1), xieta(2));
    fprintf('  Query point match: target=(%.4f,%.4f), solved=(%.4f,%.4f), error=%.2e um\n', ...
        pt.rq*1e6, pt.zq*1e6, rzActual(1)*1e6, rzActual(2)*1e6, norm(rzActual-[pt.rq,pt.zq])*1e6);
    fprintf('  OLD (broken, centroid): rr=%.4f tt=%.4f zz=%.4f rz=%.4f\n', oldSigma(1), oldSigma(2), oldSigma(3), oldSigma(4));
    fprintf('  NEW (position-aware):   rr=%.4f tt=%.4f zz=%.4f rz=%.4f P=%.4f\n', ...
        newSigma(1), newSigma(2), newSigma(3), newSigma(4), newP);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
