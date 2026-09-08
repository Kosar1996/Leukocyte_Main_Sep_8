clc;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
k = out.stopStep;
fl = out.fluidHist{k};
par = out.par;

points = struct( ...
    'L1', [3.1050e-6, 3.4500e-6], ...
    'L2', [3.1220e-6, 3.5000e-6], ...
    'E1', [3.2030e-6, 3.4500e-6], ...
    'E2', [3.2040e-6, 3.5000e-6]);

[pCellF, sigmaCellF, centerF] = recover_fluid_nodes_pressure_stress_Q4(fl.meshF, fl.ur2D, fl.uz2D, par.mu, fl.pCell);
Fp   = scatteredInterpolant(centerF(:,1), centerF(:,2), pCellF, 'linear', 'nearest');
Frr  = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,1), 'linear', 'nearest');
Ftt  = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,2), 'linear', 'nearest');
Fzz  = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,3), 'linear', 'nearest');
Frz  = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,4), 'linear', 'nearest');

meshF2 = add_fluid_nodes(fl.meshF);
mraw = fl.meshF;
function [i,j] = ig_ij(mraw, rq, z)
Nr = mraw.Nr; Nz = mraw.Nz;
Zcol = mraw.Zp(1,:);
[~,j0] = min(abs(Zcol - z)); j = min(max(j0,1),Nz-1);
if Zcol(j) > z && j>1, j=j-1; end
Rcol = mraw.Rp(:,j);
i0 = find(Rcol <= rq, 1, 'last'); if isempty(i0), i0=1; end
i = min(max(i0,1),Nr-1);
end

labels = fieldnames(points);
fprintf('%4s | %10s %10s %10s %10s %10s | %10s %10s %10s %10s %10s\n', ...
    'Pt','OLDrr','OLDtt','OLDzz','OLDrz','OLD_P','NEWrr','NEWtt','NEWzz','NEWrz','NEW_P');
for li = 1:numel(labels)
    lab = labels{li};
    rz = points.(lab);
    rHM = rz(1); zHM = rz(2);
    oldRR = Frr(rHM,zHM); oldTT = Ftt(rHM,zHM); oldZZ = Fzz(rHM,zHM); oldRZ = Frz(rHM,zHM); oldP = Fp(rHM,zHM);
    [ig,jg] = ig_ij(mraw, rHM, zHM);
    [sigVec, newP] = locate_and_interp_fluid_stress(mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, rHM, zHM, ig, jg);
    % sigVec = [rr, tt, zz, rz] per locate_and_interp_fluid_stress's own convention
    fprintf('%4s | %10.4f %10.4f %10.4f %10.4f %10.4f | %10.4f %10.4f %10.4f %10.4f %10.4f\n', ...
        lab, oldRR, oldTT, oldZZ, oldRZ, oldP, sigVec(1), sigVec(2), sigVec(3), sigVec(4), newP);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
