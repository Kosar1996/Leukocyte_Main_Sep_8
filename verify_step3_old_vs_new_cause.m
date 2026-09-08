clc;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
k = out.stopStep;
fl = out.fluidHist{k};
par = out.par;

rHM = 3.1050e-6; zHM = 3.4500e-6;

fprintf('=== OLD METHOD (recover_fluid_nodes_pressure_stress_Q4.m, buggy centroid-only) ===\n');
[pCellF, sigmaCellF, centerF] = recover_fluid_nodes_pressure_stress_Q4(fl.meshF, fl.ur2D, fl.uz2D, par.mu, fl.pCell);
Fp  = scatteredInterpolant(centerF(:,1), centerF(:,2), pCellF, 'linear', 'nearest');
Frr = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,1), 'linear', 'nearest');
fprintf('  At heatmap coord (r=%.4f, z=%.4f) um: sigma_rr=%.4f  P=%.4f\n', rHM*1e6, zHM*1e6, Frr(rHM,zHM), Fp(rHM,zHM));

fprintf('\n=== NEW METHOD (locate_and_interp_fluid_stress.m, position-aware fix) ===\n');
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
[ig,jg] = ig_ij(mraw, rHM, zHM);
[sigVec, pNew] = locate_and_interp_fluid_stress(mraw, meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell, rHM, zHM, ig, jg);
fprintf('  At heatmap coord (r=%.4f, z=%.4f) um: sigma_rr=%.4f  P=%.4f\n', rHM*1e6, zHM*1e6, sigVec(1), pNew);

fprintf('\nSMOKE_TEST_STATUS: OK\n');
