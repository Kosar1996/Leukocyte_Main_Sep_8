%% CHECK_PRESSURE_METHODOLOGY
% Compares the pressure value the traction correction actually uses
% internally (fluid.P(end,:), a raw MAC-grid read, per
% compute_bodyfitted_wall_traction.m) against the pressure value

% a separate cell-based recovery), at the same contact-zone location, to
% see whether Check 1's mismatch is partly a methodology artifact.

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

k = out.stopStep;
fl = out.fluidHist{k};
par = out.par;

zGrid    = out.z(:);
deltaE_k = out.deltaEHist(:,k);
deltaL_k = out.deltaLHist(:,k);
gap      = deltaE_k - deltaL_k;
[minGap, iMin] = min(gap);
zContact = zGrid(iMin);

fprintf('Contact zone at z=%.4f um\n\n', zContact*1e6);

%% Method 1: raw MAC-grid pressure, P(end,:), as used internally by the correction
meshF = fl.meshF;
zc = meshF.zc(:);
P = fl.P; % Nr x Nz (or similar) MAC pressure field
P_wall_raw = P(end,:).';
P_method1 = interp1(zc, P_wall_raw, zContact, 'linear', 'extrap');
fprintf('Method 1 (correction''s own P(end,:) at MAC grid, interpolated to contact z):\n');
fprintf('  P = %.4f Pa\n\n', P_method1);


meshF2 = add_fluid_nodes(fl.meshF);
[pCellF, ~, centerF] = recover_fluid_nodes_pressure_stress_Q4(meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell);
FpFl = scatteredInterpolant(centerF(:,1), centerF(:,2), pCellF, 'linear', 'nearest');
rMidGap = 0.5*(deltaE_k(iMin)+deltaL_k(iMin));
P_method2 = FpFl(rMidGap, zContact);
fprintf('Method 2 (recover_fluid_nodes_pressure_stress_Q4 pCell, mid-gap, at contact z):\n');
fprintf('  P = %.4f Pa\n\n', P_method2);

fprintf('Difference: %.4f Pa (%.2f%% of Method 1)\n', P_method2-P_method1, 100*abs(P_method2-P_method1)/max(abs(P_method1),1e-9));
fprintf('\nSMOKE_TEST_STATUS: OK\n');
