%% STEP1_FIND_THIN_REGION

% adds: pick it AWAY from where the endothelium boundary is problematic
% (z > ~3.9-4.0 um, close to the mesh's truncated edge) to remove that
% confound. Checks the FULL gap(z) profile, not just the global minimum,
% to see whether a second, independent thin region exists elsewhere.

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

zGrid = out.z(:);
k = out.stopStep;
deltaE_k = out.deltaEHist(:,k);
deltaL_k = out.deltaLHist(:,k);
gap = deltaE_k - deltaL_k;

meshEmaxZ = max(out.meshE.nodes(:,2)); % actual, undeformed reference extent -- deformed extent close to this
fprintf('Endothelium mesh z-range (reference): [%.4f, %.4f] um\n', min(out.meshE.nodes(:,2))*1e6, meshEmaxZ*1e6);

fprintf('\nFull gap(z) profile at step %d:\n', k);
fprintf('%8s %10s %14s\n', 'z[um]', 'gap[nm]', 'note');
for i = 1:numel(zGrid)
    note = '';
    if zGrid(i) > meshEmaxZ - 0.3e-6
        note = '<-- near/past problematic endothelium edge';
    end
    fprintf('%8.3f %10.2f  %s\n', zGrid(i)*1e6, gap(i)*1e9, note);
end

% Find local minima away from the problematic edge region
safeIdx = zGrid < (meshEmaxZ - 0.5e-6) & zGrid > min(out.meshE.nodes(:,2)) + 0.5e-6;
[gapSorted, order] = sort(gap(safeIdx));
zSafe = zGrid(safeIdx);
fprintf('\n=== Thinnest gap locations AWAY from the problematic edge (z < %.3fum, z > %.3fum) ===\n', ...
    (meshEmaxZ-0.5e-6)*1e6, (min(out.meshE.nodes(:,2))+0.5e-6)*1e6);
for i = 1:min(10, numel(gapSorted))
    fprintf('  z=%.4f um, gap=%.2f nm\n', zSafe(order(i))*1e6, gapSorted(i)*1e9);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
