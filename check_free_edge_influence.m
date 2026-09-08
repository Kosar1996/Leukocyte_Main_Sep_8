%% CHECK_FREE_EDGE_INFLUENCE
% Tests whether the unconstrained z-end (z~4.04 um, free/traction-free by
% omission) is measurably pulling sigma_zz toward zero near the contact
% zone in later steps, as the contact point drifts closer to that edge
% (0.44 um away at step 17, down to 0.29 um by step 40).
%
% Method: at the last step, sample sigma_zz along the interface from the
% contact zone (z~3.75) out to the mesh edge (z~4.04), and see whether it
% shows a clear decay-to-zero trend as z approaches the free edge --
% the expected signature of a real free-surface boundary effect.

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

k = out.stopStep;
st = out.stateHist{k};
par = out.par;
dtStep = out.dtHist(k);

stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, par);
rE = out.meshE.nodes(:,1) + st.uE(1:2:end);
zE = out.meshE.nodes(:,2) + st.uE(2:2:end);
FzzE = scatteredInterpolant(rE, zE, stressE.sigma_zz, 'linear', 'nearest');

zGrid = out.z(:);
deltaE_k = out.deltaEHist(:,k);
zEdge = max(zE);

zSample = linspace(3.6e-6, zEdge, 15)';
fprintf('Step %d: mesh z-edge at %.4f um. Sampling sigma_zz from contact zone to edge.\n\n', k, zEdge*1e6);
fprintf('%8s %10s | %10s\n', 'z[um]', 'dist2edge', 'sigma_zz');
for i = 1:numel(zSample)
    zz = zSample(i);
    rWall = interp1(zGrid, deltaE_k, zz, 'linear', 'extrap');
    szz = FzzE(rWall - 0.02e-6, zz);
    fprintf('%8.4f %10.4f | %10.4f\n', zz*1e6, (zEdge-zz)*1e6, szz);
end

fprintf('\nIf sigma_zz trends toward 0 specifically as z approaches the edge, that is the\n');
fprintf('signature of a real free-boundary effect reaching into the region we care about.\n');
fprintf('If it stays roughly flat/noisy with no trend toward the edge, the free boundary is\n');
fprintf('probably not the driver and can be documented rather than fixed.\n');
fprintf('\nSMOKE_TEST_STATUS: OK\n');
