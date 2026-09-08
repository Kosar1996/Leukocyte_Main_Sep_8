%% CHECK_TRACTION_DEFINITIONS_AGREE
% Tests the hypothesis that the traction the correction is DRIVEN against
% (compute_bodyfitted_wall_traction.m: raw pressure, unrotated) differs
% from the traction the %mismatch check JUDGES it against
% (compute_interface_traction_mismatch.m: full fluid stress tensor,
% rotated into the true local normal/tangent direction accounting for
% wall slope). If these disagree, especially at the sloped contact zone,
% that's the actual reason the correction can converge to its own
% equilibrium while still reading as a large mismatch.

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

zGrid    = out.z(:);
deltaE_k = out.deltaEHist(:,k);
deltaL_k = out.deltaLHist(:,k);
gap      = deltaE_k - deltaL_k;
[minGap, iMin] = min(gap);
zContact = zGrid(iMin);

%% Method A: what the CORRECTION is actually driven against
[~, trE_A] = compute_bodyfitted_wall_traction(fl.meshF, fl, par);
normalA = interp1(trE_A.z, trE_A.normal, zGrid, 'linear', 'extrap');
tangentA = interp1(trE_A.z, trE_A.tangent, zGrid, 'linear', 'extrap');

%% Method B: what the %mismatch CHECK actually judges it against
meshF2 = add_fluid_nodes(fl.meshF);
[~, sigmaCellF, centerF] = recover_fluid_nodes_pressure_stress_Q4(meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell);
FfluidRR = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,1), 'linear', 'nearest');
FfluidZZ = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,3), 'linear', 'nearest');
FfluidRZ = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,4), 'linear', 'nearest');

dz_fd = 1e-8;
deltaEofz = @(z) interp1(zGrid, deltaE_k, z, 'linear', 'extrap');
slopeEofz = @(z) (deltaEofz(z+dz_fd) - deltaEofz(z-dz_fd)) / (2*dz_fd);

epsFrac = 0.1;
zSample = linspace(min(zGrid)+0.3e-6, max(zGrid)-0.3e-6, 20)';
normalB = nan(size(zSample));
tangentB = nan(size(zSample));
slopeAtZ = nan(size(zSample));

for i = 1:numel(zSample)
    z = zSample(i);
    rEwall = deltaEofz(z);
    rLwall = interp1(zGrid, deltaL_k, z, 'linear', 'extrap');
    gapHere = rEwall - rLwall;
    eps_ = epsFrac * gapHere;

    slopeE = slopeEofz(z);
    slopeAtZ(i) = slopeE;
    nE = [1, -slopeE] / norm([1, -slopeE]);

    rq = rEwall - eps_; % just inside the fluid
    sigF = [FfluidRR(rq,z), FfluidRZ(rq,z); FfluidRZ(rq,z), FfluidZZ(rq,z)];
    that = [-nE(2); nE(1)];
    t = sigF * nE(:);
    normalB(i) = t.' * nE(:);
    tangentB(i) = t.' * that;
end

normalA_atSample = interp1(zGrid, normalA, zSample, 'linear', 'extrap');
tangentA_atSample = interp1(zGrid, tangentA, zSample, 'linear', 'extrap');

fprintf('Contact zone (min gap) at z=%.3f um\n\n', zContact*1e6);
fprintf('%8s %8s | %10s %10s | %10s %10s | %9s %9s\n', ...
    'z[um]', 'slope', 'normalA', 'normalB', 'tangentA', 'tangentB', 'diffN%', 'diffT%');
for i = 1:numel(zSample)
    diffN = 100*abs(normalA_atSample(i)-normalB(i))/max(abs(normalB(i)),1e-9);
    diffT = 100*abs(tangentA_atSample(i)-tangentB(i))/max(abs(tangentB(i)),1e-9);
    marker = '';
    if abs(zSample(i)-zContact) < 0.3e-6
        marker = '  <-- near contact zone';
    end
    fprintf('%8.3f %8.4f | %10.4f %10.4f | %10.4f %10.4f | %8.1f%% %8.1f%%%s\n', ...
        zSample(i)*1e6, slopeAtZ(i), normalA_atSample(i), normalB(i), ...
        tangentA_atSample(i), tangentB(i), diffN, diffT, marker);
end
fprintf('\nMethod A = compute_bodyfitted_wall_traction.m (what drives the correction)\n');
fprintf('Method B = compute_interface_traction_mismatch.m style, rotated for slope (what the %%mismatch check uses)\n');
fprintf('Large disagreement, especially where |slope| is large, would confirm the correction is driven\n');
fprintf('against a simpler target than what it is actually judged against.\n');
fprintf('\nSMOKE_TEST_STATUS: OK\n');
