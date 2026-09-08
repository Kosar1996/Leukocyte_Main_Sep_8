%% STEP3_REDO_HEATMAP_MATRIX
% Corrected Step 3: build an actual regular "heatmap matrix" (r,z grid,
% each point classified and evaluated as leukocyte/endothelium/fluid
% depending on position, matching the style of the original annotated
% plot and my earlier zoomed-contour script), then find the closest
% HEATMAP grid point to each of the 4 step-2 interface points, and check
% whether that point's domain classification makes physical sense given


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
FrrE = scatteredInterpolant(rE, zE, stressE.sigma_rr, 'linear', 'nearest');
FzzE = scatteredInterpolant(rE, zE, stressE.sigma_zz, 'linear', 'nearest');
FrzE = scatteredInterpolant(rE, zE, stressE.sigma_rz, 'linear', 'nearest');
FrrL = scatteredInterpolant(rL, zL, stressL.sigma_rr, 'linear', 'nearest');
FzzL = scatteredInterpolant(rL, zL, stressL.sigma_zz, 'linear', 'nearest');
FrzL = scatteredInterpolant(rL, zL, stressL.sigma_rz, 'linear', 'nearest');

meshF2 = add_fluid_nodes(fl.meshF);
[pCellF, sigmaCellF, centerF] = recover_fluid_nodes_pressure_stress_Q4(meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell);
FrrFl = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,1), 'linear', 'nearest');
FzzFl = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,3), 'linear', 'nearest');
FrzFl = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,4), 'linear', 'nearest');
FpFl  = scatteredInterpolant(centerF(:,1), centerF(:,2), pCellF,          'linear', 'nearest');

zGrid = out.z(:);
deltaE_k = out.deltaEHist(:,k);
deltaL_k = out.deltaLHist(:,k);
deltaLofz = @(z) interp1(zGrid, deltaL_k, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(zGrid, deltaE_k, z, 'linear', 'extrap');

%% Build the regular heatmap matrix, zoomed to the thin-fluid region
zHM = linspace(3.30e-6, 3.65e-6, 351)';   % ~1nm spacing in z
rHM = linspace(3.00e-6, 3.30e-6, 301)';   % ~1nm spacing in r
[RRm, ZZm] = meshgrid(rHM, zHM);
nZ = numel(zHM); nR = numel(rHM);

VALrr = nan(nZ, nR); VALzz = nan(nZ, nR); VALrz = nan(nZ, nR); VALp = nan(nZ,nR);
CLASS = strings(nZ, nR); % "leuko", "endo", "fluid"

for i = 1:nZ
    z = zHM(i);
    rWallL = deltaLofz(z);
    rWallE = deltaEofz(z);
    for j = 1:nR
        r = rHM(j);
        if r <= rWallL
            if z >= min(zL) && z <= max(zL)
                CLASS(i,j) = "leuko";
                VALrr(i,j) = FrrL(r,z); VALzz(i,j) = FzzL(r,z); VALrz(i,j) = FrzL(r,z);
            end
        elseif r >= rWallE
            if z >= min(zE) && z <= max(zE)
                CLASS(i,j) = "endo";
                VALrr(i,j) = FrrE(r,z); VALzz(i,j) = FzzE(r,z); VALrz(i,j) = FrzE(r,z);
            end
        else
            CLASS(i,j) = "fluid";
            VALrr(i,j) = FrrFl(r,z); VALzz(i,j) = FzzFl(r,z); VALrz(i,j) = FrzFl(r,z); VALp(i,j) = FpFl(r,z);
        end
    end
end

%% For each of the 4 interface points, find the closest HEATMAP grid coordinate
zPoints = [3.45e-6, 3.50e-6];
labels = {'L1','L2'; 'E1','E2'};

epsFrac = 0.1; % same default as compute_interface_traction_mismatch.m / step 2

fprintf('================ LEUKOCYTE-FLUID INTERFACE (heatmap matrix) ================\n');
for i = 1:numel(zPoints)
    z = zPoints(i);
    rWallL_i = deltaLofz(z);
    rWallE_i = deltaEofz(z);
    gapHere = rWallE_i - rWallL_i;
    rWall = rWallL_i + epsFrac*gapHere; % the ACTUAL fluid-side sample point step 2 used, not the bare wall coordinate
    D2 = (ZZm-z).^2 + (RRm-rWall).^2;
    D2(CLASS=="") = inf; % skip undefined cells
    [minD2, lin] = min(D2(:));
    [ii,jj] = ind2sub(size(D2), lin);
    fprintf('\n--- %s: fluid-side sample point used by the correction (r=%.4f, z=%.4f) um, wall at r=%.4f ---\n', ...
        labels{1,i}, rWall*1e6, z*1e6, rWallL_i*1e6);
    fprintf('  Closest heatmap coordinate: (r=%.4f, z=%.4f) um, distance=%.2f nm, classified as: %s\n', ...
        RRm(ii,jj)*1e6, ZZm(ii,jj)*1e6, sqrt(minD2)*1e9, CLASS(ii,jj));
    if CLASS(ii,jj) ~= "fluid"
        fprintf('  --> ALERT: closest heatmap coordinate to the FLUID interface is classified as %s, not fluid.\n', CLASS(ii,jj));
    end
    fprintf('  sigma_rr=%.4f, sigma_zz=%.4f, sigma_rz=%.4f Pa', VALrr(ii,jj), VALzz(ii,jj), VALrz(ii,jj));
    if CLASS(ii,jj)=="fluid"
        fprintf(', P=%.4f Pa', VALp(ii,jj));
    end
    fprintf('\n');
end

fprintf('\n================ ENDOTHELIUM-FLUID INTERFACE (heatmap matrix) ================\n');
for i = 1:numel(zPoints)
    z = zPoints(i);
    rWallL_i = deltaLofz(z);
    rWallE_i = deltaEofz(z);
    gapHere = rWallE_i - rWallL_i;
    rWall = rWallE_i - epsFrac*gapHere; % fluid-side sample point, offset INTO the fluid from the endothelium wall
    D2 = (ZZm-z).^2 + (RRm-rWall).^2;
    D2(CLASS=="") = inf;
    [minD2, lin] = min(D2(:));
    [ii,jj] = ind2sub(size(D2), lin);
    fprintf('\n--- %s: fluid-side sample point used by the correction (r=%.4f, z=%.4f) um, wall at r=%.4f ---\n', ...
        labels{2,i}, rWall*1e6, z*1e6, rWallE_i*1e6);
    fprintf('  Closest heatmap coordinate: (r=%.4f, z=%.4f) um, distance=%.2f nm, classified as: %s\n', ...
        RRm(ii,jj)*1e6, ZZm(ii,jj)*1e6, sqrt(minD2)*1e9, CLASS(ii,jj));
    if CLASS(ii,jj) ~= "fluid"
        fprintf('  --> ALERT: closest heatmap coordinate to the FLUID interface is classified as %s, not fluid.\n', CLASS(ii,jj));
    end
    fprintf('  sigma_rr=%.4f, sigma_zz=%.4f, sigma_rz=%.4f Pa', VALrr(ii,jj), VALzz(ii,jj), VALrz(ii,jj));
    if CLASS(ii,jj)=="fluid"
        fprintf(', P=%.4f Pa', VALp(ii,jj));
    end
    fprintf('\n');
end

%% Also check the boundary itself: scan along z at the interface r, see if
% classification flips inconsistently right at the wall (a direct test of
% "on the heatmap it shows as a fluid point but its coordinate suggests
% it should be a solid point")
fprintf('\n================ Classification consistency scan near the leukocyte wall ================\n');
fprintf('%8s %10s %10s %10s\n', 'z[um]', 'rWallL', 'nearestR', 'class');
for i = 1:20:nZ
    z = zHM(i);
    rWallL = deltaLofz(z);
    [~, jNear] = min(abs(rHM - rWallL));
    fprintf('%8.4f %10.4f %10.4f %10s\n', z*1e6, rWallL*1e6, rHM(jNear)*1e6, CLASS(i,jNear));
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
