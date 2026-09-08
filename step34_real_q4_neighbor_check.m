%% STEP34_REAL_Q4_NEIGHBOR_CHECK
% Corrected version of Steps 3-4. Previous version used scatteredInterpolant
% (Delaunay triangulation, always 3 vertices in 2D) to approximate "the
% points used for interpolation" -- but that is NOT what
% recover_fluid_nodes_pressure_stress_Q4.m actually does. That function
% computes stress at each fluid element from that element's own 4 CORNER

% diagram (1 solid point, 4 neighbor arrows), not a 3-vertex triangle.
%
% This script finds the REAL enclosing Q4 element for each of the 4
% correction sample points -- via structured-grid bracketing on the
% actual meshF.Rp/meshF.Zp body-fitted grid, the same grid
% add_fluid_nodes.m turns into meshF.elems -- and checks whether each of
% the 4 real corner nodes is genuinely inside the fluid phase at its own z.

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

meshF = add_fluid_nodes(fl.meshF);   % adds .nodes (Rp,Zp) and .elems (4-node connectivity)
mraw = fl.meshF;

zGrid = out.z(:);
deltaE_k = out.deltaEHist(:,k);
deltaL_k = out.deltaLHist(:,k);
deltaLofz = @(z) interp1(zGrid, deltaL_k, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(zGrid, deltaE_k, z, 'linear', 'extrap');

epsFrac = 0.1;
zPoints = [3.45e-6, 3.50e-6];

function report_point(label, rq, zq, mraw, meshF, deltaLofz, deltaEofz)
    Nr = mraw.Nr; Nz = mraw.Nz;
    Zcol = mraw.Zp(1,:);              % z depends only on column j (structured)
    [~, j0] = min(abs(Zcol - zq));
    j = min(max(j0,1), Nz-1);         % left column index of bracketing pair
    if Zcol(j) > zq && j > 1, j = j-1; end

    Rcol = mraw.Rp(:,j);              % radial profile at column j (monotonic in i)
    i0 = find(Rcol <= rq, 1, 'last');
    if isempty(i0), i0 = 1; end
    i = min(max(i0,1), Nr-1);

    % element index matching add_fluid_nodes.m's e = (j-1)*(Nr-1) + i ordering
    e = (j-1)*(Nr-1) + i;
    conn = meshF.elems(e,:);
    corners = meshF.nodes(conn,:);    % 4 corner (r,z) pairs, real Q4 element

    fprintf('\n--- %s: query (r=%.4f, z=%.4f) um -> element (i=%d,j=%d) ---\n', ...
        label, rq*1e6, zq*1e6, i, j);
    anyAlert = false;
    for c = 1:4
        rc = corners(c,1); zc = corners(c,2);
        rWallL = deltaLofz(zc); rWallE = deltaEofz(zc);
        inGap = rc > rWallL && rc < rWallE;
        tag = 'OK (inside fluid gap)';
        if ~inGap
            tag = 'ALERT: outside supposed phase!';
            anyAlert = true;
        end
        fprintf('  corner %d: (r=%.4f, z=%.4f) um -- rWallL=%.4f, rWallE=%.4f at that z -- %s\n', ...
            c, rc*1e6, zc*1e6, rWallL*1e6, rWallE*1e6, tag);
    end
    if ~anyAlert
        fprintf('  --> All 4 real Q4 corner nodes confirmed inside the true fluid gap. No alert.\n');
    else
        fprintf('  --> ALERT triggered for %s.\n', label);
    end
end

fprintf('================ LEUKOCYTE-FLUID INTERFACE (real Q4 element, 4 corners) ================\n');
labelsL = {'L1','L2'};
for idx = 1:numel(zPoints)
    z = zPoints(idx);
    rWallL_i = deltaLofz(z); rWallE_i = deltaEofz(z);
    gapHere = rWallE_i - rWallL_i;
    rq = rWallL_i + epsFrac*gapHere;
    report_point(labelsL{idx}, rq, z, mraw, meshF, deltaLofz, deltaEofz);
end

fprintf('\n================ ENDOTHELIUM-FLUID INTERFACE (real Q4 element, 4 corners) ================\n');
labelsE = {'E1','E2'};
for idx = 1:numel(zPoints)
    z = zPoints(idx);
    rWallL_i = deltaLofz(z); rWallE_i = deltaEofz(z);
    gapHere = rWallE_i - rWallL_i;
    rq = rWallE_i - epsFrac*gapHere;
    report_point(labelsE{idx}, rq, z, mraw, meshF, deltaLofz, deltaEofz);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
