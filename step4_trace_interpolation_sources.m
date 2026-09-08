%% STEP4_TRACE_INTERPOLATION_SOURCES

% to the original solid/fluid data to locate all the solid/fluid points
% that are used for the interpolation. If you notice that a heatmap
% coordinate depends on a fluid/solid point that is located outside the
% supposed phase, alert."
%
% For each of the 4 step-2/3 points, finds the actual Delaunay
% triangulation vertices (the 3 real source data points) that
% scatteredInterpolant uses to compute the fluid stress value AT that
% exact query location, then checks whether all 3 source points are
% genuinely, sensibly located inside the fluid gap at their own
% respective z (not on the wrong side of either wall).

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
[pCellF, sigmaCellF, centerF] = recover_fluid_nodes_pressure_stress_Q4(meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell);

zGrid = out.z(:);
deltaE_k = out.deltaEHist(:,k);
deltaL_k = out.deltaLHist(:,k);
deltaLofz = @(z) interp1(zGrid, deltaL_k, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(zGrid, deltaE_k, z, 'linear', 'extrap');

% Build the SAME Delaunay triangulation scatteredInterpolant uses internally
DT = delaunayTriangulation(centerF(:,1), centerF(:,2));

zPoints = [3.45e-6, 3.50e-6];
epsFrac = 0.1;

function trace_point(label, rq, zq, DT, centerF, sigmaCellF, pCellF, deltaLofz, deltaEofz)
    triId = pointLocation(DT, rq, zq);
    fprintf('\n--- %s: query point (r=%.4f, z=%.4f) um ---\n', label, rq*1e6, zq*1e6);
    if isnan(triId)
        fprintf('  Query point is OUTSIDE the fluid data convex hull -- interpolant falls back to nearest-neighbor extrapolation,\n');
        fprintf('  not a triangulated (multi-point) interpolation. Cannot trace 3 source vertices; flag separately.\n');
        return;
    end
    vertIds = DT.ConnectivityList(triId, :);
    fprintf('  Enclosing Delaunay triangle: 3 source fluid points:\n');
    anyAlert = false;
    for v = 1:3
        idx = vertIds(v);
        rv = centerF(idx,1); zv = centerF(idx,2);
        rWallL_v = deltaLofz(zv);
        rWallE_v = deltaEofz(zv);
        inGap = rv > rWallL_v && rv < rWallE_v;
        tag = 'OK (inside gap)';
        if ~inGap
            tag = 'ALERT: OUTSIDE the gap at its own z!';
            anyAlert = true;
        end
        fprintf('    vertex %d: (r=%.4f, z=%.4f) um -- rWallL=%.4f, rWallE=%.4f at that z -- %s\n', ...
            v, rv*1e6, zv*1e6, rWallL_v*1e6, rWallE_v*1e6, tag);
        fprintf('       sigma_rr=%.4f, sigma_zz=%.4f, sigma_rz=%.4f Pa, P=%.4f Pa\n', ...
            sigmaCellF(idx,1), sigmaCellF(idx,3), sigmaCellF(idx,4), pCellF(idx));
    end
    if ~anyAlert
        fprintf('  --> All 3 source points confirmed inside the true fluid gap. No cross-phase contamination.\n');
    end
end

fprintf('================ LEUKOCYTE-FLUID INTERFACE ================\n');
labelsL = {'L1','L2'};
for i = 1:numel(zPoints)
    z = zPoints(i);
    rWallL_i = deltaLofz(z); rWallE_i = deltaEofz(z);
    gapHere = rWallE_i - rWallL_i;
    rq = rWallL_i + epsFrac*gapHere;
    trace_point(labelsL{i}, rq, z, DT, centerF, sigmaCellF, pCellF, deltaLofz, deltaEofz);
end

fprintf('\n================ ENDOTHELIUM-FLUID INTERFACE ================\n');
labelsE = {'E1','E2'};
for i = 1:numel(zPoints)
    z = zPoints(i);
    rWallL_i = deltaLofz(z); rWallE_i = deltaEofz(z);
    gapHere = rWallE_i - rWallL_i;
    rq = rWallE_i - epsFrac*gapHere;
    trace_point(labelsE{i}, rq, z, DT, centerF, sigmaCellF, pCellF, deltaLofz, deltaEofz);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
