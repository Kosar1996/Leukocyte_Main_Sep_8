%% INVESTIGATE_STUCK_INTERFACE_DOF
% Looks for the anomaly the other session's investigation flagged as its
% "current lead": one specific interface degree of freedom where the
% applied traction changes measurably between correction passes, but the
% displacement comes back bit-for-bit identical every time.
%
% Method: run the same partitioned correction pass apply_bodyfitted_MAC_
% traction_correction_feedback.m uses (solid re-solve under current fluid
% traction, relax, fluid re-solve), but instrumented to save the FULL
% uE/uL and tractionE/tractionL at every single pass -- not just the
% aggregate %mismatch summary. Then, per interface node, compare:
%   - how much the applied traction at that node's z-location changed
%     pass-to-pass
%   - how much that node's displacement (r and z DOFs separately) changed
%     pass-to-pass
% A node with large traction change but ~zero displacement change across
% MULTIPLE passes (not just one) is the signature being hunted for.

clc; close all;
cd(fileparts(mfilename('fullpath')));

savedRunFile = 'out_1D_t10_for_review.mat';
S = load(savedRunFile);
runFieldName = fieldnames(S);
out = S.(runFieldName{1});

SE = load(out.par.endotheliumPrestressFile);
SL = load(out.par.leukocytePrestressFile);
baseE = SE.baseE;
baseL = SL.baseL;

parL = out.par;
if isfield(out.par, 'EL'), parL.Ee = out.par.EL; end
if isfield(out.par, 'nuL'), parL.nuE = out.par.nuL; end
if isfield(out.par, 'GL')
    parL.Ge = out.par.GL;
elseif isfield(parL, 'Ee') && isfield(parL, 'nuE')
    parL.Ge = parL.Ee/(2*(1+parL.nuE));
end
if isfield(out.par, 'KL')
    parL.Ke = out.par.KL;
elseif isfield(parL, 'Ee') && isfield(parL, 'nuE')
    parL.Ke = parL.Ee/(3*(1-2*parL.nuE));
end
if isfield(out.par, 'etaL'), parL.etaE = out.par.etaL; end
if isfield(out.par, 'useViscoelasticLeukocyte')
    parL.useViscoelasticEndothelium = out.par.useViscoelasticLeukocyte;
end

k = out.stopStep;
z = out.z;
par = out.par;
par.dt = out.dtHist(k);

% Same "genuine room to improve" scenario as demo_traction_correction_
% feedback.m's Scenario B: older solid shape paired with the later step's
% fluid traction, so the correction has real work to do across several
% passes (needed to see a "stuck across multiple passes" pattern, not
% just a one-off).
old = out.stateHist{k-1};
state = out.stateHist{k-1};
fluid = out.fluidHist{k};

nPasses = 8;
relax = 1.0;
if isfield(par,'bodyFittedTractionCorrectionRelax') && isfinite(par.bodyFittedTractionCorrectionRelax)
    relax = min(1.0, max(0.0, par.bodyFittedTractionCorrectionRelax));
end

useRLoutInner = use_RLout_fluid_interface_for_solid_leukocyte(par);
hasL = ~useRLoutInner && ~isempty(out.meshL) && ~isempty(out.interfaceL) && ...
    isfield(state,'uL') && ~isempty(state.uL);

nE = numel(state.uE);
nL = numel(state.uL);
uEHist = nan(nE, nPasses+1);
uLHist = nan(nL, nPasses+1);
uEHist(:,1) = state.uE;
uLHist(:,1) = state.uL;

% Traction sampled at each interface mesh node's own z-location, per pass.
idsE = out.interfaceE(:);
idsL = out.interfaceL(:);
zE_nodes = out.meshE.nodes(idsE,2);   % reference z of endothelium interface nodes
zL_nodes = out.meshL.nodes(idsL,2);   % reference z of leukocyte interface nodes

tnEHist = nan(numel(idsE), nPasses+1);
ttEHist = nan(numel(idsE), nPasses+1);
tnLHist = nan(numel(idsL), nPasses+1);
ttLHist = nan(numel(idsL), nPasses+1);

fprintf('Running %d instrumented correction passes...\n', nPasses);
for ic = 1:nPasses
    if ~isfield(fluid,'tractionE') || ~isfield(fluid,'tractionL')
        [fluid.tractionL, fluid.tractionE] = compute_bodyfitted_wall_traction(fluid.meshF, fluid, par);
    end

    % Record the CURRENT applied traction (before this pass), sampled at
    % each interface node's z.
    tnEHist(:,ic) = interp1(fluid.tractionE.z(:), fluid.tractionE.normal(:), zE_nodes, 'linear', 'extrap');
    ttEHist(:,ic) = interp1(fluid.tractionE.z(:), fluid.tractionE.tangent(:), zE_nodes, 'linear', 'extrap');
    if hasL
        tnLHist(:,ic) = interp1(fluid.tractionL.z(:), fluid.tractionL.normal(:), zL_nodes, 'linear', 'extrap');
        ttLHist(:,ic) = interp1(fluid.tractionL.z(:), fluid.tractionL.tangent(:), zL_nodes, 'linear', 'extrap');
    end

    uEold = state.uE;
    uEcorr = solve_finite_def_solid(out.meshE, old.uE, fluid.tractionE, ...
        out.interfaceE, baseE, par.supportE, par, state.uE);
    state.uE = uEold + relax * (uEcorr - uEold);
    [state.deltaE, state.UwE] = monolithic_interface_kinematics_value_only( ...
        out.meshE, state.uE, old.uE, out.interfaceE, z, par);

    if hasL
        uLold = state.uL;
        uLcorr = solve_finite_def_solid(out.meshL, old.uL, fluid.tractionL, ...
            out.interfaceL, baseL, par.supportL, parL, state.uL);
        state.uL = uLold + relax * (uLcorr - uLold);
        [state.deltaL, state.UwL] = monolithic_interface_kinematics_value_only( ...
            out.meshL, state.uL, old.uL, out.interfaceL, z, par);
    end

    state = attach_physical_solid_interface_fields( ...
        state, old, out.meshE, out.interfaceE, out.meshL, out.interfaceL, z, par);

    [fluid, okFluid, reasonFluid] = solve_selected_poststep_fluid(z, old, state, par);
    if ~okFluid
        fprintf('  pass %d: fluid re-solve failed (%s), stopping early.\n', ic, reasonFluid);
        break;
    end

    uEHist(:,ic+1) = state.uE;
    uLHist(:,ic+1) = state.uL;
    fprintf('  pass %d done: max|duE|=%.3e, max|duL|=%.3e\n', ic, ...
        max(abs(uEHist(:,ic+1)-uEHist(:,ic))), max(abs(uLHist(:,ic+1)-uLHist(:,ic))));
end

% Also record the FINAL traction (after the last pass), so we have
% nPasses+1 traction samples to match nPasses+1 displacement samples.
if ~isfield(fluid,'tractionE') || ~isfield(fluid,'tractionL')
    [fluid.tractionL, fluid.tractionE] = compute_bodyfitted_wall_traction(fluid.meshF, fluid, par);
end
tnEHist(:,nPasses+1) = interp1(fluid.tractionE.z(:), fluid.tractionE.normal(:), zE_nodes, 'linear', 'extrap');
ttEHist(:,nPasses+1) = interp1(fluid.tractionE.z(:), fluid.tractionE.tangent(:), zE_nodes, 'linear', 'extrap');
if hasL
    tnLHist(:,nPasses+1) = interp1(fluid.tractionL.z(:), fluid.tractionL.normal(:), zL_nodes, 'linear', 'extrap');
    ttLHist(:,nPasses+1) = interp1(fluid.tractionL.z(:), fluid.tractionL.tangent(:), zL_nodes, 'linear', 'extrap');
end

%% Per-node analysis: leukocyte interface
fprintf('\n===================== LEUKOCYTE interface nodes =====================\n');
find_stuck_nodes(idsL, zL_nodes, out.meshL.nodes(idsL,1), uLHist, tnLHist, ttLHist, 'Leukocyte');

%% Per-node analysis: endothelium interface
fprintf('\n===================== ENDOTHELIUM interface nodes =====================\n');
find_stuck_nodes(idsE, zE_nodes, out.meshE.nodes(idsE,1), uEHist, tnEHist, ttEHist, 'Endothelium');

save('stuck_dof_investigation.mat', 'uEHist', 'uLHist', 'tnEHist', 'ttEHist', 'tnLHist', 'ttLHist', ...
    'idsE', 'idsL', 'zE_nodes', 'zL_nodes');
fprintf('\nSaved raw per-pass history to stuck_dof_investigation.mat\n');

function find_stuck_nodes(ids, zNodes, rNodes, uHist, tnHist, ttHist, label)
nNodes = numel(ids);
nPasses1 = size(uHist,2);

maxDispDeltaR = nan(nNodes,1);
maxDispDeltaZ = nan(nNodes,1);
maxTracDeltaN = nan(nNodes,1);
maxTracDeltaT = nan(nNodes,1);

for i = 1:nNodes
    rDof = 2*ids(i) - 1;
    zDof = 2*ids(i);
    rSeries = uHist(rDof, :);
    zSeries = uHist(zDof, :);
    if any(isnan(rSeries)) || any(isnan(zSeries))
        continue;
    end
    maxDispDeltaR(i) = max(abs(diff(rSeries)));
    maxDispDeltaZ(i) = max(abs(diff(zSeries)));
    maxTracDeltaN(i) = max(abs(diff(tnHist(i,1:nPasses1))));
    maxTracDeltaT(i) = max(abs(diff(ttHist(i,1:nPasses1))));
end

% "Stuck" = traction changed by a meaningful amount (> tracThresh, in Pa)
% but displacement barely moved at all (< dispThresh, in meters) across
% EVERY consecutive pair of passes (max|diff| small means it never moved
% much on any single pass, not just net-zero over the whole run).
tracThresh = 1.0;      % Pa -- "meaningfully changing" load
dispThresh = 1e-13;    % m  -- effectively frozen (near machine precision
                        % relative to displacement scales in this problem)

fprintf('%-6s %-10s %-10s %-14s %-14s %-14s %-14s\n', ...
    'node', 'z [um]', 'r [um]', 'max|dTn| Pa', 'max|dTt| Pa', 'max|dur| m', 'max|duz| m');
stuckCount = 0;
for i = 1:nNodes
    flagR = maxTracDeltaN(i) > tracThresh && maxDispDeltaR(i) < dispThresh;
    flagZ = maxTracDeltaT(i) > tracThresh && maxDispDeltaZ(i) < dispThresh;
    if flagR || flagZ
        stuckCount = stuckCount + 1;
        fprintf('%-6d %-10.4f %-10.4f %-14.4g %-14.4g %-14.4g %-14.4g  <== STUCK%s%s\n', ...
            ids(i), zNodes(i)*1e6, rNodes(i)*1e6, maxTracDeltaN(i), maxTracDeltaT(i), ...
            maxDispDeltaR(i), maxDispDeltaZ(i), string_or(flagR,' (radial)',''), string_or(flagZ,' (axial)',''));
    end
end

if stuckCount == 0
    fprintf('No stuck nodes found on %s side by these thresholds (tracThresh=%.2f Pa, dispThresh=%.2g m).\n', ...
        label, tracThresh, dispThresh);
    fprintf('Top 5 nodes by max traction change, for reference:\n');
    [~, ord] = sort(max(maxTracDeltaN, maxTracDeltaT), 'descend', 'MissingPlacement','last');
    for kk = 1:min(5,numel(ord))
        i = ord(kk);
        fprintf('  node %-6d z=%-8.4f r=%-8.4f max|dTn|=%-10.4g max|dTt|=%-10.4g max|dur|=%-10.4g max|duz|=%-10.4g\n', ...
            ids(i), zNodes(i)*1e6, rNodes(i)*1e6, maxTracDeltaN(i), maxTracDeltaT(i), maxDispDeltaR(i), maxDispDeltaZ(i));
    end
else
    fprintf('\n%d stuck node(s) found on %s side.\n', stuckCount, label);
end
end

function out = string_or(cond, a, b)
if cond, out = a; else, out = b; end
end
