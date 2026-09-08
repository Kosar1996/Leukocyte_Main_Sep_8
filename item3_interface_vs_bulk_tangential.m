%% ITEM3_INTERFACE_VS_BULK_TANGENTIAL

% "examine whether the interface tangential stresses agree closely with
% the tangential stresses in the bulk fluid and solid nodes closest to
% the interface."
%
% Compares sigma_rz AT the interface (what the correction actually
% samples) against sigma_rz in the BULK -- fluid at the gap midpoint
% (as far from either wall as this gap allows), solid several element
% widths deep into the material -- at the same z. Tested at the
% narrowest-gap location (step 40, z~3.75um, where her concern about a
% shrinking gap applies most directly) and, for comparison, at the wider
% z=3.45-3.50um points already used in the Issue #2 exercise.

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
FrzE = scatteredInterpolant(rE, zE, stressE.sigma_rz, 'linear', 'nearest');
FrzL = scatteredInterpolant(rL, zL, stressL.sigma_rz, 'linear', 'nearest');

meshF2 = add_fluid_nodes(fl.meshF);
[~, sigmaCellF, centerF] = recover_fluid_nodes_pressure_stress_Q4(meshF2, fl.ur2D, fl.uz2D, par.mu, fl.pCell);
FrzFl = scatteredInterpolant(centerF(:,1), centerF(:,2), sigmaCellF(:,4), 'linear', 'nearest');

zGrid = out.z(:);
deltaE_k = out.deltaEHist(:,k);
deltaL_k = out.deltaLHist(:,k);
gap = deltaE_k - deltaL_k;
[minGap, iMin] = min(gap);
zContact = zGrid(iMin);

epsFrac = 0.1;
% "bulk" solid depth: a fixed physical distance well beyond the epsFrac
% offset, representative of "a few element widths in", not interface-adjacent
bulkSolidDepth = 0.15e-6;

function report_row(label, z, deltaLofz, deltaEofz, FrzL, FrzE, FrzFl, epsFrac, bulkSolidDepth)
    rWallL = deltaLofz(z); rWallE = deltaEofz(z);
    gapHere = rWallE - rWallL;
    epsHere = epsFrac * gapHere;

    % Fluid: interface-adjacent (leukocyte side) vs bulk (gap midpoint)
    rzFl_interfaceL = FrzFl(rWallL + epsHere, z);
    rzFl_bulk       = FrzFl(0.5*(rWallL+rWallE), z);  % gap midpoint = fluid "bulk"

    % Fluid: interface-adjacent (endothelium side)
    rzFl_interfaceE = FrzFl(rWallE - epsHere, z);

    % Solid: interface-adjacent vs bulk (deeper into the material)
    rzL_interface = FrzL(rWallL - epsHere, z);
    rzL_bulk      = FrzL(max(rWallL - bulkSolidDepth, 0), z);
    rzE_interface = FrzE(rWallE + epsHere, z);
    rzE_bulk      = FrzE(rWallE + bulkSolidDepth, z);

    fprintf('\n--- %s: z=%.4f um, gap=%.2f nm ---\n', label, z*1e6, gapHere*1e9);
    fprintf('  LEUKOCYTE side: solid interface=%.4f Pa, solid bulk (%.0fnm deep)=%.4f Pa (diff=%.4f Pa)\n', ...
        rzL_interface, bulkSolidDepth*1e9, rzL_bulk, rzL_interface-rzL_bulk);
    fprintf('                  fluid interface=%.4f Pa, fluid bulk (gap midpoint)=%.4f Pa (diff=%.4f Pa)\n', ...
        rzFl_interfaceL, rzFl_bulk, rzFl_interfaceL-rzFl_bulk);
    fprintf('  ENDOTHELIUM side: solid interface=%.4f Pa, solid bulk (%.0fnm deep)=%.4f Pa (diff=%.4f Pa)\n', ...
        rzE_interface, bulkSolidDepth*1e9, rzE_bulk, rzE_interface-rzE_bulk);
    fprintf('                    fluid interface=%.4f Pa, fluid bulk (gap midpoint)=%.4f Pa (diff=%.4f Pa)\n', ...
        rzFl_interfaceE, rzFl_bulk, rzFl_interfaceE-rzFl_bulk);
end

deltaLofz = @(z) interp1(zGrid, deltaL_k, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(zGrid, deltaE_k, z, 'linear', 'extrap');

fprintf('================ Narrowest-gap location (her concern applies most directly) ================\n');
report_row('Narrowest gap (zContact)', zContact, deltaLofz, deltaEofz, FrzL, FrzE, FrzFl, epsFrac, bulkSolidDepth);

fprintf('\n================ Comparison points: z=3.45-3.50um (wider, from the Issue #2 exercise) ================\n');
report_row('z=3.45um', 3.45e-6, deltaLofz, deltaEofz, FrzL, FrzE, FrzFl, epsFrac, bulkSolidDepth);
report_row('z=3.50um', 3.50e-6, deltaLofz, deltaEofz, FrzL, FrzE, FrzFl, epsFrac, bulkSolidDepth);

fprintf('\nSMOKE_TEST_STATUS: OK\n');
