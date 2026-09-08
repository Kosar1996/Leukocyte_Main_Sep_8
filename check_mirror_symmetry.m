function check_mirror_symmetry(out, plotstep, zCenter, dzList)
%CHECK_MIRROR_SYMMETRY

% pressure difference (pIn=pOut=0 in this base case), the squeeze is
% geometrically symmetric about the leukocyte's waist, so:
%   sigma_rr, sigma_tt, sigma_zz  should be MIRROR-symmetric about that
%       z-plane: value at (zCenter+dz) should equal value at (zCenter-dz).
%   sigma_rz (shear) should be ANTI-symmetric: value at (zCenter+dz)
%       should equal MINUS the value at (zCenter-dz), i.e. their SUM
%       should be ~0.
%
% CHANGED from the first version:
%   1. zCenter, if not supplied, is now computed EXACTLY from the
%      leukocyte's own deformed mesh (midpoint of its z-extent) instead
%      of being eyeballed off a contour plot. You can still pass your own
%      zCenter to override this.
%   2. Each domain (fluid, endothelium, leukocyte) is now sampled at its
%      OWN appropriate radius, since they occupy very different radial
%      ranges (endothelium sits near the outer wall, ~r=4um; leukocyte
%      interior is smaller r; fluid is the thin gap between them). Using
%      one single radius for all three (as the first version did) is why
%      the endothelium rows came back all-NaN last time.
%
%   check_mirror_symmetry(out, plotstep)
%   check_mirror_symmetry(out, plotstep, zCenter)
%   check_mirror_symmetry(out, plotstep, zCenter, dzList)
%
% dzList: vector of z-offsets (m) to test, default [1 2 3]*1e-6.

if nargin < 4 || isempty(dzList)
    dzList = [1 2 3] * 1e-6;
end

uL0 = out.stateHist{plotstep}.uL;
rL0 = out.meshL.nodes(:,1) + uL0(1:2:end);
zL0 = out.meshL.nodes(:,2) + uL0(2:2:end);

if nargin < 3 || isempty(zCenter)
    zCenter = 0.5 * (min(zL0) + max(zL0));
    fprintf('zCenter not supplied -- using leukocyte deformed-mesh midpoint: %.4f um\n', ...
        zCenter*1e6);
end

fluid = out.fluidHist{plotstep};
mu = out.par.mu;

meshF = add_fluid_nodes(fluid.meshF);
[~, sigmaCellFluid, centerFluid] = recover_fluid_nodes_pressure_stress_Q4( ...
    meshF, fluid.ur2D, fluid.uz2D, mu, fluid.pCell);

stressE = recover_nodal_stress_axisym(out.meshE, out.stateHist{plotstep}.uE, out.par);
parL = out.par;
if isfield(out.par, 'GL') && isfinite(out.par.GL), parL.Ge = out.par.GL; end
if isfield(out.par, 'KL') && isfinite(out.par.KL), parL.Ke = out.par.KL; end
stressL = recover_nodal_stress_axisym(out.meshL, out.stateHist{plotstep}.uL, parL);

uE = out.stateHist{plotstep}.uE;
rE = out.meshE.nodes(:,1) + uE(1:2:end);
zE = out.meshE.nodes(:,2) + uE(2:2:end);

rL = rL0;
zL = zL0;

interpF = @(x, y, v) scatteredInterpolant(x(:), y(:), v(:), 'linear', 'none');

Ffluid = struct( ...
    'rr', interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,1)), ...
    'tt', interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,2)), ...
    'zz', interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,3)), ...
    'rz', interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,4)));

FE = struct( ...
    'rr', interpF(rE, zE, stressE.sigma_rr), ...
    'tt', interpF(rE, zE, stressE.sigma_tt), ...
    'zz', interpF(rE, zE, stressE.sigma_zz), ...
    'rz', interpF(rE, zE, stressE.sigma_rz));

FL = struct( ...
    'rr', interpF(rL, zL, stressL.sigma_rr), ...
    'tt', interpF(rL, zL, stressL.sigma_tt), ...
    'zz', interpF(rL, zL, stressL.sigma_zz), ...
    'rz', interpF(rL, zL, stressL.sigma_rz));

% Per-domain sample radius: midpoint of each domain's own actual r-range,
% so we're guaranteed to land inside real data instead of guessing one
% radius for all three very differently-located domains.
rFluidQuery = mean([min(centerFluid(:,1)), max(centerFluid(:,1))]);
rEQuery     = mean([min(rE), max(rE)]);
rLQuery     = mean([min(rL), max(rL)]);

fprintf('\n=== Mirror symmetry check: step %d (t=%.4g s), zCenter=%.4f um ===\n', ...
    plotstep, out.t(plotstep), zCenter*1e6);

fprintf('\n--- fluid (r=%.3f um) ---\n', rFluidQuery*1e6);
run_checks(Ffluid, rFluidQuery, zCenter, dzList);

fprintf('\n--- endothelium (r=%.3f um) ---\n', rEQuery*1e6);
run_checks(FE, rEQuery, zCenter, dzList);

fprintf('\n--- leukocyte (r=%.3f um) ---\n', rLQuery*1e6);
run_checks(FL, rLQuery, zCenter, dzList);

end

function run_checks(F, r, zCenter, dzList)
for dz = dzList
    zp = zCenter + dz;
    zm = zCenter - dz;
    report_component('rr', F.rr, r, zp, zm, dz, false);
    report_component('tt', F.tt, r, zp, zm, dz, false);
    report_component('zz', F.zz, r, zp, zm, dz, false);
    report_component('rz', F.rz, r, zp, zm, dz, true);
end
end

function report_component(label, F, r, zp, zm, dz, isShear)
vp = F(r, zp);
vm = F(r, zm);
if ~isfinite(vp) || ~isfinite(vm)
    fprintf('  dz=%.2fum  sigma_%-3s : +side or -side outside data range (NaN) -- try a smaller dz\n', ...
        dz*1e6, label);
    return;
end
if isShear
    resid = vp + vm;   % should be ~0 for antisymmetric shear
    fprintf('  dz=%.2fum  sigma_%-3s : +side=%9.3f  -side=%9.3f  SUM (want ~0) = %9.3f\n', ...
        dz*1e6, label, vp, vm, resid);
else
    resid = vp - vm;   % should be ~0 for symmetric normal stress
    fprintf('  dz=%.2fum  sigma_%-3s : +side=%9.3f  -side=%9.3f  DIFF (want ~0) = %9.3f\n', ...
        dz*1e6, label, vp, vm, resid);
end
end
