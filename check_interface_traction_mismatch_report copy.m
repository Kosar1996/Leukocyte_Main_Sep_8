function report = check_interface_traction_mismatch_report(out, plotstep, opts)
%CHECK_INTERFACE_TRACTION_MISMATCH_REPORT

% traction on the two sides of the interface (leukocyte|fluid and
% fluid|endothelium), rotate both into the local normal/tangential frame
% of the interface curve r = delta(z), and report the mismatch
%   en(s) = tn_solid(s) - tn_fluid(s)
%   et(s) = tt_solid(s) - tt_fluid(s)
% as a function of interface arc length s, along with max/median/90th
% percentile of |en| and |et|.
%
% This reuses the exact rotation/interpolation machinery already in
% check_interface_traction_continuity.m (same stress recovery, same
% viscoelastic solid stress via recover_nodal_stress_axisym_viscoelastic.m,
% same local-slope normal-vector construction), but:
%   1. evaluates on a DENSE z-grid instead of 5 spot checks, so percentile
%      statistics are meaningful;
%   2. converts z to physical arc length s along each side's own deformed
%      interface curve;
%   3. returns numeric arrays + summary stats instead of only printing;
%   4. optionally plots en(s), et(s) for both sides.
%
% Usage:
%   report = check_interface_traction_mismatch_report(out, plotstep)
%   report = check_interface_traction_mismatch_report(out, plotstep, opts)
%
% opts (all optional):
%   opts.nQuery   number of interface sample points (default 61)
%   opts.epsFrac  fraction of local gap width used to offset the sample
%                 point off the wall on each side (default 0.1, same
%                 default as check_interface_traction_continuity.m)
%   opts.trimFrac fraction of the z-range trimmed off each end before
%                 sampling, to stay inside the fluid mesh's interpolation
%                 domain (default 0.05)
%   opts.makePlot true/false (default true)
%   opts.label    string used in plot titles / printed header, e.g.
%                 '1D hybrid, step 1' (default '')
%
% report fields (each a struct with sub-structs .leuko and .endo):
%   report.leuko.s / .en / .et         arc length [m] and mismatch [Pa] arrays
%   report.leuko.stats.maxAbsEn / .medianAbsEn / .p90AbsEn
%   report.leuko.stats.maxAbsEt / .medianAbsEt / .p90AbsEt
%   report.endo.*                      same fields for the endothelium side
%   report.t, report.plotstep, report.label

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts,'nQuery')   || isempty(opts.nQuery),   opts.nQuery   = 61;   end
if ~isfield(opts,'epsFrac')  || isempty(opts.epsFrac),  opts.epsFrac  = 0.1;  end
if ~isfield(opts,'trimFrac') || isempty(opts.trimFrac), opts.trimFrac = 0.05; end
if ~isfield(opts,'makePlot') || isempty(opts.makePlot), opts.makePlot = true; end
if ~isfield(opts,'label'),    opts.label = '';                                end

fluid = out.fluidHist{plotstep};
mu = out.par.mu;
meshF = add_fluid_nodes(fluid.meshF);
[~, sigmaCellFluid, centerFluid] = recover_fluid_nodes_pressure_stress_Q4( ...
    meshF, fluid.ur2D, fluid.uz2D, mu, fluid.pCell);

st = out.stateHist{plotstep};
if ~isfield(st, 'uEPrev') || ~isfield(st, 'uLPrev')
    error(['check_interface_traction_mismatch_report: out.stateHist{%d} is missing ', ...
        'uEPrev/uLPrev -- cannot compute the viscous stress contribution.'], plotstep);
end
if isfield(out, 'dtHist') && numel(out.dtHist) >= plotstep && isfinite(out.dtHist(plotstep))
    dtStep = out.dtHist(plotstep);
else
    dtStep = out.par.dt;
    warning('out.dtHist not available for step %d; falling back to par.dt.', plotstep);
end

stressE = recover_nodal_stress_axisym_viscoelastic( ...
    out.meshE, st.uE, st.uEPrev, dtStep, out.par);

parL = out.par;
if isfield(out.par, 'GL')  && isfinite(out.par.GL),  parL.Ge = out.par.GL;  end
if isfield(out.par, 'KL')  && isfinite(out.par.KL),  parL.Ke = out.par.KL;  end
if isfield(out.par, 'etaL') && isfinite(out.par.etaL), parL.etaE = out.par.etaL; end
stressL = recover_nodal_stress_axisym_viscoelastic( ...
    out.meshL, st.uL, st.uLPrev, dtStep, parL);

uE = st.uE;
rE = out.meshE.nodes(:,1) + uE(1:2:end);
zE = out.meshE.nodes(:,2) + uE(2:2:end);
uL = st.uL;
rL = out.meshL.nodes(:,1) + uL(1:2:end);
zL = out.meshL.nodes(:,2) + uL(2:2:end);

interpF = @(x,y,v) scatteredInterpolant(x(:), y(:), v(:), 'linear', 'none');
FfluidRR = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,1));
FfluidZZ = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,3));
FfluidRZ = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,4));

FLrr = interpF(rL, zL, stressL.sigma_rr);
FLzz = interpF(rL, zL, stressL.sigma_zz);
FLrz = interpF(rL, zL, stressL.sigma_rz);
FErr = interpF(rE, zE, stressE.sigma_rr);
FEzz = interpF(rE, zE, stressE.sigma_zz);
FErz = interpF(rE, zE, stressE.sigma_rz);

mesh = fluid.meshF;
if ~all(isfield(mesh, {'zc','deltaL_c','deltaE_c'}))
    error('check_interface_traction_mismatch_report: needs fluid.meshF.zc/deltaL_c/deltaE_c.');
end

zLo = max(min(mesh.zc), min(zL));
zHi = min(max(mesh.zc), max(zL));
span = zHi - zLo;
zLoQ = zLo + opts.trimFrac * span;
zHiQ = zHi - opts.trimFrac * span;
zQuery = linspace(zLoQ, zHiQ, opts.nQuery);

deltaLofz = @(z) interp1(mesh.zc, mesh.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh.zc, mesh.deltaE_c, z, 'linear', 'extrap');

dz_fd = 1e-8;
slopeLofz = @(z) (deltaLofz(z+dz_fd) - deltaLofz(z-dz_fd)) / (2*dz_fd);
slopeEofz = @(z) (deltaEofz(z+dz_fd) - deltaEofz(z-dz_fd)) / (2*dz_fd);

nQ = numel(zQuery);
enL = nan(nQ,1); etL = nan(nQ,1);
enE = nan(nQ,1); etE = nan(nQ,1);
slopeLv = nan(nQ,1); slopeEv = nan(nQ,1);

for k = 1:nQ
    z = zQuery(k);
    rLwall = deltaLofz(z);
    rEwall = deltaEofz(z);
    gap = rEwall - rLwall;
    eps_ = opts.epsFrac * gap;

    slopeL = slopeLofz(z);
    slopeLv(k) = slopeL;
    nL = [1, -slopeL] / norm([1, -slopeL]);

    slopeE = slopeEofz(z);
    slopeEv(k) = slopeE;
    nE = [1, -slopeE] / norm([1, -slopeE]);

    rq  = rLwall + eps_;
    sigF = [FfluidRR(rq,z), FfluidRZ(rq,z); FfluidRZ(rq,z), FfluidZZ(rq,z)];
    rq2 = rLwall - eps_;
    sigS = [FLrr(rq2,z), FLrz(rq2,z); FLrz(rq2,z), FLzz(rq2,z)];
    [tnF, ttF] = traction_components(sigF, nL);
    [tnS, ttS] = traction_components(sigS, nL);
    enL(k) = tnS - tnF;
    etL(k) = ttS - ttF;

    rq  = rEwall - eps_;
    sigF = [FfluidRR(rq,z), FfluidRZ(rq,z); FfluidRZ(rq,z), FfluidZZ(rq,z)];
    rq2 = rEwall + eps_;
    sigS = [FErr(rq2,z), FErz(rq2,z); FErz(rq2,z), FEzz(rq2,z)];
    [tnF, ttF] = traction_components(sigF, nE);
    [tnS, ttS] = traction_components(sigS, nE);
    enE(k) = tnS - tnF;
    etE(k) = ttS - ttF;
end

sL = arc_length_from_slope(zQuery, slopeLv);
sE = arc_length_from_slope(zQuery, slopeEv);

report = struct();
report.plotstep = plotstep;
report.label = opts.label;
if isfield(out,'t') && numel(out.t) >= plotstep
    report.t = out.t(plotstep);
else
    report.t = NaN;
end

report.leuko = pack_side(sL, enL, etL, zQuery(:), deltaLofz(zQuery(:)));
report.endo  = pack_side(sE, enE, etE, zQuery(:), deltaEofz(zQuery(:)));

fprintf('\n=== Interface traction mismatch report: %s step %d (t=%.4g s) ===\n', ...
    opts.label, plotstep, report.t);
print_side('Leukocyte | fluid', report.leuko);
print_side('fluid | Endothelium', report.endo);

if opts.makePlot
    figure('Name', sprintf('Interface traction mismatch %s step %d', opts.label, plotstep));
    subplot(2,1,1);
    plot(sL*1e6, enL, 'o-', 'DisplayName', 'e_n leukocyte side'); hold on;
    plot(sE*1e6, enE, 's-', 'DisplayName', 'e_n endothelium side');
    yline(0, 'k:');
    xlabel('interface arc length s [\mum]'); ylabel('e_n = t_n^{solid} - t_n^{fluid}  [Pa]');
    legend('Location','best'); grid on;
    title(sprintf('Normal traction mismatch, %s step %d (t=%.4g s)', opts.label, plotstep, report.t));

    subplot(2,1,2);
    plot(sL*1e6, etL, 'o-', 'DisplayName', 'e_t leukocyte side'); hold on;
    plot(sE*1e6, etE, 's-', 'DisplayName', 'e_t endothelium side');
    yline(0, 'k:');
    xlabel('interface arc length s [\mum]'); ylabel('e_t = t_t^{solid} - t_t^{fluid}  [Pa]');
    legend('Location','best'); grid on;
    title('Tangential traction mismatch');
end

end

function side = pack_side(s, en, et, z, rWall)
side = struct();
side.s = s(:);
side.en = en(:);
side.et = et(:);
side.z = z(:);
side.rWall = rWall(:);

absEn = abs(en(isfinite(en)));
absEt = abs(et(isfinite(et)));

side.stats = struct();
if isempty(absEn)
    side.stats.maxAbsEn = NaN; side.stats.medianAbsEn = NaN; side.stats.p90AbsEn = NaN;
else
    side.stats.maxAbsEn = max(absEn);
    side.stats.medianAbsEn = median(absEn);
    side.stats.p90AbsEn = prctile(absEn, 90);
end
if isempty(absEt)
    side.stats.maxAbsEt = NaN; side.stats.medianAbsEt = NaN; side.stats.p90AbsEt = NaN;
else
    side.stats.maxAbsEt = max(absEt);
    side.stats.medianAbsEt = median(absEt);
    side.stats.p90AbsEt = prctile(absEt, 90);
end

[~, iEn] = max(abs(en(:)));
[~, iEt] = max(abs(et(:)));
side.stats.locMaxEn = [side.rWall(iEn), side.z(iEn)];
side.stats.locMaxEt = [side.rWall(iEt), side.z(iEt)];
end

function print_side(name, side)
fprintf('  %-20s |en|: max=%9.4f  median=%9.4f  p90=%9.4f  Pa   (peak at r=%.3f, z=%.3f um)\n', ...
    name, side.stats.maxAbsEn, side.stats.medianAbsEn, side.stats.p90AbsEn, ...
    side.stats.locMaxEn(1)*1e6, side.stats.locMaxEn(2)*1e6);
fprintf('  %-20s |et|: max=%9.4f  median=%9.4f  p90=%9.4f  Pa   (peak at r=%.3f, z=%.3f um)\n', ...
    '', side.stats.maxAbsEt, side.stats.medianAbsEt, side.stats.p90AbsEt, ...
    side.stats.locMaxEt(1)*1e6, side.stats.locMaxEt(2)*1e6);
end

function s = arc_length_from_slope(zQuery, slope)
% Cumulative physical arc length along a curve r = delta(z), given the
% local slope d(delta)/dz sampled at zQuery. s(1) = 0.
ds = sqrt(1 + slope(:).^2);
dz = diff(zQuery(:));
segMean = 0.5 * (ds(1:end-1) + ds(2:end));
s = [0; cumsum(segMean .* dz)];
end

function [tn, tt] = traction_components(sig, n)
n = n(:);
that = [-n(2); n(1)];
t = sig * n;
tn = t.' * n;
tt = t.' * that;
end
