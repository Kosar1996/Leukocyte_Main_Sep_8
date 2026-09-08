clc;
S = load('out_debug_radialaxial.mat');
out = S.out;
k = out.stopStep;
st = out.stateHist{k};
par = out.par;
dtStep = out.dtHist(k);
fl = out.fluidHist{k};

parLmismatch = par;
if isfield(par,'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par,'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par,'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

mismatchOpts = struct('nQuery', 31, 'epsFrac', 0.1, 'trimFrac', 0.05);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);

% Interpolate tractionE.tangent (used to LOAD the solid solve) onto the
% same z query points used by the mismatch measurement (ttFluid).
tanFromTractionE = interp1(fl.tractionE.z(:), fl.tractionE.tangent(:), endo.z(:), 'linear', 'extrap');
normFromTractionE = interp1(fl.tractionE.z(:), fl.tractionE.normal(:), endo.z(:), 'linear', 'extrap');

fprintf('%-8s %-16s %-16s %-12s | %-16s %-16s\n', ...
    'z[um]', 'ttFluid(mismatch)', 'tangent(loadedE)', '%diff', 'tnFluid(mismatch)', 'normal(loadedE)');
for i = 1:5:numel(endo.z)
    tPct = abs(endo.ttFluid(i) - tanFromTractionE(i)) / max(abs(tanFromTractionE(i)), 1e-6) * 100;
    fprintf('%-8.3f %-16.4f %-16.4f %-12.2f | %-16.4f %-16.4f\n', ...
        endo.z(i)*1e6, endo.ttFluid(i), tanFromTractionE(i), tPct, endo.tnFluid(i), normFromTractionE(i));
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
