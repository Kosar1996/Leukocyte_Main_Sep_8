clc;
S = load('out_pure2DMAC_signtest_FIXED.mat');
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

tanFromTractionE = interp1(fl.tractionE.z(:), fl.tractionE.tangent(:), endo.z(:), 'linear', 'extrap');
tanFromTractionL = interp1(fl.tractionL.z(:), fl.tractionL.tangent(:), leuko.z(:), 'linear', 'extrap');

fprintf('=== ENDOTHELIUM tangent: loaded (fluid.tractionE.tangent = -tauE) vs measured (endo.ttFluid) ===\n');
fprintf('%-8s %-16s %-16s %-12s\n', 'z[um]', 'ttFluid(meas)', 'tangent(loaded)', '%diff');
for i = 1:5:numel(endo.z)
    pct = abs(endo.ttFluid(i) - tanFromTractionE(i)) / max(abs(tanFromTractionE(i)), 1e-6) * 100;
    fprintf('%-8.3f %-16.4f %-16.4f %-12.2f\n', endo.z(i)*1e6, endo.ttFluid(i), tanFromTractionE(i), pct);
end

fprintf('\n=== LEUKOCYTE tangent: loaded (fluid.tractionL.tangent = +tauL) vs measured (leuko.ttFluid) ===\n');
fprintf('%-8s %-16s %-16s %-12s\n', 'z[um]', 'ttFluid(meas)', 'tangent(loaded)', '%diff');
for i = 1:5:numel(leuko.z)
    pct = abs(leuko.ttFluid(i) - tanFromTractionL(i)) / max(abs(tanFromTractionL(i)), 1e-6) * 100;
    fprintf('%-8.3f %-16.4f %-16.4f %-12.2f\n', leuko.z(i)*1e6, leuko.ttFluid(i), tanFromTractionL(i), pct);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
