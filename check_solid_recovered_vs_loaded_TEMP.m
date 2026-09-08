clc;
S = load('out_pure2DMAC_fixedpoint.mat');
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

normFromTractionE = interp1(fl.tractionE.z(:), fl.tractionE.normal(:), endo.z(:), 'linear', 'extrap');

[maxPct, iMax] = max(endo.pctEn);
fprintf('Worst EnE query point: index %d, z=%.4f um, pctEn=%.2f%%\n', iMax, endo.z(iMax)*1e6, maxPct);
fprintf('  en (tnS - tnF)      = %.6e Pa\n', endo.en(iMax));
fprintf('  tnFluid (measured)  = %.6e Pa\n', endo.tnFluid(iMax));
fprintf('  tnS (solid, inferred = en+tnFluid) = %.6e Pa\n', endo.en(iMax) + endo.tnFluid(iMax));
fprintf('  normal(loaded, fluid.tractionE)    = %.6e Pa\n', normFromTractionE(iMax));
fprintf('  rWall (endo wall radius at this z) = %.6e m\n', endo.rWall(iMax));
fprintf('  slope (endo wall dr/dz)            = %.6e\n', endo.slope(iMax));

fprintf('\nFull table:\n');
fprintf('%-8s %-14s %-14s %-14s %-14s\n', 'z[um]', 'en', 'tnFluid', 'tnS(inferred)', 'loaded');
for i = 1:numel(endo.z)
    tnS = endo.en(i) + endo.tnFluid(i);
    fprintf('%-8.3f %-14.4f %-14.4f %-14.4f %-14.4f\n', ...
        endo.z(i)*1e6, endo.en(i), endo.tnFluid(i), tnS, normFromTractionE(i));
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
