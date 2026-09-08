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

fprintf('fl fieldnames:\n');
disp(fieldnames(fl));

if ~isfield(fl, 'tractionL')
    fprintf('NO tractionL field found on fl -- leukocyte side loaded differently, need to inspect solver.\n');
    fprintf('SMOKE_TEST_STATUS: NO_TRACTIONL\n');
    return;
end

disp('tractionL fields:');
disp(fieldnames(fl.tractionL));

tanFromTractionL = interp1(fl.tractionL.z(:), fl.tractionL.tangent(:), leuko.z(:), 'linear', 'extrap');
normFromTractionL = interp1(fl.tractionL.z(:), fl.tractionL.normal(:), leuko.z(:), 'linear', 'extrap');

fprintf('\n%-8s %-16s %-16s %-12s | %-16s %-16s\n', ...
    'z[um]', 'ttFluid(mismatch)', 'tangent(loadedL)', '%diff', 'tnFluid(mismatch)', 'normal(loadedL)');
for i = 1:5:numel(leuko.z)
    tPct = abs(leuko.ttFluid(i) - tanFromTractionL(i)) / max(abs(tanFromTractionL(i)), 1e-6) * 100;
    fprintf('%-8.3f %-16.4f %-16.4f %-12.2f | %-16.4f %-16.4f\n', ...
        leuko.z(i)*1e6, leuko.ttFluid(i), tanFromTractionL(i), tPct, leuko.tnFluid(i), normFromTractionL(i));
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
