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

fprintf('%-8s %-14s %-14s %-10s | %-14s %-14s %-10s\n', ...
    'z[um]', 'endo.et', 'endo.ttFluid', 'pctEt', 'leuko.et', 'leuko.ttFluid', 'pctEt(L)');
for i = 1:numel(endo.z)
    fprintf('%-8.3f %-14.6e %-14.6e %-10.2f | %-14.6e %-14.6e %-10.2f\n', ...
        endo.z(i)*1e6, endo.et(i), endo.ttFluid(i), endo.pctEt(i), ...
        leuko.et(i), leuko.ttFluid(i), leuko.pctEt(i));
end

fprintf('\nendo.ttFluid sign changes at index: ');
disp(find(diff(sign(endo.ttFluid)) ~= 0).');
fprintf('leuko.ttFluid sign changes at index: ');
disp(find(diff(sign(leuko.ttFluid)) ~= 0).');

[maxEt, iE] = max(endo.pctEt);
[maxLt, iL] = max(leuko.pctEt);
fprintf('\nmax endo pctEt = %.2f at z=%.3f um (ttFluid=%.4e, et=%.4e)\n', maxEt, endo.z(iE)*1e6, endo.ttFluid(iE), endo.et(iE));
fprintf('max leuko pctEt = %.2f at z=%.3f um (ttFluid=%.4e, et=%.4e)\n', maxLt, leuko.z(iL)*1e6, leuko.ttFluid(iL), leuko.et(iL));

fprintf('\nSMOKE_TEST_STATUS: OK\n');
