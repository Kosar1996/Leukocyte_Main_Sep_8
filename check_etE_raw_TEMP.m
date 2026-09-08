clc;
S = load('out_dtsmall_2steps.mat');
out = S.out;
k = out.stopStep;
st = out.stateHist{k};
par = out.par;
dtStep = out.dtHist(k);

parLmismatch = par;
if isfield(par,'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par,'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par,'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

fl = out.fluidHist{k};
mismatchOpts = struct('nQuery', 31, 'epsFrac', 0.1, 'trimFrac', 0.05);

[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);

[maxPctEt, iMax] = max(endo.pctEt);
fprintf('Endothelium tangential mismatch:\n');
fprintf('  max pctEt = %.2f%% at query index %d\n', maxPctEt, iMax);
fprintf('  raw solid-fluid mismatch (et, Pa) at that point: need to inspect endo fields\n');
disp(fieldnames(endo));

% Try to find the raw absolute mismatch and fluid traction arrays
if isfield(endo, 'et')
    fprintf('  et(iMax) [Pa] = %.6e\n', endo.et(iMax));
end
if isfield(endo, 'ttFluid')
    fprintf('  ttFluid(iMax) [Pa] = %.6e\n', endo.ttFluid(iMax));
end

fprintf('\npctEt distribution: min=%.2f median=%.2f p90=%.2f max=%.2f\n', ...
    min(endo.pctEt), median(endo.pctEt), prctile(endo.pctEt,90), max(endo.pctEt));

fprintf('\nSMOKE_TEST_STATUS: OK\n');
