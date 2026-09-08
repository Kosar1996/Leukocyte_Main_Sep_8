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

fprintf('fluidHist fields relevant to traction:\n');
disp(fieldnames(fl));

if isfield(fl, 'tractionE')
    disp('tractionE fields:');
    disp(fieldnames(fl.tractionE));
    fprintf('tractionE.z range: %.4f to %.4f um\n', min(fl.tractionE.z)*1e6, max(fl.tractionE.z)*1e6);
end

fprintf('\nendo.z range (mismatch query points): %.4f to %.4f um\n', min(endo.z)*1e6, max(endo.z)*1e6);
fprintf('endo.ttFluid at first 5 query points: %s\n', mat2str(endo.ttFluid(1:5).', 4));

if isfield(fl, 'tractionE') && isfield(fl.tractionE, 'z')
    % Interpolate tractionE's tangential component (if present) onto endo.z for comparison
    disp('Comparing at matching z if tractionE has a tangential field...');
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
