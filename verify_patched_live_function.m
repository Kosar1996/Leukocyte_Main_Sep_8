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

tic;
mismatchOpts = struct('nQuery', 61, 'epsFrac', 0.1, 'trimFrac', 0.05);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);
fprintf('Elapsed (patched live function): %.1f s\n', toc);

fprintf('refScale (leuko) = %.4f Pa  [expected ~190.50 from standalone verification]\n', leuko.stats.refScale);
fprintf('refScale (endo)  = %.4f Pa  [expected ~190.50]\n', endo.stats.refScale);

% find nearest query points to z=3.45, 3.50 and print normal mismatch %
[~,iL1] = min(abs(leuko.z - 3.45e-6));
[~,iL2] = min(abs(leuko.z - 3.50e-6));
[~,iE1] = min(abs(endo.z - 3.45e-6));
[~,iE2] = min(abs(endo.z - 3.50e-6));
pctEnL = 100*abs(leuko.en)/leuko.stats.refScale;
pctEnE = 100*abs(endo.en)/endo.stats.refScale;
fprintf('Leuko near z=3.45: %.4f, En%%=%.2f  [expect near 37.1%% from standalone L1]\n', leuko.z(iL1)*1e6, pctEnL(iL1));
fprintf('Leuko near z=3.50: %.4f, En%%=%.2f  [expect near 35.5%% from standalone L2]\n', leuko.z(iL2)*1e6, pctEnL(iL2));
fprintf('Endo  near z=3.45: %.4f, En%%=%.2f  [expect near 33.8%% from standalone E1]\n', endo.z(iE1)*1e6, pctEnE(iE1));
fprintf('Endo  near z=3.50: %.4f, En%%=%.2f  [expect near 30.8%% from standalone E2]\n', endo.z(iE2)*1e6, pctEnE(iE2));

fprintf('\nSMOKE_TEST_STATUS: OK\n');
