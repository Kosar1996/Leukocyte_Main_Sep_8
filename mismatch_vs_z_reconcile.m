%% MISMATCH_VS_Z_RECONCILE

% at the current/final state -- not per correction pass. Check whether
% L1's reported 36.1% is consistent with (an entry in / at most the max
% of) the full mismatch-vs-z profile from the SAME run, and check the
% trend vs z (increase-then-decrease near the thin region, or the
% opposite, monotonic-with-gap-thinning trend).

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

mismatchOpts = struct('nQuery', 61, 'epsFrac', 0.1, 'trimFrac', 0.05);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);

pctEnL = 100*abs(leuko.en)/leuko.stats.refScale;
pctEtL = 100*abs(leuko.et)/leuko.stats.refScale;
pctEnE = 100*abs(endo.en)/endo.stats.refScale;
pctEtE = 100*abs(endo.et)/endo.stats.refScale;

fprintf('================ LEUKOCYTE interface: %%mismatch vs z (%d query points, this run, final state) ================\n', numel(leuko.z));
fprintf('%8s %8s %8s\n', 'z[um]', 'En[%]', 'Et[%]');
for i = 1:numel(leuko.z)
    fprintf('%8.4f %8.2f %8.2f\n', leuko.z(i)*1e6, pctEnL(i), pctEtL(i));
end
fprintf('\nMax En (leukocyte, whole interface) = %.2f%% at z=%.4f um\n', max(pctEnL), leuko.z(find(pctEnL==max(pctEnL),1))*1e6);
fprintf('Max Et (leukocyte, whole interface) = %.2f%% at z=%.4f um\n', max(pctEtL), leuko.z(find(pctEtL==max(pctEtL),1))*1e6);

fprintf('\n================ ENDOTHELIUM interface: %%mismatch vs z ================\n');
fprintf('%8s %8s %8s\n', 'z[um]', 'En[%]', 'Et[%]');
for i = 1:numel(endo.z)
    fprintf('%8.4f %8.2f %8.2f\n', endo.z(i)*1e6, pctEnE(i), pctEtE(i));
end
fprintf('\nMax En (endothelium, whole interface) = %.2f%% at z=%.4f um\n', max(pctEnE), endo.z(find(pctEnE==max(pctEnE),1))*1e6);
fprintf('Max Et (endothelium, whole interface) = %.2f%% at z=%.4f um\n', max(pctEtE), endo.z(find(pctEtE==max(pctEtE),1))*1e6);

fprintf('\n================ Reconcile against report L1=36.1%%, L2=35.5%%, E1=32.5%%, E2=30.7%% (normal mismatch) ================\n');
[~, iZ345] = min(abs(leuko.z - 3.45e-6));
[~, iZ350] = min(abs(leuko.z - 3.50e-6));
fprintf('Leukocyte list, nearest to z=3.45um: z=%.4f, En=%.2f%% (report L1 said 36.1%%)\n', leuko.z(iZ345)*1e6, pctEnL(iZ345));
fprintf('Leukocyte list, nearest to z=3.50um: z=%.4f, En=%.2f%% (report L2 said 35.5%%)\n', leuko.z(iZ350)*1e6, pctEnL(iZ350));
[~, iZ345E] = min(abs(endo.z - 3.45e-6));
[~, iZ350E] = min(abs(endo.z - 3.50e-6));
fprintf('Endothelium list, nearest to z=3.45um: z=%.4f, En=%.2f%% (report E1 said 32.5%%)\n', endo.z(iZ345E)*1e6, pctEnE(iZ345E));
fprintf('Endothelium list, nearest to z=3.50um: z=%.4f, En=%.2f%% (report E2 said 30.7%%)\n', endo.z(iZ350E)*1e6, pctEnE(iZ350E));

fprintf('\nSMOKE_TEST_STATUS: OK\n');
