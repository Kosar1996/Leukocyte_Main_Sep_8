clc;
file = 'out_2D_unified_domain_40step.mat';
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

fprintf('Fields present: %s\n', strjoin(fieldnames(leuko.stats), ', '));
fprintf('Leuko: pctMaxEn=%.2f%% pctP90En=%.2f%% pctMedianEn=%.2f%%\n', leuko.stats.pctMaxEn, leuko.stats.pctP90En, leuko.stats.pctMedianEn);
fprintf('Leuko: pctMaxEt=%.2f%% pctP90Et=%.2f%% pctMedianEt=%.2f%%\n', leuko.stats.pctMaxEt, leuko.stats.pctP90Et, leuko.stats.pctMedianEt);
fprintf('Endo:  pctMaxEn=%.2f%% pctP90En=%.2f%% pctMedianEn=%.2f%%\n', endo.stats.pctMaxEn, endo.stats.pctP90En, endo.stats.pctMedianEn);
fprintf('Endo:  pctMaxEt=%.2f%% pctP90Et=%.2f%% pctMedianEt=%.2f%%\n', endo.stats.pctMaxEt, endo.stats.pctP90Et, endo.stats.pctMedianEt);
fprintf('\nSMOKE_TEST_STATUS: OK\n');
