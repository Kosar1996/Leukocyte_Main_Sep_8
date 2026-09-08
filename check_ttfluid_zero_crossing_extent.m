clc; close all;
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

fprintf('Leukocyte side: %d query points\n', numel(leuko.z));
fprintf('  ttFluid range: [%.4f, %.4f] Pa\n', min(leuko.ttFluid), max(leuko.ttFluid));
nNearZeroL = sum(abs(leuko.ttFluid) < 1.0);
fprintf('  Points with |ttFluid| < 1.0 Pa: %d of %d (%.0f%%)\n', nNearZeroL, numel(leuko.ttFluid), 100*nNearZeroL/numel(leuko.ttFluid));
signChangesL = sum(diff(sign(leuko.ttFluid)) ~= 0);
fprintf('  Sign changes across z (zero-crossings): %d\n\n', signChangesL);

fprintf('Endothelium side: %d query points\n', numel(endo.z));
fprintf('  ttFluid range: [%.4f, %.4f] Pa\n', min(endo.ttFluid), max(endo.ttFluid));
nNearZeroE = sum(abs(endo.ttFluid) < 1.0);
fprintf('  Points with |ttFluid| < 1.0 Pa: %d of %d (%.0f%%)\n', nNearZeroE, numel(endo.ttFluid), 100*nNearZeroE/numel(endo.ttFluid));
signChangesE = sum(diff(sign(endo.ttFluid)) ~= 0);
fprintf('  Sign changes across z (zero-crossings): %d\n', signChangesE);

fig = figure('Position',[100 100 900 600]);
subplot(2,1,1);
plot(leuko.z*1e6, leuko.ttFluid, '-o','MarkerSize',3); hold on; yline(0,'k--');
xlabel('z [\mum]'); ylabel('ttFluid (leuko side) [Pa]'); grid on;
title('Local tangential fluid traction vs z -- leukocyte side');
subplot(2,1,2);
plot(endo.z*1e6, endo.ttFluid, '-o','MarkerSize',3); hold on; yline(0,'k--');
xlabel('z [\mum]'); ylabel('ttFluid (endo side) [Pa]'); grid on;
title('Local tangential fluid traction vs z -- endothelium side');
saveas(fig, 'ttfluid_vs_z_check.png');

fprintf('\nSMOKE_TEST_STATUS: OK\n');
