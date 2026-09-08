%% PLOT_MISMATCH_VS_Z_NARROW
% Visualizes %mismatch vs z across the FULL narrow-domain interface
% (not just the 4 thin-region points used in the report), showing the
% low-at-edges / high-plateau-in-the-bulk pattern.

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

reportZ = [3.45, 3.50, 3.45, 3.50];
reportVal = [36.1, 35.5, 32.5, 30.7];
reportLabel = {'L1','L2','E1','E2'};

fig = figure('Position',[100 100 1000 700]);

subplot(2,1,1);
plot(leuko.z*1e6, pctEnL, 'b-o', 'MarkerSize',3, 'DisplayName','Normal (En)'); hold on;
plot(leuko.z*1e6, pctEtL, 'r-s', 'MarkerSize',3, 'DisplayName','Tangential (Et)');
plot(reportZ(1:2), reportVal(1:2), 'k*', 'MarkerSize',14, 'LineWidth',2, 'DisplayName','Report points (L1,L2)');
text(reportZ(1), reportVal(1)+2, 'L1', 'FontWeight','bold');
text(reportZ(2), reportVal(2)+2, 'L2', 'FontWeight','bold');
xlabel('z [\mum]'); ylabel('% mismatch');
title('Leukocyte-fluid interface: %mismatch vs z (full interface, narrow domain)');
legend('Location','best'); grid on;

subplot(2,1,2);
plot(endo.z*1e6, pctEnE, 'b-o', 'MarkerSize',3, 'DisplayName','Normal (En)'); hold on;
plot(endo.z*1e6, pctEtE, 'r-s', 'MarkerSize',3, 'DisplayName','Tangential (Et)');
plot(reportZ(3:4), reportVal(3:4), 'k*', 'MarkerSize',14, 'LineWidth',2, 'DisplayName','Report points (E1,E2)');
text(reportZ(3), reportVal(3)+2, 'E1', 'FontWeight','bold');
text(reportZ(4), reportVal(4)+2, 'E2', 'FontWeight','bold');
xlabel('z [\mum]'); ylabel('% mismatch');
title('Endothelium-fluid interface: %mismatch vs z (full interface, narrow domain)');
legend('Location','best'); grid on;

saveas(fig, 'mismatch_vs_z_narrow_domain.png');
fprintf('Saved mismatch_vs_z_narrow_domain.png\n');
fprintf('SMOKE_TEST_STATUS: OK\n');
