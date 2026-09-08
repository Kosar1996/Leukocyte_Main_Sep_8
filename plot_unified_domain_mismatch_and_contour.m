%% PLOT_UNIFIED_DOMAIN_MISMATCH_AND_CONTOUR
% Visual companion to regenerate_slide4_unified.m: (1) %mismatch vs z
% across the full interface, both sides, unified-domain final state --
% same style as the original report's full-interface plot, so it's
% directly comparable; (2) fluid pressure contour over the unified
% domain, so the >80% mismatch finding can be seen spatially, not just
% as point numbers in a table.

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

mismatchOpts = struct('nQuery', 121, 'epsFrac', 0.1, 'trimFrac', 0.02);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);

pctEnL = 100*abs(leuko.en)/leuko.stats.refScale;
pctEtL = 100*abs(leuko.et)/leuko.stats.refScale;
pctEnE = 100*abs(endo.en)/endo.stats.refScale;
pctEtE = 100*abs(endo.et)/endo.stats.refScale;

zPoints = [3.45, 3.50];

fig1 = figure('Position',[100 100 900 700]);
subplot(2,1,1);
plot(leuko.z*1e6, pctEnL, 'b-o','MarkerSize',3); hold on;
plot(leuko.z*1e6, pctEtL, 'r-s','MarkerSize',3);
[~,i1]=min(abs(leuko.z*1e6-zPoints(1))); [~,i2]=min(abs(leuko.z*1e6-zPoints(2)));
plot(leuko.z(i1)*1e6, pctEnL(i1), 'k*','MarkerSize',14,'LineWidth',2);
plot(leuko.z(i2)*1e6, pctEnL(i2), 'k*','MarkerSize',14,'LineWidth',2);
xlabel('z [\mum]'); ylabel('% mismatch');
title('Leukocyte-fluid interface: %mismatch vs z (unified domain, step 40)');
legend('Normal (En)','Tangential (Et)','Report points (L1,L2)','Location','best');
grid on;

subplot(2,1,2);
plot(endo.z*1e6, pctEnE, 'b-o','MarkerSize',3); hold on;
plot(endo.z*1e6, pctEtE, 'r-s','MarkerSize',3);
[~,i1]=min(abs(endo.z*1e6-zPoints(1))); [~,i2]=min(abs(endo.z*1e6-zPoints(2)));
plot(endo.z(i1)*1e6, pctEnE(i1), 'k*','MarkerSize',14,'LineWidth',2);
plot(endo.z(i2)*1e6, pctEnE(i2), 'k*','MarkerSize',14,'LineWidth',2);
xlabel('z [\mum]'); ylabel('% mismatch');
title('Endothelium-fluid interface: %mismatch vs z (unified domain, step 40)');
legend('Normal (En)','Tangential (Et)','Report points (E1,E2)','Location','best');
grid on;

saveas(fig1, 'unified_domain_mismatch_vs_z.png');
fprintf('Saved unified_domain_mismatch_vs_z.png\n');

%% Fluid pressure contour over the unified domain
mraw = fl.meshF;
fig2 = figure('Position',[100 100 900 500]);
pgrid = reshape(fl.pCell, mraw.Nr, mraw.Nz);
pcolor(mraw.Zp*1e6, mraw.Rp*1e6, pgrid); shading interp; colorbar;
xlabel('z [\mum]'); ylabel('r [\mum]');
title('Fluid pressure contour, unified domain (z=[-2,6]), step 40');
hold on;
zLine = linspace(min(mraw.zc), max(mraw.zc), 300);
deltaLofz = @(z) interp1(mraw.zc, mraw.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mraw.zc, mraw.deltaE_c, z, 'linear', 'extrap');
plot(zLine*1e6, deltaLofz(zLine)*1e6, 'w-', 'LineWidth', 1.5);
plot(zLine*1e6, deltaEofz(zLine)*1e6, 'w-', 'LineWidth', 1.5);

saveas(fig2, 'unified_domain_pressure_contour.png');
fprintf('Saved unified_domain_pressure_contour.png\n');

fprintf('\nSMOKE_TEST_STATUS: OK\n');
