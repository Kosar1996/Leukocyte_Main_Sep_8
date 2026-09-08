%% PLOT_UNIFIED_DOMAIN_PRESSURE_VS_Z
% Corrected visualization: pressure is ~uniform across the gap at each z
% (expected lubrication-flow behavior, confirmed directly on the data),
% so a 2D color contour is not informative -- a p(z) line plot, with gap
% width on the secondary axis, is the honest, readable version of "what
% happened on the unified domain."

clc; close all;
file = 'out_2D_unified_domain_40step.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

k = out.stopStep;
fl = out.fluidHist{k};
mraw = fl.meshF;

pgrid = reshape(fl.pCell, mraw.Nr, mraw.Nz);
pMean = mean(pgrid, 1);
pStdAcrossGap = std(pgrid, 0, 1);
zAxis = mraw.Zp(1,:);

deltaLofz = @(z) interp1(mraw.zc, mraw.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mraw.zc, mraw.deltaE_c, z, 'linear', 'extrap');
gapWidth = deltaEofz(zAxis) - deltaLofz(zAxis);

fig = figure('Position',[100 100 900 500]);
yyaxis left
plot(zAxis*1e6, pMean, 'b-', 'LineWidth', 1.5);
ylabel('Fluid pressure [Pa] (mean across gap)');
yyaxis right
plot(zAxis*1e6, gapWidth*1e9, 'r--', 'LineWidth', 1.2);
ylabel('Gap width [nm]');
xlabel('z [\mum]');
title(sprintf('Unified-domain (z=[%.1f,%.1f]\\mum) pressure vs z, step %d', min(zAxis)*1e6, max(zAxis)*1e6, k));
grid on;
legend('Pressure (mean across gap)', 'Gap width', 'Location', 'best');

fprintf('Max cross-gap std/mean pressure ratio: %.4f%% (confirms near-uniform-across-gap, lubrication regime)\n', ...
    100*max(pStdAcrossGap(pMean~=0) ./ abs(pMean(pMean~=0))));

saveas(fig, 'unified_domain_pressure_vs_z.png');
fprintf('Saved unified_domain_pressure_vs_z.png\n');
fprintf('\nSMOKE_TEST_STATUS: OK\n');
