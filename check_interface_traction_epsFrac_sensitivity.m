function result = check_interface_traction_epsFrac_sensitivity(out, plotstep, opts)
%CHECK_INTERFACE_TRACTION_EPSFRAC_SENSITIVITY

% mismatch come from the force or area? dig into it more."
%
% check_interface_traction_mismatch_report.m samples the fluid-side and
% solid-side stress at a query point offset from the wall by
% opts.epsFrac * localGapWidth on each side (not exactly at the shared
% interface). If the reported mismatch were a genuine force discontinuity
% at the interface, it should not depend much on that offset. If it's
% mostly a sampling-distance artifact (each side's own stress field has a
% real spatial gradient near the wall, so comparing points further apart
% compares increasingly different physical locations), the reported
% mismatch should grow with epsFrac.
%
% This function reruns check_interface_traction_mismatch_report.m across
% a sweep of epsFrac values, fits max|en|/max|et| (each side) linearly vs
% epsFrac, and extrapolates to epsFrac=0 (the true wall) to split the
% reported mismatch at any given epsFrac into:
%   - a "genuine" part (the eps=0 intercept -- does not vanish as the
%     sampling offset shrinks to zero, so it's not explained by the
%     sampling methodology)
%   - an "offset-artifact" part (reported value minus the intercept)
%
% IMPORTANT CAVEAT (kept in the output, not just this comment): eps=0 is
% never actually sampled -- it is an extrapolation of the fitted line, not
% a direct measurement, because sampling exactly at the shared boundary is
% ill-defined for this interpolation scheme. The R^2 of the linear fit is
% returned alongside every extrapolated value so this can be judged
% instead of asserted.
%
% Usage:
%   result = check_interface_traction_epsFrac_sensitivity(out, plotstep)
%   result = check_interface_traction_epsFrac_sensitivity(out, plotstep, opts)
%
% opts (all optional):
%   opts.epsFracList  vector of epsFrac values to sweep
%                     (default [0.005 0.01 0.02 0.05 0.1 0.2 0.3])
%   opts.refEpsFrac   the epsFrac value that reported numbers are usually
%                     quoted at, used to compute the artifact fraction
%                     (default 0.1, matching check_interface_traction_
%                     mismatch_report.m's own default)
%   opts.label        string used in plot titles (default '')
%   opts.makePlot     true/false (default true)
%
% result fields (each of leuko_en/leuko_et/endo_en/endo_et is a struct):
%   .epsFracList          the swept epsFrac values
%   .<side>_<comp>.values         max|.| at each epsFrac
%   .<side>_<comp>.slope          Pa per unit epsFrac
%   .<side>_<comp>.intercept      extrapolated value at epsFrac=0 [Pa]
%   .<side>_<comp>.R2             linear fit quality, 0-1
%   .<side>_<comp>.reportedAtRef  value at the epsFrac closest to opts.refEpsFrac
%   .<side>_<comp>.artifactFrac   (reportedAtRef - intercept) / reportedAtRef

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts,'epsFracList') || isempty(opts.epsFracList)
    opts.epsFracList = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3];
end
if ~isfield(opts,'refEpsFrac') || isempty(opts.refEpsFrac)
    opts.refEpsFrac = 0.1;
end
if ~isfield(opts,'label'), opts.label = ''; end
if ~isfield(opts,'makePlot') || isempty(opts.makePlot), opts.makePlot = true; end

epsFracList = opts.epsFracList(:);
nE = numel(epsFracList);

leuko_en = zeros(nE,1); leuko_et = zeros(nE,1);
endo_en  = zeros(nE,1); endo_et  = zeros(nE,1);

for i = 1:nE
    rep = check_interface_traction_mismatch_report(out, plotstep, ...
        struct('makePlot', false, 'epsFrac', epsFracList(i), 'label', ''));
    leuko_en(i) = rep.leuko.stats.maxAbsEn;
    leuko_et(i) = rep.leuko.stats.maxAbsEt;
    endo_en(i)  = rep.endo.stats.maxAbsEn;
    endo_et(i)  = rep.endo.stats.maxAbsEt;
end

result = struct();
result.epsFracList = epsFracList;
result.leuko_en = fit_and_pack(epsFracList, leuko_en, opts.refEpsFrac);
result.leuko_et = fit_and_pack(epsFracList, leuko_et, opts.refEpsFrac);
result.endo_en  = fit_and_pack(epsFracList, endo_en,  opts.refEpsFrac);
result.endo_et  = fit_and_pack(epsFracList, endo_et,  opts.refEpsFrac);

fprintf('\n=== Interface traction mismatch: epsFrac sensitivity, %s step %d ===\n', ...
    opts.label, plotstep);
fprintf(['  (extrapolated eps=0 value is a FIT-BASED ESTIMATE, not a direct\n', ...
    '   measurement -- see .R2 for how well the linear trend is supported\n', ...
    '   over the tested epsFrac range)\n\n']);
print_row('Leukocyte |en|', result.leuko_en, opts.refEpsFrac);
print_row('Leukocyte |et|', result.leuko_et, opts.refEpsFrac);
print_row('Endothelium |en|', result.endo_en, opts.refEpsFrac);
print_row('Endothelium |et|', result.endo_et, opts.refEpsFrac);

if opts.makePlot
    figure('Name', sprintf('Traction mismatch epsFrac sensitivity %s step %d', opts.label, plotstep));
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    plot_panel(epsFracList, leuko_en, result.leuko_en, 'Leukocyte |e_n|');
    plot_panel(epsFracList, leuko_et, result.leuko_et, 'Leukocyte |e_t|');
    plot_panel(epsFracList, endo_en,  result.endo_en,  'Endothelium |e_n|');
    plot_panel(epsFracList, endo_et,  result.endo_et,  'Endothelium |e_t|');
    sgtitle(sprintf(['%s step %d: max interface traction mismatch vs sampling offset epsFrac.\n', ...
        'Dashed line = linear fit; circle at epsFrac=0 = extrapolated wall value (estimate, not measured).'], ...
        opts.label, plotstep));
end

end

function s = fit_and_pack(epsFracList, values, refEpsFrac)
p = polyfit(epsFracList, values, 1);
fitted = polyval(p, epsFracList);
ssRes = sum((values - fitted).^2);
ssTot = sum((values - mean(values)).^2);
if ssTot > 0
    R2 = 1 - ssRes/ssTot;
else
    R2 = NaN;
end

[~, iRef] = min(abs(epsFracList - refEpsFrac));
reportedAtRef = values(iRef);
intercept = p(2);

s = struct();
s.values = values;
s.slope = p(1);
s.intercept = intercept;
s.R2 = R2;
s.reportedAtRef = reportedAtRef;
if reportedAtRef ~= 0
    s.artifactFrac = (reportedAtRef - intercept) / reportedAtRef;
else
    s.artifactFrac = NaN;
end
end

function print_row(name, s, refEpsFrac)
fprintf('  %-18s reported@eps=%.2g: %8.3f Pa | extrapolated eps=0: %8.3f Pa | R^2=%.4f | %.1f%% offset-dependent\n', ...
    name, refEpsFrac, s.reportedAtRef, s.intercept, s.R2, 100*s.artifactFrac);
end

function plot_panel(epsFracList, values, s, titleStr)
nexttile;
plot(epsFracList, values, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.2 0.4 0.8]); hold on;
xFit = linspace(0, max(epsFracList), 50);
yFit = s.slope*xFit + s.intercept;
plot(xFit, yFit, 'k--', 'LineWidth', 1.2);
plot(0, s.intercept, 'ro', 'MarkerSize', 9, 'LineWidth', 1.5);
text(0, s.intercept, sprintf('  %.1f Pa (R^2=%.3f)', s.intercept, s.R2), ...
    'Color', 'r', 'FontSize', 9, 'VerticalAlignment', 'bottom');
xlabel('epsFrac (fraction of local gap width)');
ylabel('max |mismatch| [Pa]');
title(sprintf('%s (%.0f%% offset-dependent @ eps=0.1)', titleStr, 100*s.artifactFrac));
grid on;
end
