%% ENDOTHELIUM_NORMAL_VS_TANGENTIAL_OVER_TIME
% Same prestress-domination hypothesis that explains the far-field finding
% (sigma_rr/sigma_tt persist with distance, sigma_zz/sigma_rz decay), but
% checked over TIME instead of over r: does the endothelium's own NORMAL
% stress at the interface (sigma_rr, what dominates the traction-correction
% "EnE" mismatch) stay nearly frozen across real accepted steps, while the
% TANGENTIAL/shear stress (sigma_rz, what dominates "EtE") visibly moves?
%
% Uses the already-saved 10-step hybrid run -- no solver rerun needed.

clc; close all;
cd(fileparts(mfilename('fullpath')));

file = 'out_1D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

nSteps = out.stopStep;
parE = out.par;
nZ = 15;

sigRRp90 = nan(nSteps,1);
sigRZp90 = nan(nSteps,1);

for k = 1:nSteps
    st = out.stateHist{k};
    useVisco = isfield(st, 'uEPrev') && ~isempty(st.uEPrev) && ...
        isfield(out, 'dtHist') && numel(out.dtHist) >= k && isfinite(out.dtHist(k));
    if useVisco
        stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, out.dtHist(k), parE);
    else
        stressE = recover_nodal_stress_axisym(out.meshE, st.uE, parE);
    end

    Rnod = out.meshE.nodes(:,1);
    Znod = out.meshE.nodes(:,2);
    rDef = Rnod + st.uE(1:2:end);
    zDef = Znod + st.uE(2:2:end);

    deltaE_k = out.deltaEHist(:,k);   % wall radius per z at this step
    zGrid = out.z(:);
    zQuery = linspace(min(zGrid), max(zGrid), nZ)';
    rWallQuery = interp1(zGrid, deltaE_k, zQuery, 'linear', 'extrap');

    Frr = scatteredInterpolant(rDef, zDef, stressE.sigma_rr, 'linear', 'nearest');
    Frz = scatteredInterpolant(rDef, zDef, stressE.sigma_rz, 'linear', 'nearest');

    rrVals = abs(Frr(rWallQuery, zQuery));
    rzVals = abs(Frz(rWallQuery, zQuery));

    sigRRp90(k) = prctile(rrVals, 90);
    sigRZp90(k) = prctile(rzVals, 90);
end

steps = (1:nSteps)';
pctRR = 100*(sigRRp90 - sigRRp90(1)) / max(abs(sigRRp90(1)), 1e-9);
pctRZ = 100*(sigRZp90 - sigRZp90(1)) / max(abs(sigRZp90(1)), 1e-9);

fprintf('%5s | %12s %12s | %9s %9s\n', 'step', 'p90|sigRR|', 'p90|sigRZ|', '%%chgRR', '%%chgRZ');
for k = 1:nSteps
    fprintf('%5d | %12.4f %12.4f | %9.2f %9.2f\n', k, sigRRp90(k), sigRZp90(k), pctRR(k), pctRZ(k));
end

fig = figure('Position',[100 100 900 400],'Color','w');
tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile(tl);
plot(ax1, steps, sigRRp90, 'o-', 'LineWidth', 1.5, 'DisplayName', '\sigma_{rr} (normal)'); hold(ax1, 'on');
plot(ax1, steps, sigRZp90, 's-', 'LineWidth', 1.5, 'DisplayName', '\sigma_{rz} (shear)');
xlabel(ax1, 'accepted step'); ylabel(ax1, 'p90 |\sigma| at interface [Pa]');
title(ax1, 'Raw magnitude vs step'); legend(ax1, 'Location', 'best'); grid(ax1, 'on');

ax2 = nexttile(tl);
plot(ax2, steps, pctRR, 'o-', 'LineWidth', 1.5, 'DisplayName', '\sigma_{rr} (normal)'); hold(ax2, 'on');
plot(ax2, steps, pctRZ, 's-', 'LineWidth', 1.5, 'DisplayName', '\sigma_{rz} (shear)');
yline(ax2, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax2, 'accepted step'); ylabel(ax2, '%% change from step 1');
title(ax2, 'Relative change vs step'); legend(ax2, 'Location', 'best'); grid(ax2, 'on');

sgtitle(fig, sprintf('%s: endothelium interface stress over time', file), 'Interpreter', 'none');

outPng = fullfile(fileparts(mfilename('fullpath')), 'endothelium_normal_vs_tangential_over_time.png');
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('\nSaved: %s\n', outPng);
fprintf('If sigma_rr stays much flatter (in %% terms) than sigma_rz across these 10 real steps,\n');
fprintf('that is independent support for the same prestress-domination mechanism seen in EnE.\n');
