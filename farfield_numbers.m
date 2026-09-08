%% FARFIELD_NUMBERS
% Pulls actual numeric endothelium stress values (not just color) at a
% few r-points out to the true outer boundary, to check whether stress

% runs directly. No rerun of the solver needed.

clc; close all;
cd(fileparts(mfilename('fullpath')));

files = {'out_2_new_t1.mat', 'out_1D_t10_for_review.mat'};
labels = {'2D', 'Hybrid'};
rQuery = [4, 6, 8, 10, 12, 15];
rPlot = linspace(4, 15, 60);
zQuery = 0;

figure('Name', 'Far-field stress vs r', 'Position', [100 100 1000 700]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
axRR = nexttile; hold(axRR, 'on'); title(axRR, '\sigma_{rr}'); xlabel(axRR,'r [\mum]'); ylabel(axRR,'Pa'); grid(axRR,'on');
axZZ = nexttile; hold(axZZ, 'on'); title(axZZ, '\sigma_{zz}'); xlabel(axZZ,'r [\mum]'); ylabel(axZZ,'Pa'); grid(axZZ,'on');
axTT = nexttile; hold(axTT, 'on'); title(axTT, '\sigma_{\theta\theta}'); xlabel(axTT,'r [\mum]'); ylabel(axTT,'Pa'); grid(axTT,'on');
axRZ = nexttile; hold(axRZ, 'on'); title(axRZ, '\sigma_{rz}'); xlabel(axRZ,'r [\mum]'); ylabel(axRZ,'Pa'); grid(axRZ,'on');
lineStyles = {'-o', '--s'};

for fi = 1:numel(files)
    S = load(files{fi});
    fn = fieldnames(S);
    out = S.(fn{1});
    if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

    k = out.stopStep;
    st = out.stateHist{k};

    parE = out.par;
    useVisco = isfield(st, 'uEPrev') && ~isempty(st.uEPrev) && ...
        isfield(out, 'dtHist') && numel(out.dtHist) >= k && isfinite(out.dtHist(k));

    if useVisco
        dtStep = out.dtHist(k);
        stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, parE);
    else
        stressE = recover_nodal_stress_axisym(out.meshE, st.uE, parE);
    end

    Rnod = out.meshE.nodes(:,1);
    Znod = out.meshE.nodes(:,2);
    rDef = Rnod + st.uE(1:2:end);
    zDef = Znod + st.uE(2:2:end);

    fprintf('\n============ %s (%s, step %d, t=%.4e s, viscoRecovery=%d) ============\n', labels{fi}, files{fi}, k, out.t(k), useVisco);
    fprintf('%6s | %10s %10s %10s %10s\n', 'r[um]', 'sigma_rr', 'sigma_zz', 'sigma_tt', 'sigma_rz');

    Frr = scatteredInterpolant(rDef*1e6, zDef*1e6, stressE.sigma_rr, 'linear', 'nearest');
    Fzz = scatteredInterpolant(rDef*1e6, zDef*1e6, stressE.sigma_zz, 'linear', 'nearest');
    Ftt = scatteredInterpolant(rDef*1e6, zDef*1e6, stressE.sigma_tt, 'linear', 'nearest');
    Frz = scatteredInterpolant(rDef*1e6, zDef*1e6, stressE.sigma_rz, 'linear', 'nearest');

    for r = rQuery
        vrr = Frr(r, zQuery);
        vzz = Fzz(r, zQuery);
        vtt = Ftt(r, zQuery);
        vrz = Frz(r, zQuery);
        fprintf('%6.1f | %10.4f %10.4f %10.4f %10.4f\n', r, vrr, vzz, vtt, vrz);
    end

    plot(axRR, rPlot, Frr(rPlot, zQuery*ones(size(rPlot))), lineStyles{fi}, 'DisplayName', labels{fi}, 'MarkerIndices', 1:10:numel(rPlot));
    plot(axZZ, rPlot, Fzz(rPlot, zQuery*ones(size(rPlot))), lineStyles{fi}, 'DisplayName', labels{fi}, 'MarkerIndices', 1:10:numel(rPlot));
    plot(axTT, rPlot, Ftt(rPlot, zQuery*ones(size(rPlot))), lineStyles{fi}, 'DisplayName', labels{fi}, 'MarkerIndices', 1:10:numel(rPlot));
    plot(axRZ, rPlot, Frz(rPlot, zQuery*ones(size(rPlot))), lineStyles{fi}, 'DisplayName', labels{fi}, 'MarkerIndices', 1:10:numel(rPlot));
end

yline(axRR, 0, 'k:', 'HandleVisibility', 'off');
yline(axZZ, 0, 'k:', 'HandleVisibility', 'off');
yline(axTT, 0, 'k:', 'HandleVisibility', 'off');
yline(axRZ, 0, 'k:', 'HandleVisibility', 'off');
legend(axRR, 'Location', 'best', 'Interpreter', 'none');
sgtitle('Endothelium stress vs r at z=0 -- does it decay to zero by r=15?');

fprintf('\nDone. sigma_zz/sigma_rz should be near-zero from r~6-8 out. sigma_rr/sigma_tt will NOT be near-zero -- that''s the real finding.\n');
