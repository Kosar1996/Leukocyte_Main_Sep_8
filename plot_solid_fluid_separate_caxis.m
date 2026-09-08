function plot_solid_fluid_separate_caxis(out, plotstep)
%PLOT_SOLID_FLUID_SEPARATE_CAXIS

% with its own independent (auto) color scale, instead of one shared
% colorbar. In the combined plot_select_native2d_stress.m figure the
% solid values (much larger magnitude) can visually dominate and hide
% fluid-side structure, since both are drawn on the same caxis range.
%
% Produces two 2x2 figures for the given step:
%   Figure 1: SOLID  sigma_rr, sigma_zz, sigma_tt, sigma_rz  (own caxis each)
%   Figure 2: FLUID   sigma_rr, sigma_zz, sigma_tt, sigma_rz  (own caxis each)
%
% Usage:
%   plot_solid_fluid_separate_caxis(out, plotstep)

fluid = out.fluidHist{plotstep};
if isfield(fluid, 'sigmaCellNode') && ~isempty(fluid.sigmaCellNode)
    fluidSigma  = fluid.sigmaCellNode;
    fluidCenter = fluid.centerNode;
elseif isfield(fluid, 'sigmaCell') && ~isempty(fluid.sigmaCell)
    fluidSigma  = fluid.sigmaCell;
    fluidCenter = fluid.cellCenter;
else
    warning('No recovered fluid stress field found.');
    fluidSigma = []; fluidCenter = [];
end
% storage order [srr, stt, szz, srz] -> display order [rr, zz, tt, rz]
fluidColumns = [1, 3, 2, 4];

st = out.stateHist{plotstep};
stressE = recover_nodal_stress_axisym(out.meshE, st.uE, out.par);

parL = out.par;
if isfield(out.par, 'GL') && isfinite(out.par.GL), parL.Ge = out.par.GL; end
if isfield(out.par, 'KL') && isfinite(out.par.KL), parL.Ke = out.par.KL; end
stressL = recover_nodal_stress_axisym(out.meshL, st.uL, parL);

componentTitles = {'\sigma_{rr}', '\sigma_{zz}', '\sigma_{\theta\theta}', '\sigma_{rz}'};
solidStressE = {stressE.sigma_rr, stressE.sigma_zz, stressE.sigma_tt, stressE.sigma_rz};
solidStressL = {stressL.sigma_rr, stressL.sigma_zz, stressL.sigma_tt, stressL.sigma_rz};

% ---------------- SOLID-ONLY figure ----------------
figure('Name', sprintf('SOLID stress only, step %d', plotstep));
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:4
    ax = nexttile;
    set(ax, 'FontSize', 14);
    hold(ax, 'on');
    plot_select_nodal_stress_contour(out.meshE, st.uE, solidStressE{k});
    plot_select_nodal_stress_contour(out.meshL, st.uL, solidStressL{k});
    draw_select_fluid_boundaries(ax, out, plotstep);
    colorbar(ax);  % independent auto caxis per tile (default behavior)
    title(ax, sprintf('SOLID %s, t = %.4g s', componentTitles{k}, out.t(plotstep)));
    xlabel(ax, 'r [\mum]'); ylabel(ax, 'z [\mum]');
    axis(ax, 'equal'); axis(ax, 'tight'); xlim(ax, [0 4]); box(ax, 'on');
end
sgtitle('Solid-only elastic stress, independent color scale per tile');

% ---------------- FLUID-ONLY figure ----------------
figure('Name', sprintf('FLUID stress only, step %d', plotstep));
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:4
    ax = nexttile;
    set(ax, 'FontSize', 14);
    hold(ax, 'on');
    if is_hybrid_fluid_plot(out) && ~isempty(fluidCenter)
        selectedFluidStress = fluidSigma(:, fluidColumns(k));
        local_contour_masked_field(ax, out, plotstep, fluidCenter, selectedFluidStress, 48);
    end
    draw_select_fluid_boundaries(ax, out, plotstep);
    colorbar(ax);
    title(ax, sprintf('FLUID %s, t = %.4g s', componentTitles{k}, out.t(plotstep)));
    xlabel(ax, 'r [\mum]'); ylabel(ax, 'z [\mum]');
    axis(ax, 'equal'); axis(ax, 'tight'); xlim(ax, [0 4]); box(ax, 'on');
end
sgtitle('Fluid-only stress (pressure+viscous), independent color scale per tile');

end

function local_contour_masked_field(ax, out, plotstep, fluidCenter, fluidField, nLevels)
% Same regular-grid + solid-mask approach as
% plot_select_native2d_stress.m's contour_masked_regular_field_select,
% duplicated here (as a local function) so this file can run standalone.
valid = isfinite(fluidCenter(:,1)) & isfinite(fluidCenter(:,2)) & isfinite(fluidField(:));
if nnz(valid) < 3
    warning('Not enough finite fluid field points to plot.');
    return;
end
rVals = fluidCenter(valid,1);
zVals = fluidCenter(valid,2);
fVals = fluidField(valid);

nrPlot = 260; nzPlot = 260;
rGrid = linspace(min(rVals), max(rVals), nrPlot);
zGrid = linspace(min(zVals), max(zVals), nzPlot);
[RGrid, ZGrid] = meshgrid(rGrid, zGrid);

F = scatteredInterpolant(rVals(:), zVals(:), fVals(:), 'linear', 'none');
FGrid = F(RGrid, ZGrid);

statePlot = state_for_plot_at_step(out, plotstep);
solidMask = false(size(RGrid));
if isfield(out, 'meshL') && isfield(statePlot, 'uL') && ~isempty(out.meshL) && ~isempty(statePlot.uL)
    solidMask = solidMask | local_mask(out.meshL, statePlot.uL, RGrid, ZGrid);
end
if isfield(out, 'meshE') && isfield(statePlot, 'uE') && ~isempty(out.meshE) && ~isempty(statePlot.uE)
    solidMask = solidMask | local_mask(out.meshE, statePlot.uE, RGrid, ZGrid);
end
FGrid(solidMask) = NaN;

finiteVals = FGrid(isfinite(FGrid));
if isempty(finiteVals)
    warning('No finite fluid field values remain after solid masking.');
    return;
end
if max(finiteVals) > min(finiteVals)
    contourf(ax, RGrid*1e6, ZGrid*1e6, FGrid, nLevels, 'LineColor', 'none');
else
    contourf(ax, RGrid*1e6, ZGrid*1e6, FGrid, 1, 'LineColor', 'none');
end
end

function mask = local_mask(mesh, u, RGrid, ZGrid)
mask = false(size(RGrid));
if isempty(mesh) || isempty(u) || ~isfield(mesh,'nodes') || ~isfield(mesh,'conn') || ...
        numel(u) < 2*size(mesh.nodes,1)
    return;
end
rDef = mesh.nodes(:,1) + u(1:2:end);
zDef = mesh.nodes(:,2) + u(2:2:end);
for e = 1:size(mesh.conn,1)
    ids = mesh.conn(e,:);
    rv = rDef(ids); zv = zDef(ids);
    inBox = RGrid >= min(rv) & RGrid <= max(rv) & ZGrid >= min(zv) & ZGrid <= max(zv);
    if any(inBox(:))
        localMask = false(size(mask));
        localMask(inBox) = inpolygon(RGrid(inBox), ZGrid(inBox), rv, zv);
        mask = mask | localMask;
    end
end
end
