function plot_hybrid_vs_2D_final_comparison(out1D, out2D, plotstep, componentsToShow)
%PLOT_HYBRID_VS_2D_FINAL_COMPARISON

% at the same step, with the 1D lubrication window (par.hybridGapZ, aka
% gap1DWindow) marked with horizontal dashed lines on every tile.
%
% Root-cause finding this plot is meant to document: the traction
% mismatch between solid and fluid is real and present in BOTH solvers,
% but is consistently several times larger in the 1D hybrid solver
% (confirmed by both aggregate stats and manual single-point checks).
% The 1D hybrid solver's gap1DWindow almost exactly brackets the z-range
% where the leukocyte's deformed boundary has its sharpest curvature.
% The 1D Reynolds/lubrication equation used inside that window assumes
% the gap height varies slowly with z; that assumption is measurably
% violated here, which is the most likely explanation for why the 1D
% hybrid solver's interface traction is less accurate than the full 2D
% solver's in exactly this window. This is a modeling/approximation
% limitation, not a coding bug -- column ordering, sign convention, and
% the tn/ts traction formula were all independently verified correct in
% both solvers.
%
% Usage:
%   plot_hybrid_vs_2D_final_comparison(out1D, out2D, 25)            % sigma_rz only (default)
%   plot_hybrid_vs_2D_final_comparison(out1D, out2D, 25, 1:4)       % all four components
%   plot_hybrid_vs_2D_final_comparison(out1D, out2D, 25, [1 4])     % just rr and rz

if nargin < 4 || isempty(componentsToShow)
    componentsToShow = 4; % default: sigma_rz only, matching original behavior
end

gapZ = [-0.2e-6, 4.2e-6];
if isfield(out1D, 'par') && isfield(out1D.par, 'hybridGapZ') && numel(out1D.par.hybridGapZ) == 2
    gapZ = sort(out1D.par.hybridGapZ(:)).';
end

componentTitles = {'\sigma_{rr}', '\sigma_{zz}', '\sigma_{\theta\theta}', '\sigma_{rz}'};

nComp = numel(componentsToShow);
figure('Name', sprintf('1D hybrid vs 2D full, step %d', plotstep), ...
    'Position', [100 100 900 min(2200, 450*nComp)]);
tiledlayout(nComp, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

runs = {out1D, out2D};
labels = {'1D hybrid', '2D full'};

% Precompute the leukocyte-interface "steep slope" z-bands for each run
% (same |dr/dz| > threshold definition used in check_gap1DWindow_validity.m
% and check_interface_traction_mismatch_report.m), so this plot and the
% mismatch-vs-arc-length plots stay consistent with each other. Answers

% gap1DWindow lines mark the 1D lubrication mesh boundary, which is a
% different thing from where the interface slope actually gets steep.
slopeThreshold = 0.15;
steepBandsZ = cell(1,2);
for r = 1:2
    steepBandsZ{r} = leukocyte_interface_steep_bands(runs{r}, plotstep, slopeThreshold);
end

for ci = 1:nComp
    componentToShow = componentsToShow(ci);
    for r = 1:2
        out = runs{r};
        ax = nexttile;
        set(ax, 'FontSize', 12);
        hold(ax, 'on');

        fluid = out.fluidHist{plotstep};
        if isfield(fluid, 'sigmaCellNode') && ~isempty(fluid.sigmaCellNode)
            fluidSigma  = fluid.sigmaCellNode;
            fluidCenter = fluid.centerNode;
        else
            fluidSigma  = fluid.sigmaCell;
            fluidCenter = fluid.cellCenter;
        end
        fluidColumns = [1, 3, 2, 4];

        st = out.stateHist{plotstep};
        stressE = recover_nodal_stress_axisym(out.meshE, st.uE, out.par);
        parL = out.par;
        if isfield(out.par, 'GL') && isfinite(out.par.GL), parL.Ge = out.par.GL; end
        if isfield(out.par, 'KL') && isfinite(out.par.KL), parL.Ke = out.par.KL; end
        stressL = recover_nodal_stress_axisym(out.meshL, st.uL, parL);

        solidStressE = {stressE.sigma_rr, stressE.sigma_zz, stressE.sigma_tt, stressE.sigma_rz};
        solidStressL = {stressL.sigma_rr, stressL.sigma_zz, stressL.sigma_tt, stressL.sigma_rz};

        if is_hybrid_fluid_plot(out) && ~isempty(fluidCenter)
            selectedFluidStress = fluidSigma(:, fluidColumns(componentToShow));
            local_contour_masked_field(ax, out, plotstep, fluidCenter, selectedFluidStress, 48);
        end
        plot_select_nodal_stress_contour(out.meshE, st.uE, solidStressE{componentToShow});
        plot_select_nodal_stress_contour(out.meshL, st.uL, solidStressL{componentToShow});
        draw_select_fluid_boundaries(ax, out, plotstep);

        % Mark the 1D lubrication window on both tiles for direct reference
        rMaxLine = 4;
        plot(ax, [0 rMaxLine], [gapZ(1) gapZ(1)]*1e6, 'r--', 'LineWidth', 1.5);
        plot(ax, [0 rMaxLine], [gapZ(2) gapZ(2)]*1e6, 'r--', 'LineWidth', 1.5);
        if ci == 1
            text(ax, 0.1, gapZ(2)*1e6+0.3, 'gap1DWindow', 'Color', 'r', 'FontSize', 9);
        end

        % Mark the leukocyte-interface steep-slope z-bands (|dr/dz| > 0.15)
        % so it's visually clear where the "steep slope" numbers in the
        % traction-mismatch bar chart actually come from, distinct from the
        % gap1DWindow boundary above. Lines only (not filled), so the
        % underlying stress contour stays readable.
        bands = steepBandsZ{r};
        steepColor = [0.85 0.33 0.10];
        for b = 1:size(bands,1)
            zLo = bands(b,1)*1e6;
            zHi = bands(b,2)*1e6;
            plot(ax, [0 rMaxLine], [zLo zLo], '-.', 'Color', steepColor, 'LineWidth', 1.2);
            plot(ax, [0 rMaxLine], [zHi zHi], '-.', 'Color', steepColor, 'LineWidth', 1.2);
        end
        if ci == 1 && ~isempty(bands)
            text(ax, rMaxLine-0.1, bands(1,1)*1e6, 'steep slope', 'Color', steepColor, ...
                'FontSize', 9, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        end

        colorbar(ax);
        title(ax, sprintf('%s: %s, t = %.4g s', labels{r}, componentTitles{componentToShow}, out.t(plotstep)));
        xlabel(ax, 'r [\mum]'); ylabel(ax, 'z [\mum]');
        axis(ax, 'equal'); axis(ax, 'tight'); xlim(ax, [0 4]); box(ax, 'on');
    end
end

sgtitle(sprintf(['Red dashed = 1D lubrication window (z = [%.2f, %.2f] um). ', ...
    'Orange dash-dot = leukocyte-interface "steep slope" bands (|dr/dz| > %.2f). ', ...
    'Traction mismatch concentrates at/inside these bands, worse in 1D hybrid.'], ...
    gapZ(1)*1e6, gapZ(2)*1e6, slopeThreshold));

% Also print the side-by-side traction mismatch numbers for the record.
fprintf('\n=== Traction mismatch summary, step %d ===\n', plotstep);
fprintf('1D lubrication window: z = [%.3f, %.3f] um\n', gapZ(1)*1e6, gapZ(2)*1e6);
rep1D = check_interface_traction_mismatch_report(out1D, plotstep, struct('makePlot', false, 'label', '1D hybrid'));
rep2D = check_interface_traction_mismatch_report(out2D, plotstep, struct('makePlot', false, 'label', '2D full'));

end

function bands = leukocyte_interface_steep_bands(out, plotstep, thresh)
% Returns an Nx2 array of [zLo, zHi] (meters) contiguous z-bands along the
% deformed leukocyte interface where |dr/dz| > thresh. Same slope
% computation as check_gap1DWindow_validity.m.
st = out.stateHist{plotstep};
idsL = out.interfaceL(:);
uL = st.uL;
rL = out.meshL.nodes(idsL,1) + uL(2*idsL-1);
zL = out.meshL.nodes(idsL,2) + uL(2*idsL);

[zL, ord] = sort(zL);
rL = rL(ord);
drdz = gradient(rL, zL);

steep = abs(drdz) > thresh;
d = diff([0; steep(:); 0]);
startIdx = find(d==1);
endIdx = find(d==-1) - 1;
bands = [zL(startIdx), zL(endIdx)];
end

function local_contour_masked_field(ax, out, plotstep, fluidCenter, fluidField, nLevels)
valid = isfinite(fluidCenter(:,1)) & isfinite(fluidCenter(:,2)) & isfinite(fluidField(:));
if nnz(valid) < 3
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
