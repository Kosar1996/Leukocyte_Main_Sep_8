% %% RUN_INPUT_FULL2D_PRESSURE
% % Input deck for the hybrid pressure calculation:
% %   1D lubrication mesh in the leukocyte/endothelium gap,
% %   2D body-fitted MAC/Stokes mesh in the upstream/downstream exterior.
% %
% % This file keeps run_input.m unchanged and starts a separate case that uses
% % the mixed-dimensional fluid solve inside:
% %   softlube_run_case_global_coupled(cfg)
% %
% % Important distinction:
% %   - The gap interval uses the Reynolds/lubrication p(z) solve.
% %   - The exterior intervals are solved as true 2D body-fitted MAC/Stokes
% %     blocks, matched to the 1D pressure at z = -0.2 and z = 4.2 um.
% %   - The final composite pressure field is stored as P(r,z) in out.PHist,
% %     on the physical grid out.RPHist/out.ZPHist.
% %
% % To run:
% %   cd('/Users/shu/Softlubrication')
% %   run_input_full2D_pressure
% %
% % Main hybrid outputs:
% %   out.PHist(:,:,end)       composite hybrid pressure P(r,z) at final step
% %   out.RPHist(:,:,end)      physical r locations for P
% %   out.ZPHist(:,:,end)      physical z locations for P
% %   out.urCHist(:,:,end)     radial velocity at pressure-cell centers
% %   out.uzCHist(:,:,end)     axial velocity at pressure-cell centers
% %   out.native2D             compact final-step hybrid 2D fields
% %   out.hybridMesh           final 1D-gap / 2D-exterior mesh description
% 
% clc;
% clearvars;
% 
% %% 1. Paths
% softlubeDir = '';
% addpath(softlubeDir);
% 
% %% 2. Time Controls
% dt = 2e-4;           % time step size [s]
% nSteps = 60;          % requested time steps
% tEnd = nSteps * dt;
% 
% useAdaptiveTimeStep = true;
% dtMin = dt / 1024;
% 
% %% 3. Global Domain and Pressure Boundary Conditions
% % The hybrid pressure uses a 1D gap solve and 2D exterior body-fitted
% % r-z blocks. The axial ends use P = 0. In the exterior MAC/Stokes blocks,
% % r = 0 is the axis outside the leukocyte and r = 15 um is the exterior wall
% % outside the endothelium.
% zMin = -6e-6;        % -6 um
% zMax = 10e-6;        % 10 um
% rOuter = 15e-6;      % 15 um
% 
% pAtZMin = 0;         % P = 0 at z = -6 um [Pa]
% pAtZMax = 0;         % P = 0 at z = 10 um [Pa]
% axisPressureBC = 'dPdr=0';
% 
% %% 4. Fluid Mesh Resolution
% % Hybrid split. The gap interval is handled by the 1D lubrication mesh; the
% % exterior intervals are handled by the 2D body-fitted MAC/Stokes mesh.
% gap1DWindow = [-0.2e-6, 4.2e-6];
% 
% % Resolve the 1D lubrication region finely and keep the exterior coarse.
% finePressureWindow = gap1DWindow;
% coarseDz = 1.0e-6;
% fineDz = 0.05e-6;
% 
% % 2D MAC cells across the radial exterior fluid blocks. The radial grid is
% % nonuniform below, with extra resolution in the thin-gap band r = 2..4 um.
% NrExterior2D = 20;
% NrFluid2D = NrExterior2D;
% nrEnv = str2double(getenv('SOFTLUBE_FULL2D_NR'));
% if isfinite(nrEnv) && nrEnv >= 5
%     NrFluid2D = round(nrEnv);
%     NrExterior2D = NrFluid2D;
% end
% radialFineWindow = [2e-6, 4e-6];
% radialFineWeight = 4.0;
% 
% %% 5. Plotting and Output
% plotNative2DPressure = true;
% plotNative2DVelocity = false;
% plotFirstStepHybridMesh = true;
% plotGlobalDomainSchematic = true;
% 
% % Keep the old solver plot block off. The local plot at the bottom shows
% % only the composite hybrid result from out.PHist.
% makeLegacyPlots = false;
% 
% saveOutput = false;
% outputFile = fullfile(softlubeDir, 'simulation_output_full2D_pressure.mat');
% 
% closeFiguresAtStart = true;
% printEvery = 1;
% 
% %% 6. Build the Base Case
% cfg = struct();
% 
% % Prestressed initial geometries. Leukocyte geometry/mesh settings are read
% % inside the solver from solid_leu_P600.mat and are not repeated here.
% cfg.geometry.endotheliumPrestressFile = fullfile(softlubeDir, 'solid_endo_P300.mat');
% cfg.geometry.leukocytePrestressFile = fullfile(softlubeDir, 'solid_leuko_P6002.mat');
% cfg.geometry.usePrestressedLeukocyteIC = true;
% cfg.geometry.useLeukocyteCenterline = true;
% 
% cfg.fluid.zMin = zMin;
% cfg.fluid.zMax = zMax;
% cfg.fluid.mu = 1.2e-3;
% cfg.fluid.slipE = 0e-9;
% cfg.fluid.slipL = 0e-9;
% cfg.fluid.NrFluid2D = NrFluid2D;
% cfg.fluid.Nr = NrFluid2D;
% cfg.fluid.NrExterior2D = NrExterior2D;
% cfg.fluid.bodyFittedRadialFineWindow = radialFineWindow;
% cfg.fluid.bodyFittedRadialFineWeight = radialFineWeight;
% cfg.fluid.innerBoundary = 'deformed_leukocyte';
% cfg.fluid.useSlopeAwareWallKinematics = true;
% 
% cfg.bc.pIn = pAtZMin;
% cfg.bc.pOut = pAtZMax;
% cfg.bc.UwL = 0;
% cfg.bc.UwE = 0;
% 
% cfg.solid.noLeukocyte = false;
% cfg.solid.leukocyte.enabled = true;
% cfg.solid.leukocyte.deformable = true;
% cfg.solid.leukocyte.EL = 200;
% cfg.solid.leukocyte.nuL = 0.46;
% cfg.solid.leukocyte.useViscoelasticLeukocyte = true;
% cfg.solid.leukocyte.etaL = 1;
% cfg.solid.leukocyte.supportL = 'axis';
% 
% cfg.solid.endothelium.Ee = 500;
% cfg.solid.endothelium.nuE = 0.46;
% cfg.solid.endothelium.useViscoelasticEndothelium = true;
% cfg.solid.endothelium.etaE = 1;
% 
% cfg.interface.useExactDeformedInterface = true;
% cfg.interface.useExactDeformedInterfaceInMonolithic = true;
% cfg.interface.useReferenceZForTractionMapping = false;
% cfg.interface.zeroLeukocyteTractionOutsideOverlap = false;
% cfg.interface.leukocytePressureSupportZ = [zMin, zMax];
% cfg.interface.warnPrestressLoadMismatch = false;
% 
% cfg.numerics.dt = dt;
% cfg.numerics.tEnd = tEnd;
% cfg.numerics.minGap = 0.001e-6;
% cfg.numerics.maxNewtonMono = 100;
% cfg.numerics.tolNewtonMono = 1e-6;
% cfg.numerics.monoSolidAbsTol = 1e-7;
% cfg.numerics.monoFluidAbsTol = 1e-8;
% cfg.numerics.maxFunctionEvaluationsTwoSolid = 400;
% cfg.numerics.enablePressureJumpRetry = false;
% cfg.numerics.maxPressureJumpAbs = inf;
% cfg.numerics.maxPressureJump2DAbs = inf;
% cfg.numerics.useObjectiveKelvinVoigt = true;
% cfg.numerics.useDirectFluidSolve = true;
% cfg.numerics.useFsolveMono = true;
% cfg.numerics.enableAdaptiveTimeStep = useAdaptiveTimeStep;
% cfg.numerics.dtMin = dtMin;
% cfg.numerics.dtRetryFactor = 0.5;
% cfg.numerics.dtGrowFactor = 1.25;
% cfg.numerics.maxTimeStepRetries = 18;
% cfg.numerics.fsolveDisplay = 'off';
% cfg.numerics.diagnosticsEnabled = true;
% 
% cfg.output.saveOutput = saveOutput;
% cfg.output.outputFile = outputFile;
% cfg.output.makePlots = makeLegacyPlots;
% cfg.output.plotFinalGlobalPressureContour = false;
% cfg.output.plotGlobal2DPressureTractionComparison = false;
% cfg.output.printEvery = printEvery;
% cfg.output.storeFull2DFluidHist = false;
% cfg.output.store2DFluidArrays = true;
% cfg.output.store2DPhysicalGridHist = true;
% cfg.output.store2DSpeedHist = true;
% 
% cfg.runtime.useEnvironment = false;
% cfg.ui.closeFigures = closeFiguresAtStart;
% cfg.parOverrides = struct();
% 
% %% 7. Global Axial Mesh and Exterior Boundary Embedding
% % The flag name is historical. Here it is kept on to reuse the global axial
% % mesh and the outside-solid radius rules:
% %   leukocyte absent: fluid inner boundary becomes r = 0
% %   endothelium absent: fluid outer boundary becomes r = 15 um
% % The composite hybrid pressure output is out.PHist, not out.global1D.
% cfg.global1D.zDomain = [zMin, zMax];
% cfg.global1D.rDomain = [0, rOuter];
% cfg.global1D.axisPressureBC = axisPressureBC;
% cfg.global1D.outerPressureBC = 0;
% 
% cfg.parOverrides.useGlobal1DPressure = true;
% cfg.parOverrides.global1DOuterRadius = rOuter;
% cfg.parOverrides.global1DAxisRadius = 1e-9;
% cfg.parOverrides.useGlobal1DCoarseEdgeMesh = true;
% cfg.parOverrides.global1DFineWindow = finePressureWindow;
% cfg.parOverrides.global1DCoarseDz = coarseDz;
% cfg.parOverrides.global1DFineDz = fineDz;
% cfg.parOverrides.global1DPlotNr = 161;
% 
% %% 8. Hybrid Gap-1D / Exterior-2D Fluid Module
% cfg.fluid.useFull2DFluid = true;
% cfg.fluid.useHybridGap1DExterior2DFluid = true;
% cfg.fluid.hybridGapZ = gap1DWindow;
% cfg.fluid.hybridExteriorMinCells = 4;
% cfg.fluid.hybridExteriorFailMode = 'warn';
% cfg.fluid.hybridUseExterior2DPressureVector = true;
% cfg.fluid.useBodyFittedMACFluid = true;
% cfg.fluid.useBodyFittedMACTractionCorrection = true;
% cfg.fluid.maxBodyFittedTractionCorrections = 1;
% cfg.fluid.bodyFittedTractionCorrectionRelax = 0.05;
% cfg.fluid.bodyFittedTractionCorrectionFailMode = 'warn';
% cfg.fluid.resolveFluidAfterTractionCorrection = true;
% cfg.fluid.useBodyFittedMACTractionInSolid = true;
% 
% cfg.parOverrides.useFull2DFluid = true;
% cfg.parOverrides.useHybridGap1DExterior2DFluid = true;
% cfg.parOverrides.hybridGapZ = gap1DWindow;
% cfg.parOverrides.NrExterior2D = NrExterior2D;
% cfg.parOverrides.bodyFittedRadialFineWindow = radialFineWindow;
% cfg.parOverrides.bodyFittedRadialFineWeight = radialFineWeight;
% cfg.parOverrides.hybridExteriorMinCells = 4;
% cfg.parOverrides.hybridExteriorFailMode = 'warn';
% cfg.parOverrides.hybridUseExterior2DPressureVector = true;
% cfg.parOverrides.useBodyFittedMACFluid = true;
% cfg.parOverrides.useBodyFittedMACTractionCorrection = true;
% cfg.parOverrides.maxBodyFittedTractionCorrections = 1;
% cfg.parOverrides.bodyFittedTractionCorrectionRelax = 0.05;
% cfg.parOverrides.bodyFittedTractionCorrectionFailMode = 'warn';
% cfg.parOverrides.resolveFluidAfterTractionCorrection = true;
% cfg.parOverrides.useBodyFittedMACTractionInSolid = true;
% cfg.parOverrides.NrFluid2D = NrFluid2D;
% cfg.parOverrides.Nr = NrFluid2D;
% 
% % Turn off the earlier projected-pressure bridge so this file reports the
% % hybrid pressure field cleanly.
% cfg.parOverrides.useGlobal2DPressureTraction = false;
% 
% %% 9. Runtime
% if plotNative2DPressure || plotNative2DVelocity || ...
%         plotFirstStepHybridMesh || plotGlobalDomainSchematic
%     set(0, 'DefaultFigureVisible', 'on');
% end
% 
% fprintf('\nRunning Softlubrication hybrid gap-1D / exterior-2D pressure case\n');
% fprintf('   dt                  = %.6e s\n', dt);
% fprintf('   tEnd                = %.6e s\n', tEnd);
% fprintf('   requested steps     = %.0f\n', nSteps);
% fprintf('   z-domain            = [%.3f, %.3f] um\n', zMin*1e6, zMax*1e6);
% fprintf('   exterior radius     = %.3f um\n', rOuter*1e6);
% fprintf('   1D gap window       = [%.3f, %.3f] um\n', gap1DWindow(1)*1e6, gap1DWindow(2)*1e6);
% fprintf('   exterior 2D radial Nr = %d\n\n', NrExterior2D);
% 
% out = softlube_run_case_global_coupled(cfg);

%% 10. Compact Final Native 2D Fields
out.native2D = collect_final_native2d(out);
out.firstStepHybridMesh = collect_hybrid_mesh_at_step( ...
    out, 1, gap1DWindow, NrExterior2D);
out.hybridMesh = collect_final_hybrid_mesh(out, gap1DWindow, NrExterior2D);

if plotFirstStepHybridMesh
    plot_hybrid_mesh_at_step(out, out.firstStepHybridMesh, 1);
end

if plotNative2DPressure
    plot_final_native2d_pressure(out);
end

if plotNative2DVelocity
    plot_final_native2d_velocity(out);
end

if plotGlobalDomainSchematic
    out.postPlots.globalDomain = softlube_plot_global1D_domain(cfg, out);
end

%% 11. Summary
fprintf('\nHybrid gap-1D / exterior-2D pressure run finished\n');
fprintf('   stopStep = %d\n', out.stopStep);
fprintf('   final t  = %.6e s\n', out.t(end));

if isfield(out.native2D, 'P') && any(isfinite(out.native2D.P(:)))
    fprintf('   native P range = [%.6e, %.6e] Pa\n', ...
        min(out.native2D.P(:), [], 'omitnan'), ...
        max(out.native2D.P(:), [], 'omitnan'));
end

if isfield(out, 'p2DMaxHist') && any(isfinite(out.p2DMaxHist))
    fprintf('   final max |P2D| = %.6e Pa\n', out.p2DMaxHist(end));
end

%% Local Plot Helpers
function native2D = collect_final_native2d(out)
    native2D = struct();
    if ~isfield(out, 'stopStep') || out.stopStep < 1 || ...
            ~isfield(out, 'PHist') || isempty(out.PHist)
        return;
    end

    k = out.stopStep;
    native2D.P = out.PHist(:,:,k);
    native2D.R = out.RPHist(:,:,k);
    native2D.Z = out.ZPHist(:,:,k);
    native2D.ur = out.urCHist(:,:,k);
    native2D.uz = out.uzCHist(:,:,k);
    native2D.speed = out.speedCHist(:,:,k);
    native2D.note = ['Hybrid pressure field: 1D lubrication in the gap ', ...
        'and body-fitted 2D MAC/Stokes in the exterior blocks. ', ...
        'Coordinates are physical r,z locations of pressure-cell centers.'];
end

function hybridMesh = collect_final_hybrid_mesh(out, gap1DWindow, NrExterior2D)
    hybridMesh = struct();
    if ~isfield(out, 'stopStep') || out.stopStep < 1
        return;
    end
    hybridMesh = collect_hybrid_mesh_at_step( ...
        out, out.stopStep, gap1DWindow, NrExterior2D);
end

function hybridMesh = collect_hybrid_mesh_at_step(out, stepIdx, gap1DWindow, NrExterior2D)
    hybridMesh = struct();
    if nargin < 2 || isempty(stepIdx)
        stepIdx = 1;
    end
    stepIdx = round(stepIdx);

    if ~isfield(out, 'stopStep') || out.stopStep < stepIdx || ...
            stepIdx < 1 || ~isfield(out, 'fluidHist') || ...
            numel(out.fluidHist) < stepIdx || isempty(out.fluidHist{stepIdx})
        return;
    end

    statePlot = state_for_plot_at_step(out, stepIdx);

    opts = struct();
    opts.gapZ = gap1DWindow;
    opts.NrExterior = NrExterior2D + 1;
    hybridMesh = softlube_build_hybrid_gap1d_exterior2d_mesh( ...
        out, statePlot, out.fluidHist{stepIdx}, opts);
    hybridMesh.stepIndex = stepIdx;
    if isfield(out, 't') && numel(out.t) >= stepIdx
        hybridMesh.t = out.t(stepIdx);
    end
end

function plot_hybrid_mesh_at_step(out, hybridMesh, stepIdx)
    if nargin < 3 || isempty(stepIdx)
        stepIdx = 1;
    end
    if isempty(fieldnames(hybridMesh))
        warning('No hybrid mesh is available at step %d to plot.', stepIdx);
        return;
    end

    statePlot = state_for_plot_at_step(out, stepIdx);
    fluidPlot = [];
    if isfield(out, 'fluidHist') && numel(out.fluidHist) >= stepIdx
        fluidPlot = out.fluidHist{stepIdx};
    end

    figure('Name', sprintf('Hybrid mesh step %d', stepIdx), 'Color', 'w');
    ax = gca;
    hold(ax, 'on');
    set(ax, 'FontSize', 18);

    hLegend = gobjects(0);
    [drewRuntimeBlocks, hRuntime] = draw_runtime_hybrid_blocks(ax, fluidPlot);
    hLegend = [hLegend; hRuntime(:)];
    if ~drewRuntimeBlocks
        hUp = draw_hybrid_summary_block( ...
            ax, hybridMesh.exterior2D.upstream, ...
            [0.40 0.55 0.75], 0.45, '2D exterior mesh');
        hDown = draw_hybrid_summary_block( ...
            ax, hybridMesh.exterior2D.downstream, ...
            [0.40 0.55 0.75], 0.45, '2D exterior mesh');
        if ~isempty(hUp) && isgraphics(hUp)
            hLegend = [hLegend; hUp];
            if ~isempty(hDown) && isgraphics(hDown)
                set(hDown, 'HandleVisibility', 'off');
            end
        elseif ~isempty(hDown) && isgraphics(hDown)
            hLegend = [hLegend; hDown];
        end
    end

    hGap = draw_hybrid_gap_mesh(ax, hybridMesh);
    hLegend = [hLegend; hGap(:)];
    hBoundary = draw_fluid_boundaries_for_state(ax, out, statePlot);
    hLegend = [hLegend; hBoundary(:)];

    xlabel(ax, 'r [\mum]');
    ylabel(ax, 'z [\mum]');
    if isfield(hybridMesh, 't') && isfinite(hybridMesh.t)
        title(ax, sprintf('Hybrid mesh at first time step: step %d, t = %.4g s', ...
            stepIdx, hybridMesh.t));
    else
        title(ax, sprintf('Hybrid mesh at first time step: step %d', stepIdx));
    end
    axis(ax, 'equal');
    axis(ax, 'tight');
    box(ax, 'on');

    hLegend = hLegend(isgraphics(hLegend));
    if ~isempty(hLegend)
        legend(ax, hLegend, get(hLegend, 'DisplayName'), ...
            'Location', 'best');
    end
end

function [drew, hLegend] = draw_runtime_hybrid_blocks(ax, fluidPlot)
    drew = false;
    hLegend = gobjects(0);
    if ~isstruct(fluidPlot) || ~isfield(fluidPlot, 'meshF') || ...
            ~isfield(fluidPlot.meshF, 'blocks')
        return;
    end

    blocks = fluidPlot.meshF.blocks;
    hUp = draw_runtime_hybrid_block(ax, blocks.upstream, ...
        [0.40 0.55 0.75], 0.45, '2D exterior mesh');
    hDown = draw_runtime_hybrid_block(ax, blocks.downstream, ...
        [0.40 0.55 0.75], 0.45, '2D exterior mesh');
    hCenters = draw_runtime_pressure_centers(ax, blocks);

    drew = (~isempty(hUp) && isgraphics(hUp)) || ...
           (~isempty(hDown) && isgraphics(hDown));
    if ~isempty(hUp) && isgraphics(hUp)
        hBlock = hUp;
        if ~isempty(hDown) && isgraphics(hDown)
            set(hDown, 'HandleVisibility', 'off');
        end
    else
        hBlock = hDown;
    end

    hNew = [hBlock; hCenters];
    hNew = hNew(isgraphics(hNew));
    if ~isempty(hNew)
        hLegend = hNew(:);
    end
end

function hFirst = draw_runtime_hybrid_block(ax, block, color, lineWidth, displayName)
    hFirst = gobjects(0);
    if ~isstruct(block) || ~isfield(block, 'mesh') || ~isstruct(block.mesh)
        return;
    end

    mesh = block.mesh;
    if isfield(mesh, 'Rzf') && isfield(mesh, 'zF') && ~isempty(mesh.Rzf)
        R = mesh.Rzf;
        Z = repmat(mesh.zF(:).', size(R, 1), 1);
    elseif isfield(mesh, 'Rp') && isfield(mesh, 'Zp') && ~isempty(mesh.Rp)
        R = mesh.Rp;
        Z = mesh.Zp;
    else
        return;
    end

    hLines = plot_mesh_matrix(ax, R, Z, color, lineWidth);
    if isempty(hLines)
        return;
    end
    hFirst = hLines(1);
    set(hFirst, 'DisplayName', displayName);
    if numel(hLines) > 1
        set(hLines(2:end), 'HandleVisibility', 'off');
    end
end

function hFirst = draw_runtime_pressure_centers(ax, blocks)
    hFirst = gobjects(0);
    blockNames = {'upstream', 'downstream'};
    for i = 1:numel(blockNames)
        name = blockNames{i};
        if ~isfield(blocks, name) || ~isstruct(blocks.(name)) || ...
                ~isfield(blocks.(name), 'mesh')
            continue;
        end
        mesh = blocks.(name).mesh;
        if ~isfield(mesh, 'Rp') || ~isfield(mesh, 'Zp') || isempty(mesh.Rp)
            continue;
        end
        h = plot(ax, mesh.Rp(:)*1e6, mesh.Zp(:)*1e6, '.', ...
            'Color', [0.10 0.25 0.80], 'MarkerSize', 4, ...
            'DisplayName', 'pressure centers');
        if isempty(hFirst)
            hFirst = h;
        else
            set(h, 'HandleVisibility', 'off');
        end
    end
end

function hFirst = draw_hybrid_summary_block(ax, block, color, lineWidth, displayName)
    hFirst = gobjects(0);
    if ~isstruct(block) || ~isfield(block, 'R') || ~isfield(block, 'Z') || ...
            isempty(block.R)
        return;
    end

    hLines = plot_mesh_matrix(ax, block.R, block.Z, color, lineWidth);
    if isempty(hLines)
        return;
    end
    hFirst = hLines(1);
    set(hFirst, 'DisplayName', displayName);
    if numel(hLines) > 1
        set(hLines(2:end), 'HandleVisibility', 'off');
    end
end

function hLines = plot_mesh_matrix(ax, R, Z, color, lineWidth)
    hLines = gobjects(0);
    if isempty(R) || isempty(Z)
        return;
    end
    hRadial = plot(ax, R*1e6, Z*1e6, '-', ...
        'Color', color, 'LineWidth', lineWidth);
    hAxial = plot(ax, R.'*1e6, Z.'*1e6, '-', ...
        'Color', color, 'LineWidth', lineWidth);
    hLines = [hRadial(:); hAxial(:)];
end

function h = draw_hybrid_gap_mesh(ax, hybridMesh)
    h = gobjects(0);
    if ~isfield(hybridMesh, 'gap1D') || ~isfield(hybridMesh.gap1D, 'z') || ...
            isempty(hybridMesh.gap1D.z)
        return;
    end

    gap = hybridMesh.gap1D;
    h(end+1,1) = plot(ax, gap.rL*1e6, gap.z*1e6, '-', ...
        'Color', [0.05 0.05 0.05], 'LineWidth', 1.2, ...
        'DisplayName', '1D gap boundaries');
    hEnd = plot(ax, gap.rE*1e6, gap.z*1e6, '-', ...
        'Color', [0.05 0.05 0.05], 'LineWidth', 1.2);
    set(hEnd, 'HandleVisibility', 'off');

    h(end+1,1) = plot(ax, gap.rMid*1e6, gap.z*1e6, 'o--', ...
        'Color', [0.85 0.20 0.10], 'MarkerFaceColor', [0.85 0.20 0.10], ...
        'MarkerSize', 3.5, 'LineWidth', 0.9, ...
        'DisplayName', '1D gap nodes');

    if isfield(hybridMesh, 'gapZ') && isfield(hybridMesh, 'rDomain')
        rSpan = hybridMesh.rDomain(:).';
        for i = 1:numel(hybridMesh.gapZ)
            hCouple = plot(ax, rSpan*1e6, ...
                hybridMesh.gapZ(i)*[1 1]*1e6, ':', ...
                'Color', [0.15 0.15 0.15], 'LineWidth', 1.0);
            set(hCouple, 'HandleVisibility', 'off');
        end
    end
end

function h = draw_fluid_boundaries_for_state(ax, out, statePlot)
    h = gobjects(0);
    if ~isfield(out, 'z') || isempty(out.z)
        return;
    end
    z = out.z(:);

    if isfield(out, 'meshL') && isfield(out, 'interfaceL') && ...
            isfield(statePlot, 'uL') && ~isempty(statePlot.uL)
        [rLDef, zLDef] = local_deformed_interface_curve( ...
            out.meshL, statePlot.uL, out.interfaceL);
        h(end+1,1) = plot(ax, rLDef*1e6, zLDef*1e6, 'k-', ...
            'LineWidth', 1.4, 'DisplayName', 'solid/fluid boundaries');
    elseif isfield(statePlot, 'deltaL')
        h(end+1,1) = plot(ax, statePlot.deltaL(:)*1e6, z*1e6, 'k-', ...
            'LineWidth', 1.4, 'DisplayName', 'solid/fluid boundaries');
    end

    if isfield(out, 'meshE') && isfield(out, 'interfaceE') && ...
            isfield(statePlot, 'uE') && ~isempty(statePlot.uE)
        [rEDef, zEDef] = local_deformed_interface_curve( ...
            out.meshE, statePlot.uE, out.interfaceE);
        hE = plot(ax, rEDef*1e6, zEDef*1e6, 'k-', 'LineWidth', 1.4);
        set(hE, 'HandleVisibility', 'off');
    elseif isfield(statePlot, 'deltaE')
        hE = plot(ax, statePlot.deltaE(:)*1e6, z*1e6, 'k-', ...
            'LineWidth', 1.4);
        set(hE, 'HandleVisibility', 'off');
    end
end

function plot_final_native2d_pressure(out)
    native2D = out.native2D;
    if ~isfield(native2D, 'P') || isempty(native2D.P) || ...
            ~any(isfinite(native2D.P(:)))
        warning('No finite native 2D pressure field is available to plot.');
        return;
    end

    figure;
    ax = gca;
    set(ax, 'FontSize', 22);
    if is_hybrid_fluid_plot(out)
        contour_masked_regular_field(ax, out, native2D.P, 1.0, 48);
    else
        contourf(ax, native2D.R*1e6, native2D.Z*1e6, native2D.P, ...
            48, 'LineColor', 'none');
    end
    colorbar;
    hold(ax, 'on');
    draw_final_solid_blanks(ax, out);
    draw_fluid_boundaries(ax, out);
    xlabel(ax, 'r [\mum]');
    ylabel(ax, 'z [\mum]');
    title(ax, sprintf('Hybrid pressure P(r,z), t = %.4g s', ...
        out.t(out.stopStep)));
    axis(ax, 'equal');
    axis(ax, 'tight');
    box(ax, 'on');
end

function plot_final_native2d_velocity(out)
    native2D = out.native2D;
    if ~isfield(native2D, 'uz') || isempty(native2D.uz) || ...
            ~any(isfinite(native2D.uz(:)))
        warning('No finite native 2D velocity field is available to plot.');
        return;
    end

    figure;
    ax = gca;
    set(ax, 'FontSize', 22);
    if is_hybrid_fluid_plot(out)
        contour_masked_regular_field(ax, out, native2D.uz, 1e6, 48);
    else
        contourf(ax, native2D.R*1e6, native2D.Z*1e6, native2D.uz*1e6, ...
            48, 'LineColor', 'none');
    end
    colorbar;
    hold(ax, 'on');
    draw_final_solid_blanks(ax, out);
    draw_fluid_boundaries(ax, out);
    xlabel(ax, 'r [\mum]');
    ylabel(ax, 'z [\mum]');
    title(ax, sprintf('Hybrid axial velocity u_z(r,z), t = %.4g s', ...
        out.t(out.stopStep)));
    axis(ax, 'equal');
    axis(ax, 'tight');
    box(ax, 'on');
end

function tf = is_hybrid_fluid_plot(out)
    tf = false;
    if ~isfield(out, 'stopStep') || out.stopStep < 1 || ...
            ~isfield(out, 'fluidHist') || numel(out.fluidHist) < out.stopStep || ...
            isempty(out.fluidHist{out.stopStep})
        return;
    end
    fluidPlot = out.fluidHist{out.stopStep};
    tf = isfield(fluidPlot, 'meshType') && strcmpi(fluidPlot.meshType, 'hybrid_gap1d_exterior2d');
end

function contour_masked_regular_field(ax, out, fieldNative, valueScale, nLevels)
% Plot the hybrid field on a regular r-z image grid and mask the deformed
% solids. This avoids drawing artificial curvilinear cells that connect the
% 1D gap strip directly to the 2D exterior reservoir across solid caps.
    native2D = out.native2D;
    valid = isfinite(native2D.R) & isfinite(native2D.Z) & isfinite(fieldNative);
    if nnz(valid) < 3
        warning('Not enough finite hybrid field points to plot.');
        return;
    end

    rVals = native2D.R(valid);
    zVals = native2D.Z(valid);
    fVals = valueScale * fieldNative(valid);

    nrPlot = 260;
    nzPlot = 260;
    rGrid = linspace(min(rVals), max(rVals), nrPlot);
    zGrid = linspace(min(zVals), max(zVals), nzPlot);
    [RGrid, ZGrid] = meshgrid(rGrid, zGrid);

    F = scatteredInterpolant(rVals(:), zVals(:), fVals(:), 'linear', 'none');
    FGrid = F(RGrid, ZGrid);

    statePlot = final_state_for_plot(out);
    solidMask = deformed_solids_mask_on_grid(out, statePlot, RGrid, ZGrid);
    FGrid(solidMask) = NaN;

    finiteVals = FGrid(isfinite(FGrid));
    if isempty(finiteVals)
        warning('No finite hybrid field values remain after solid masking.');
        return;
    end
    if max(finiteVals) > min(finiteVals)
        contourf(ax, RGrid*1e6, ZGrid*1e6, FGrid, nLevels, 'LineColor', 'none');
    else
        contourf(ax, RGrid*1e6, ZGrid*1e6, FGrid, 1, 'LineColor', 'none');
    end
end

function draw_final_solid_blanks(ax, out)
    statePlot = out.state;
    if isfield(out, 'stateHist') && out.stopStep >= 1 && ...
            numel(out.stateHist) >= out.stopStep && ~isempty(out.stateHist{out.stopStep})
        statePlot = out.stateHist{out.stopStep};
    end

    if isfield(out, 'meshL') && isfield(statePlot, 'uL') && ...
            ~isempty(out.meshL) && ~isempty(statePlot.uL)
        draw_deformed_solid_patch(ax, out.meshL, statePlot.uL, [0.05 0.35 0.12]);
    end

    if isfield(out, 'meshE') && isfield(statePlot, 'uE') && ...
            ~isempty(out.meshE) && ~isempty(statePlot.uE)
        draw_deformed_solid_patch(ax, out.meshE, statePlot.uE, [0.00 0.15 0.65]);
    end
end

function draw_deformed_solid_patch(ax, mesh, u, edgeColor)
    rDef = mesh.nodes(:,1) + u(1:2:end);
    zDef = mesh.nodes(:,2) + u(2:2:end);
    patch(ax, 'Faces', mesh.conn, ...
        'Vertices', [rDef, zDef]*1e6, ...
        'FaceColor', 'w', ...
        'EdgeColor', edgeColor, ...
        'LineWidth', 0.12);
end

function draw_fluid_boundaries(ax, out)
    statePlot = final_state_for_plot(out);
    z = out.z(:);

    drewL = false;
    if isfield(out, 'meshL') && isfield(out, 'interfaceL') && ...
            isfield(statePlot, 'uL') && ~isempty(statePlot.uL)
        [rLDef, zLDef] = local_deformed_interface_curve( ...
            out.meshL, statePlot.uL, out.interfaceL);
        plot(ax, rLDef*1e6, zLDef*1e6, 'k-', 'LineWidth', 1.2);
        drewL = true;
    end
    if ~drewL && isfield(statePlot, 'deltaL')
        plot(ax, statePlot.deltaL(:)*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
    end

    drewE = false;
    if isfield(out, 'meshE') && isfield(out, 'interfaceE') && ...
            isfield(statePlot, 'uE') && ~isempty(statePlot.uE)
        [rEDef, zEDef] = local_deformed_interface_curve( ...
            out.meshE, statePlot.uE, out.interfaceE);
        plot(ax, rEDef*1e6, zEDef*1e6, 'k-', 'LineWidth', 1.2);
        drewE = true;
    end
    if ~drewE && isfield(statePlot, 'deltaE')
        plot(ax, statePlot.deltaE(:)*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
    end

    rMax = max(out.native2D.R(:), [], 'omitnan');
    zMinPlot = min(z);
    zMaxPlot = max(z);
    plot(ax, [0 rMax]*1e6, [zMinPlot zMinPlot]*1e6, 'k--', 'LineWidth', 1.0);
    plot(ax, [0 rMax]*1e6, [zMaxPlot zMaxPlot]*1e6, 'k--', 'LineWidth', 1.0);
end

function statePlot = final_state_for_plot(out)
    statePlot = state_for_plot_at_step(out, out.stopStep);
end

function statePlot = state_for_plot_at_step(out, stepIdx)
    statePlot = struct();
    if isfield(out, 'state')
        statePlot = out.state;
    end
    if isfield(out, 'stateHist') && stepIdx >= 1 && ...
            numel(out.stateHist) >= stepIdx && ~isempty(out.stateHist{stepIdx})
        statePlot = out.stateHist{stepIdx};
    end
end

function solidMask = deformed_solids_mask_on_grid(out, statePlot, RGrid, ZGrid)
    solidMask = false(size(RGrid));
    if isfield(out, 'meshL') && isfield(statePlot, 'uL') && ...
            ~isempty(out.meshL) && ~isempty(statePlot.uL)
        solidMask = solidMask | local_deformed_solid_mask( ...
            out.meshL, statePlot.uL, RGrid, ZGrid);
    end
    if isfield(out, 'meshE') && isfield(statePlot, 'uE') && ...
            ~isempty(out.meshE) && ~isempty(statePlot.uE)
        solidMask = solidMask | local_deformed_solid_mask( ...
            out.meshE, statePlot.uE, RGrid, ZGrid);
    end
end

function solidMask = local_deformed_solid_mask(mesh, u, RGrid, ZGrid)
    solidMask = false(size(RGrid));
    if isempty(mesh) || isempty(u) || ~isfield(mesh, 'nodes') || ...
            ~isfield(mesh, 'conn') || numel(u) < 2*size(mesh.nodes,1)
        return;
    end

    rDef = mesh.nodes(:,1) + u(1:2:end);
    zDef = mesh.nodes(:,2) + u(2:2:end);
    for e = 1:size(mesh.conn,1)
        ids = mesh.conn(e,:);
        rv = rDef(ids);
        zv = zDef(ids);
        inBox = RGrid >= min(rv) & RGrid <= max(rv) & ...
                ZGrid >= min(zv) & ZGrid <= max(zv);
        if any(inBox(:))
            localMask = false(size(solidMask));
            localMask(inBox) = inpolygon(RGrid(inBox), ZGrid(inBox), rv, zv);
            solidMask = solidMask | localMask;
        end
    end
end

function [rDef, zDef] = local_deformed_interface_curve(mesh, u, interfaceNodes)
    ids = interfaceNodes(:);
    rDef = mesh.nodes(ids,1) + u(2*ids - 1);
    zDef = mesh.nodes(ids,2) + u(2*ids);
    [zDef, idx] = sort(zDef);
    rDef = rDef(idx);
end
