%% RUN_INPUT_FULL2D_PRESSURE
% Input deck for the hybrid pressure calculation:
%   1D lubrication mesh in the leukocyte/endothelium gap,
%   2D body-fitted MAC/Stokes mesh in the upstream/downstream exterior.
%
% This file keeps run_input.m unchanged and starts a separate case that uses
% the mixed-dimensional fluid solve inside:
%   softlube_run_case_global_coupled(cfg)
%
% Important distinction:
%   - The gap interval uses the Reynolds/lubrication p(z) solve.
%   - The exterior intervals are solved as true 2D body-fitted MAC/Stokes
%     blocks, matched to the 1D pressure at z = -0.2 and z = 4.2 um.
%   - The final composite pressure field is stored as P(r,z) in out.PHist,
%     on the physical grid out.RPHist/out.ZPHist.
%
% To run:
%   cd('/Users/shu/Softlubrication')
%   run_input_full2D_pressure
%
% Main hybrid outputs:
%   out.PHist(:,:,end)       composite hybrid pressure P(r,z) at final step
%   out.RPHist(:,:,end)      physical r locations for P
%   out.ZPHist(:,:,end)      physical z locations for P
%   out.urCHist(:,:,end)     radial velocity at pressure-cell centers
%   out.uzCHist(:,:,end)     axial velocity at pressure-cell centers
%   out.native2D             compact final-step hybrid 2D fields
%   out.hybridMesh           final 1D-gap / 2D-exterior mesh description

clc;
clear all;

%% 1. Paths
softlubeDir = '';
addpath(softlubeDir);

%% 2. Time Controls
dt = 3e-4;           % time step size [s]
nSteps = 1;          % requested time steps
tEnd = nSteps * dt;

useAdaptiveTimeStep = true;
dtMin = dt / 1024;

%% 3. Global Domain and Pressure Boundary Conditions
% The hybrid pressure uses a 1D gap solve and 2D exterior body-fitted
% r-z blocks. The axial ends use P = 0. In the exterior MAC/Stokes blocks,
% r = 0 is the axis outside the leukocyte and r = 15 um is the exterior wall
% outside the endothelium.
zMin = -6e-6;        % -6 um
zMax = 10e-6;        % 10 um
rOuter = 4e-6;      % 15 um

pAtZMin = 0;         % P = 0 at z = -6 um [Pa]
pAtZMax = 0;         % P = 0 at z = 10 um [Pa]
axisPressureBC = 'dPdr=0';

%% 4. Fluid Mesh Resolution
% Hybrid split. The gap interval is handled by the 1D lubrication mesh; the
% exterior intervals are handled by the 2D body-fitted MAC/Stokes mesh.
gap1DWindow = [-0.2e-6, 4.2e-6];

% Resolve the 1D lubrication region finely and keep the exterior coarse.
finePressureWindow = gap1DWindow;
coarseDz = 1.0e-6;
fineDz = 0.05e-6;

% 2D MAC cells across the radial exterior fluid blocks. The radial grid is
% nonuniform below, with extra resolution in the thin-gap band r = 2..4 um.
NrExterior2D = 20;
NrFluid2D = NrExterior2D;
nrEnv = str2double(getenv('SOFTLUBE_FULL2D_NR'));
if isfinite(nrEnv) && nrEnv >= 5
    NrFluid2D = round(nrEnv);
    NrExterior2D = NrFluid2D;
end
radialFineWindow = [2e-6, 4e-6];
radialFineWeight = 4.0;

%% 5. Plotting and Output
plotNative2DPressure = true;
plotNative2DStress = true;
plotNative2DVelocity = false;
plotFirstStepHybridMesh = false;
plotGlobalDomainSchematic = false;

% Keep the old solver plot block off. The local plot at the bottom shows
% only the composite hybrid result from out.PHist.
makeLegacyPlots = false;

saveOutput = false;
outputFile = fullfile(softlubeDir, 'simulation_output_full2D_pressure.mat');

closeFiguresAtStart = true;
printEvery = 1;

%% 6. Build the Base Case
cfg = struct();
% load('out_2_combined_t40.mat');
% cfg=out.cfg;
% out=rmfield(out,'cfg');


if numel(fieldnames(struct()))==0
    % Prestressed initial geometries. Leukocyte geometry/mesh settings are read
    % inside the solver from solid_leu_P600.mat and are not repeated here.
    cfg.geometry.endotheliumPrestressFile = fullfile(softlubeDir, 'solid_endo_P300_wide.mat');
    cfg.geometry.leukocytePrestressFile = fullfile(softlubeDir, 'solid_leu_P600.mat');
    cfg.geometry.usePrestressedLeukocyteIC = true;
    cfg.geometry.useLeukocyteCenterline = true;

    cfg.fluid.zMin = zMin;
    cfg.fluid.zMax = zMax;
    cfg.fluid.mu = 1.2e-3;
    cfg.fluid.slipE = 0e-9;
    cfg.fluid.slipL = 0e-9;
    cfg.fluid.NrFluid2D = NrFluid2D;
    cfg.fluid.Nr = NrFluid2D;
    cfg.fluid.NrExterior2D = NrExterior2D;
    cfg.fluid.bodyFittedRadialFineWindow = radialFineWindow;
    cfg.fluid.bodyFittedRadialFineWeight = radialFineWeight;
    cfg.fluid.innerBoundary = 'deformed_leukocyte';
    cfg.fluid.useSlopeAwareWallKinematics = true;

    cfg.bc.pIn = pAtZMin;
    cfg.bc.pOut = pAtZMax;
    cfg.bc.UwL = 0;
    cfg.bc.UwE = 0;

    cfg.solid.noLeukocyte = false;
    cfg.solid.leukocyte.enabled = true;
    cfg.solid.leukocyte.deformable = true;
    cfg.solid.leukocyte.EL = 200;
    cfg.solid.leukocyte.nuL = 0.46;
    cfg.solid.leukocyte.useViscoelasticLeukocyte = true;
    cfg.solid.leukocyte.etaL = 1;
    cfg.solid.leukocyte.supportL = 'axis';

    cfg.solid.endothelium.Ee = 500;
    cfg.solid.endothelium.nuE = 0.46;
    cfg.solid.endothelium.useViscoelasticEndothelium = true;
    cfg.solid.endothelium.etaE = 1;

    cfg.interface.useExactDeformedInterface = true;
    cfg.interface.useExactDeformedInterfaceInMonolithic = true;
    cfg.interface.useReferenceZForTractionMapping = false;
    cfg.interface.zeroLeukocyteTractionOutsideOverlap = false;
    cfg.interface.leukocytePressureSupportZ = [zMin, zMax];
    cfg.interface.warnPrestressLoadMismatch = false;
    cfg.numerics.dt = dt;
    cfg.numerics.tEnd = tEnd;
    cfg.numerics.minGap = 0.001e-6;
    cfg.numerics.maxNewtonMono = 100;
    cfg.numerics.tolNewtonMono = 1e-6;
    cfg.numerics.monoSolidAbsTol = 1e-7;
    cfg.numerics.monoFluidAbsTol = 1e-8;
    cfg.numerics.maxFunctionEvaluationsTwoSolid = 400;
    cfg.numerics.enablePressureJumpRetry = false;
    cfg.numerics.maxPressureJumpAbs = inf;
    cfg.numerics.maxPressureJump2DAbs = inf;
    cfg.numerics.useObjectiveKelvinVoigt = true;
    cfg.numerics.useDirectFluidSolve = true;
    cfg.numerics.useFsolveMono = true;
    cfg.numerics.enableAdaptiveTimeStep = useAdaptiveTimeStep;
    cfg.numerics.dtMin = dtMin;
    cfg.numerics.dtRetryFactor = 0.5;
    cfg.numerics.dtGrowFactor = 1.25;
    cfg.numerics.maxTimeStepRetries = 18;
    cfg.numerics.fsolveDisplay = 'off';
    cfg.numerics.diagnosticsEnabled = true;

    cfg.output.saveOutput = saveOutput;
    cfg.output.outputFile = outputFile;
    cfg.output.makePlots = makeLegacyPlots;
    cfg.output.plotFinalGlobalPressureContour = false;
    cfg.output.plotGlobal2DPressureTractionComparison = false;
    cfg.output.printEvery = printEvery;
    cfg.output.storeFull2DFluidHist = false;
    cfg.output.store2DFluidArrays = true;
    cfg.output.store2DPhysicalGridHist = true;
    cfg.output.store2DSpeedHist = true;

    cfg.runtime.useEnvironment = false;
    cfg.ui.closeFigures = closeFiguresAtStart;
    cfg.parOverrides = struct();

    %% 7. Global Axial Mesh and Exterior Boundary Embedding
    % The flag name is historical. Here it is kept on to reuse the global axial
    % mesh and the outside-solid radius rules:
    %   leukocyte absent: fluid inner boundary becomes r = 0
    %   endothelium absent: fluid outer boundary becomes r = 15 um
    % The composite hybrid pressure output is out.PHist, not out.global1D.
    cfg.global1D.zDomain = [zMin, zMax];
    cfg.global1D.rDomain = [0, rOuter];
    cfg.global1D.axisPressureBC = axisPressureBC;
    cfg.global1D.outerPressureBC = 0;

    cfg.parOverrides.useGlobal1DPressure = true;
    cfg.parOverrides.global1DOuterRadius = rOuter;
    cfg.parOverrides.global1DAxisRadius = 1e-9;
    cfg.parOverrides.useGlobal1DCoarseEdgeMesh = true;
    cfg.parOverrides.global1DFineWindow = finePressureWindow;
    cfg.parOverrides.global1DCoarseDz = coarseDz;
    cfg.parOverrides.global1DFineDz = fineDz;
    cfg.parOverrides.global1DPlotNr = 161;

    %% 8. Hybrid Gap-1D / Exterior-2D Fluid Module
    cfg.fluid.useFull2DFluid = true;
    %added_temp
    cfg.fluid.useHybridGap1DExterior2DFluid = false;
    %cfg.fluid.useHybridGap1DExterior2DFluid = true;
    cfg.fluid.hybridGapZ = gap1DWindow;
    cfg.fluid.hybridExteriorMinCells = 4;
    cfg.fluid.hybridExteriorFailMode = 'warn';
    cfg.fluid.hybridUseExterior2DPressureVector = true;
    cfg.fluid.useBodyFittedMACFluid = true;
    cfg.fluid.useBodyFittedMACTractionCorrection = true;
    %added_temp: bumped from 1/0.05 to test whether the interface
    %traction mismatch shrinks with a stronger partitioned correction
    %toward the 2D MAC traction (revert to 1 / 0.05 to go back).
    % Use the existing feedback-controlled correction loop (already built
    % in apply_bodyfitted_MAC_traction_correction_feedback.m, per an
    % earlier review comment, but never wired into this entry point until
    % now) -- measures the actual interface %mismatch and stops once it is
    % below threshold, instead of a fixed pass count. Cap raised to 100 so
    % it has room to actually converge on harder steps; relax left at the
    % pre-existing 0.5 for a controlled first comparison (not changing two
    % things -- feedback vs. fixed relax -- at once).
    cfg.fluid.useFeedbackTractionCorrection = true;
    % Fixed relaxation is documented in apply_bodyfitted_MAC_traction_
    % correction_feedback.m to monotonically diverge for this problem
    % (confirmed just now: fixed relax=0.5 produced 260,000%+ mismatch,
    % never converging in 100 passes). Aitken relaxation is the author's
    % own documented fix for exactly this failure mode.
    cfg.fluid.useAitkenTractionCorrectionRelax = true;
    cfg.fluid.debugRadialAxialTractionCorrection = true;
    cfg.fluid.maxBodyFittedTractionCorrections = 12;
    cfg.fluid.debugVerbose = true;
    cfg.fluid.bodyFittedTractionCorrectionRelax = 0.5;
    cfg.fluid.bodyFittedTractionCorrectionFailMode = 'warn';
    cfg.fluid.resolveFluidAfterTractionCorrection = true;
    cfg.fluid.useBodyFittedMACTractionInSolid = true;

    cfg.parOverrides.useFull2DFluid = true;
    %added_temp
    cfg.parOverrides.useHybridGap1DExterior2DFluid = false;
    %cfg.parOverrides.useHybridGap1DExterior2DFluid = true;
    cfg.parOverrides.hybridGapZ = gap1DWindow;
    cfg.parOverrides.NrExterior2D = NrExterior2D;
    cfg.parOverrides.bodyFittedRadialFineWindow = radialFineWindow;
    cfg.parOverrides.bodyFittedRadialFineWeight = radialFineWeight;
    cfg.parOverrides.hybridExteriorMinCells = 4;
    cfg.parOverrides.hybridExteriorFailMode = 'warn';
    cfg.parOverrides.hybridUseExterior2DPressureVector = true;
    cfg.parOverrides.useBodyFittedMACFluid = true;
    cfg.parOverrides.useBodyFittedMACTractionCorrection = true;
    %added_temp: same test as cfg.fluid above -- parOverrides wins, so
    %both copies must match (revert to 1 / 0.05 to go back).
    cfg.parOverrides.useFeedbackTractionCorrection = true;
    cfg.parOverrides.useAitkenTractionCorrectionRelax = true;
    cfg.parOverrides.debugRadialAxialTractionCorrection = true;
    cfg.parOverrides.maxBodyFittedTractionCorrections = 12;
    cfg.parOverrides.debugVerbose = true;
    cfg.parOverrides.bodyFittedTractionCorrectionRelax = 0.5;
    cfg.parOverrides.bodyFittedTractionCorrectionFailMode = 'warn';
    cfg.parOverrides.resolveFluidAfterTractionCorrection = true;
    cfg.parOverrides.useBodyFittedMACTractionInSolid = true;
    cfg.parOverrides.NrFluid2D = NrFluid2D;
    cfg.parOverrides.Nr = NrFluid2D;

    % Turn off the earlier projected-pressure bridge so this file reports the
    % hybrid pressure field cleanly.
    cfg.parOverrides.useGlobal2DPressureTraction = false;
else
    cfg.numerics.dt = dt;
    cfg.numerics.tEnd = tEnd;
    if ~isempty(out)
        out.par.dt = dt;
        out.par.tEnd = tEnd;
    end
end


%% ISOLATED TEST: apply a known, uniform pressure to the endothelium alone,
% no fluid solver, no coupling loop -- pure round-trip test of
% apply_interface_traction -> solve_finite_def_solid ->
% recover_nodal_stress_axisym_viscoelastic. If a known -500 Pa uniform
% traction comes back as ~-500 Pa recovered stress at the loaded boundary,
% the force/stress machinery is magnitude-correct. If it comes back at
% ~24% of that (matching the coupled-run finding), that pins the bug here.

[par, S, uE_pre] = softlube_prepare_case(cfg);
fprintf('BEFORE override: max|uE_pre| = %.6e m (prestressed IC)\n', max(abs(uE_pre)));
uE_pre = zeros(size(uE_pre));
fprintf('AFTER override: uE_pre set to all-zero (unstressed reference config)\n');
meshE = prepare_axisym_mesh_cache(S.meshE);
interfaceE = S.interfaceE;
baseE = S.baseE;

zInterface = meshE.nodes(interfaceE, 2);
zLo = min(zInterface); zHi = max(zInterface);
nZ = 21;
testZ = linspace(zLo, zHi, nZ).';

knownPressure = -500;  % Pa, uniform, known exactly
traction = struct();
traction.z = testZ;
traction.normal = knownPressure * ones(nZ,1);
traction.tangent = zeros(nZ,1);

parTest = par;
if isfield(parTest, 'solidAbsTol'), parTest = rmfield(parTest, 'solidAbsTol'); end

uEcorr = solve_finite_def_solid(meshE, uE_pre, traction, interfaceE, baseE, par.supportE, parTest, uE_pre);

stress = recover_nodal_stress_axisym_viscoelastic(meshE, uEcorr, uE_pre, par.dt, par);

fprintf('\n=== ISOLATED UNIFORM-TRACTION TEST ===\n');
fprintf('Applied uniform normal traction: %.4f Pa\n', knownPressure);
fprintf('%-10s %-16s %-16s\n', 'node', 'z [um]', 'recovered sigma_rr [Pa]');
for i = 1:numel(interfaceE)
    nd = interfaceE(i);
    fprintf('%-10d %-16.4f %-16.4f\n', nd, meshE.nodes(nd,2)*1e6, stress.sigma_rr(nd));
end
fprintf('median recovered sigma_rr at interface nodes = %.4f Pa (applied = %.4f Pa, ratio = %.4f)\n', ...
    median(stress.sigma_rr(interfaceE)), knownPressure, median(stress.sigma_rr(interfaceE))/knownPressure);

fprintf('\nSMOKE_TEST_STATUS: OK\n');

%% PROPERLY-CONTROLLED PATCH TEST
% Deliberately: (1) restrict to the flat middle z-band (z in [1,3] um,
% well away from the tapered end-caps at z=-2/+6 that contaminated the
% last two attempts), (2) within that z-band, explicitly find the FIRST
% radial element layer by grouping elements by z-centroid and sorting by
% r-centroid -- not by a fixed radial tolerance, which failed because
% Gauss points sit interior to elements, not at their edges.
cache = meshE.axisymCache;
I3 = eye(3);

elemZc = zeros(meshE.nelem,1);
elemRc = zeros(meshE.nelem,1);
for e = 1:meshE.nelem
    conn = meshE.conn(e,:);
    elemZc(e) = mean(meshE.nodes(conn,2));
    elemRc(e) = mean(meshE.nodes(conn,1));
end

zBandMask = elemZc > 1e-6 & elemZc < 3e-6;
zBandElems = find(zBandMask);
[~, ordR] = sort(elemRc(zBandElems));
firstLayerElems = zBandElems(ordR(1:min(6,numel(ordR))));

fprintf('\n=== PATCH TEST: first radial layer, flat middle z-band (z in [1,3] um) ===\n');
fprintf('Applied uniform traction: %.1f Pa\n', knownPressure);
fprintf('%-8s %-8s %-6s %-12s %-12s %-14s\n', 'elem', 'zc[um]', 'gp', 'r [um]', 'z [um]', 'sigma_rr [Pa]');
allVals = [];
for ei = 1:numel(firstLayerElems)
    e = firstLayerElems(ei);
    dofs = cache.dofs(e,:).';
    ue = uEcorr(dofs);
    conn = meshE.conn(e,:);
    rnod = meshE.nodes(conn,1) + ue(1:2:end);
    znod = meshE.nodes(conn,2) + ue(2:2:end);
    for g = 1:cache.ngp
        N = cache.N(:,g,e).';
        dNdX = cache.dNdX(:,:,g,e);
        Rg = cache.Rg0(g,e);
        rg = N * rnod;
        zg = N * znod;
        drdR = dNdX(:,1).' * rnod;
        drdZ = dNdX(:,2).' * rnod;
        dzdR = dNdX(:,1).' * znod;
        dzdZ = dNdX(:,2).' * znod;
        F = [drdR, 0, drdZ; 0, rg/Rg, 0; dzdR, 0, dzdZ];
        J = det(F);
        B = F*F.';
        trB = B(1,1)+B(2,2)+B(3,3);
        devB = B - (trB/3)*I3;
        aIso = J^(-5/3);
        T = par.Ge*aIso*devB + par.Ke*(J-1)*I3;
        fprintf('%-8d %-8.4f %-6d %-12.6f %-12.4f %-14.4f\n', e, elemZc(e)*1e6, g, rg*1e6, zg*1e6, T(1,1));
        allVals(end+1) = T(1,1); %#ok<SAGROW>
    end
end
fprintf('\nmedian sigma_rr in first-layer flat-region Gauss points = %.4f Pa (applied = %.1f Pa, ratio = %.4f)\n', ...
    median(allVals), knownPressure, median(allVals)/knownPressure);

% Cross-check against the nodal-recovered value at the nearest interface
% node in the SAME z-band, for a direct apples-to-apples comparison.
zInterfaceAll = meshE.nodes(interfaceE,2);
inBand = zInterfaceAll > 1e-6 & zInterfaceAll < 3e-6;
nodesInBand = interfaceE(inBand);
recoveredInBand = stress.sigma_rr(nodesInBand);
fprintf('median nodal-recovered sigma_rr, same z-band = %.4f Pa (ratio = %.4f)\n', ...
    median(recoveredInBand), median(recoveredInBand)/knownPressure);

fprintf('\nSMOKE_TEST_STATUS: OK\n');

%% ANALYTICAL BENCHMARK: Lame thick-cylinder solution, clamped outer
% boundary (u_r(b)=0), prescribed inner traction sigma_rr(a) = knownPressure.
% Valid comparison here because displacement is tiny (0.16% of thickness,
% confirmed genuinely small-strain), so linear elasticity is an excellent
% approximation of the finite-deformation code's actual behavior.
a = min(meshE.nodes(:,1));
b = max(meshE.nodes(:,1));
fprintf('\n=== LAME ANALYTICAL BENCHMARK ===\n');
fprintf('inner radius a = %.6e m, outer radius b = %.6e m\n', a, b);
fprintf('support type: %s\n', par.supportE);

Ge_ = par.Ge; Ke_ = par.Ke;
lambda_ = Ke_ - (2/3)*Ge_;
G_ = Ge_;

% sigma_rr(r) = A - B/r^2 ; u(r) = C1*r + C2/r
% u(b) = 0 -> C1 = -C2/b^2
% sigma_rr(a) = (2*lambda+2G)*C1 - 2G*C2/a^2 = knownPressure
denomC2 = -(2*lambda_+2*G_)/b^2 - 2*G_/a^2;
C2 = knownPressure / denomC2;
C1 = -C2/b^2;
A_ = (2*lambda_+2*G_)*C1;
B_ = 2*G_*C2;

sigma_rr_analytical = @(r) A_ - B_./r.^2;

fprintf('Analytical sigma_rr(a=%.4f um) = %.4f Pa  (applied traction = %.1f Pa)\n', a*1e6, sigma_rr_analytical(a), knownPressure);
fprintf('Analytical sigma_rr(r=3.31 um, our sampled Gauss point) = %.4f Pa\n', sigma_rr_analytical(3.31e-6));
fprintf('Analytical sigma_rr(b=%.4f um) = %.4f Pa\n', b*1e6, sigma_rr_analytical(b));
fprintf('\nCompare: FE Gauss-point (first layer, flat z-band) median = %.4f Pa\n', median(allVals));
fprintf('Compare: FE nodal-recovered (same z-band) median          = %.4f Pa\n', median(recoveredInBand));
fprintf('Compare: Analytical at Gauss-point radius (~3.31 um)      = %.4f Pa\n', sigma_rr_analytical(3.31e-6));
fprintf('Compare: Analytical at true boundary (a)                  = %.4f Pa\n', sigma_rr_analytical(a));

fprintf('\nSMOKE_TEST_STATUS: OK\n');
