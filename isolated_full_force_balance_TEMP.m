%% Global force balance check, done CORRECTLY this time: use the FULL
% internal force (elastic + Kelvin-Voigt viscous), matching exactly what
% solve_finite_def_solid's own Newton residual uses. The earlier version
% of this check (isolated_global_force_balance_TEMP.m) used
% assemble_finite_def_axisym alone, missing Fvisc entirely -- fine when
% the old false-convergence bug meant displacement (and therefore strain
% rate / Fvisc) barely changed from uOld, but not fine now that the
% solve genuinely deforms far from uOld within one step and the split
% check (isolated_elastic_vs_visc_split_TEMP.m) showed Fvisc is NOT
% negligible (implied viscous stress contribution ~559 Pa vs applied
% -500 Pa). This version sums the TRUE total Fint (elastic+viscous) so
% the balance check is apples-to-apples with what Newton actually solved.

clc; clear all;

cfg = struct();
cfg.geometry.endotheliumPrestressFile = fullfile('', 'solid_endo_P300_wide.mat');
cfg.geometry.leukocytePrestressFile = fullfile('', 'solid_leu_P600.mat');
cfg.geometry.usePrestressedLeukocyteIC = true;
cfg.geometry.useLeukocyteCenterline = true;

cfg.fluid.zMin = -6e-6; cfg.fluid.zMax = 10e-6;
cfg.fluid.mu = 1.2e-3; cfg.fluid.slipE = 0e-9; cfg.fluid.slipL = 0e-9;
cfg.fluid.NrFluid2D = 20; cfg.fluid.Nr = 20; cfg.fluid.NrExterior2D = 20;
cfg.fluid.bodyFittedRadialFineWindow = [2e-6,4e-6];
cfg.fluid.bodyFittedRadialFineWeight = 4.0;
cfg.fluid.innerBoundary = 'deformed_leukocyte';
cfg.fluid.useSlopeAwareWallKinematics = true;

cfg.bc.pIn = 0; cfg.bc.pOut = 0; cfg.bc.UwL = 0; cfg.bc.UwE = 0;

cfg.solid.noLeukocyte = false;
cfg.solid.leukocyte.enabled = true;
cfg.solid.leukocyte.deformable = true;
cfg.solid.leukocyte.EL = 200; cfg.solid.leukocyte.nuL = 0.46;
cfg.solid.leukocyte.useViscoelasticLeukocyte = true;
cfg.solid.leukocyte.etaL = 1; cfg.solid.leukocyte.supportL = 'axis';

cfg.solid.endothelium.Ee = 500; cfg.solid.endothelium.nuE = 0.46;
cfg.solid.endothelium.useViscoelasticEndothelium = true;
cfg.solid.endothelium.etaE = 1;

cfg.interface.useExactDeformedInterface = true;
cfg.interface.useExactDeformedInterfaceInMonolithic = true;
cfg.interface.useReferenceZForTractionMapping = false;
cfg.interface.zeroLeukocyteTractionOutsideOverlap = false;
cfg.interface.leukocytePressureSupportZ = [-6e-6, 10e-6];
cfg.interface.warnPrestressLoadMismatch = false;

dt = 3e-4;
cfg.numerics.dt = dt; cfg.numerics.tEnd = dt;
cfg.numerics.minGap = 0.001e-6;
cfg.numerics.maxNewtonMono = 100; cfg.numerics.tolNewtonMono = 1e-6;
cfg.numerics.monoSolidAbsTol = 1e-7; cfg.numerics.monoFluidAbsTol = 1e-8;
cfg.numerics.maxFunctionEvaluationsTwoSolid = 400;
cfg.numerics.enablePressureJumpRetry = false;
cfg.numerics.maxPressureJumpAbs = inf; cfg.numerics.maxPressureJump2DAbs = inf;
cfg.numerics.useObjectiveKelvinVoigt = true;
cfg.numerics.useDirectFluidSolve = true; cfg.numerics.useFsolveMono = true;
cfg.numerics.enableAdaptiveTimeStep = true; cfg.numerics.dtMin = dt/1024;
cfg.numerics.dtRetryFactor = 0.5; cfg.numerics.dtGrowFactor = 1.25;
cfg.numerics.maxTimeStepRetries = 18; cfg.numerics.fsolveDisplay = 'off';
cfg.numerics.diagnosticsEnabled = true;

cfg.output.saveOutput = false;
cfg.output.makePlots = false;
cfg.output.plotFinalGlobalPressureContour = false;
cfg.output.plotGlobal2DPressureTractionComparison = false;
cfg.output.printEvery = 1;
cfg.output.storeFull2DFluidHist = false;
cfg.output.store2DFluidArrays = true;
cfg.output.store2DPhysicalGridHist = true;
cfg.output.store2DSpeedHist = true;

cfg.runtime.useEnvironment = false;
cfg.ui.closeFigures = true;
cfg.parOverrides = struct();
cfg.global1D.zDomain = [cfg.fluid.zMin, cfg.fluid.zMax];
cfg.global1D.rDomain = [0, 4e-6];
cfg.global1D.axisPressureBC = 'dPdr=0';
cfg.global1D.outerPressureBC = 0;
cfg.parOverrides.useGlobal1DPressure = true;
cfg.parOverrides.global1DOuterRadius = 4e-6;
cfg.parOverrides.global1DAxisRadius = 1e-9;
cfg.parOverrides.useGlobal1DCoarseEdgeMesh = true;
cfg.parOverrides.global1DFineWindow = [-0.2e-6,4.2e-6];
cfg.parOverrides.global1DCoarseDz = 1.0e-6;
cfg.parOverrides.global1DFineDz = 0.05e-6;
cfg.parOverrides.global1DPlotNr = 161;
cfg.fluid.useFull2DFluid = true;
cfg.fluid.useHybridGap1DExterior2DFluid = false;
cfg.fluid.hybridGapZ = [-0.2e-6,4.2e-6];
cfg.fluid.hybridExteriorMinCells = 4;
cfg.fluid.hybridExteriorFailMode = 'warn';
cfg.fluid.hybridUseExterior2DPressureVector = true;
cfg.fluid.useBodyFittedMACFluid = true;
cfg.fluid.useBodyFittedMACTractionCorrection = true;
cfg.fluid.useFeedbackTractionCorrection = true;
cfg.fluid.useAitkenTractionCorrectionRelax = true;
cfg.fluid.debugRadialAxialTractionCorrection = true;
cfg.fluid.maxBodyFittedTractionCorrections = 12;
cfg.fluid.debugVerbose = false;
cfg.fluid.bodyFittedTractionCorrectionRelax = 0.5;
cfg.fluid.bodyFittedTractionCorrectionFailMode = 'warn';
cfg.fluid.resolveFluidAfterTractionCorrection = true;
cfg.fluid.useBodyFittedMACTractionInSolid = true;
cfg.parOverrides.useFull2DFluid = true;
cfg.parOverrides.useHybridGap1DExterior2DFluid = false;
cfg.parOverrides.hybridGapZ = [-0.2e-6,4.2e-6];
cfg.parOverrides.NrExterior2D = 20;
cfg.parOverrides.bodyFittedRadialFineWindow = [2e-6,4e-6];
cfg.parOverrides.bodyFittedRadialFineWeight = 4.0;
cfg.parOverrides.hybridExteriorMinCells = 4;
cfg.parOverrides.hybridExteriorFailMode = 'warn';
cfg.parOverrides.hybridUseExterior2DPressureVector = true;
cfg.parOverrides.useBodyFittedMACFluid = true;
cfg.parOverrides.useBodyFittedMACTractionCorrection = true;
cfg.parOverrides.useFeedbackTractionCorrection = true;
cfg.parOverrides.useAitkenTractionCorrectionRelax = true;
cfg.parOverrides.debugRadialAxialTractionCorrection = true;
cfg.parOverrides.maxBodyFittedTractionCorrections = 12;
cfg.parOverrides.debugVerbose = false;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = 0.5;
cfg.parOverrides.bodyFittedTractionCorrectionFailMode = 'warn';
cfg.parOverrides.resolveFluidAfterTractionCorrection = true;
cfg.parOverrides.useBodyFittedMACTractionInSolid = true;
cfg.parOverrides.NrFluid2D = 20;
cfg.parOverrides.Nr = 20;
cfg.parOverrides.useGlobal2DPressureTraction = false;

[par, S, uE_pre] = softlube_prepare_case(cfg);
meshE = prepare_axisym_mesh_cache(S.meshE);
interfaceE = S.interfaceE;
baseE = S.baseE;

zInterface = meshE.nodes(interfaceE, 2);
zLo = min(zInterface); zHi = max(zInterface);
nZ = 21;
testZ = linspace(zLo, zHi, nZ).';

knownPressure = -500;
traction = struct();
traction.z = testZ;
traction.normal = knownPressure * ones(nZ,1);
traction.tangent = zeros(nZ,1);

parTest = par;
if isfield(parTest, 'solidAbsTol'), parTest = rmfield(parTest, 'solidAbsTol'); end
if isfield(parTest, 'solidFallbackAbsTol'), parTest = rmfield(parTest, 'solidFallbackAbsTol'); end

uEcorr = solve_finite_def_solid(meshE, uE_pre, traction, interfaceE, baseE, par.supportE, parTest, uE_pre);

ndof = size(meshE.nodes,1)*2;

FextTotal = zeros(ndof,1);
[FextTotal, ~] = apply_interface_traction(meshE, uEcorr, FextTotal, interfaceE, traction);
totalAppliedRadial = sum(FextTotal(2*interfaceE-1));

[FintElastic, ~] = assemble_finite_def_axisym(meshE, uEcorr, par);
[FintVisc, ~] = assemble_axisym_kelvin_voigt_viscous(meshE, uEcorr, uE_pre, par);
FintTotal = FintElastic + FintVisc;

totalReactionRadial_elasticOnly = sum(FintElastic(2*baseE-1));
totalReactionRadial_full = sum(FintTotal(2*baseE-1));

% Global (whole-mesh) balance over ALL dofs, elastic+viscous, radial comp.
sumFint_all_radial = sum(FintTotal(1:2:end));
sumFext_all_radial = sum(FextTotal(1:2:end));

fprintf('\n=== FULL (ELASTIC+VISCOUS) GLOBAL FORCE BALANCE CHECK ===\n');
fprintf('Total applied radial force at interface        = %.6e N\n', totalAppliedRadial);
fprintf('Total ELASTIC-ONLY reaction radial force at base = %.6e N  (ratio = %.4f)\n', ...
    totalReactionRadial_elasticOnly, totalReactionRadial_elasticOnly/totalAppliedRadial);
fprintf('Total FULL (elastic+visc) reaction radial force at base = %.6e N  (ratio = %.4f)\n', ...
    totalReactionRadial_full, totalReactionRadial_full/totalAppliedRadial);
fprintf('\nWhole-mesh (all dofs) radial sum, FULL Fint = %.6e N, Fext = %.6e N\n', ...
    sumFint_all_radial, sumFext_all_radial);
if abs(sumFext_all_radial) > 1e-20
    fprintf('Relative whole-mesh imbalance (full) = %.4f%%\n', ...
        100*abs(sumFint_all_radial - sumFext_all_radial)/abs(sumFext_all_radial));
end

sumFintElastic_all_radial = sum(FintElastic(1:2:end));
fprintf('Whole-mesh (all dofs) radial sum, ELASTIC-ONLY Fint = %.6e N\n', sumFintElastic_all_radial);
if abs(sumFext_all_radial) > 1e-20
    fprintf('Relative whole-mesh imbalance (elastic-only, WRONG comparison) = %.4f%%\n', ...
        100*abs(sumFintElastic_all_radial - sumFext_all_radial)/abs(sumFext_all_radial));
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
