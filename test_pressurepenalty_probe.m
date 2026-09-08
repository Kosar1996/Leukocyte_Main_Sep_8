%% TEST_PRESSUREPENALTY_PROBE

% pass. solve_fluid_2D_bodyfitted_MAC.m was found to unconditionally zero
% par.pressurePenalty before calling solve_stokes_bodyfitted_MAC.m, even
% though that solver already supports an artificial-compressibility
% pressure-stabilization term (-pressurePenalty*p added to each continuity
% row) meant to help the matrix conditioning problem as the gap narrows.
% That hardcoding was just changed to be overridable via
% par.fluidPressurePenalty (default 0, unchanged behavior).
%
% This script is a QUICK single-fine-step probe (not the usual 4-step
% test) run at one specific pressurePenalty value, set by PP_VALUE before
% this script is called. Purpose: sanity-check that a given magnitude (a)
% doesn't break mass conservation (dV/V, globalMass, fluxJump, sourceInt
% diagnostics) and (b) actually changes solver behavior at all, before
% committing to a full multi-step isolated comparison run. Uses the
% ORIGINAL (buggy rOuter=4e-6) starting state, unchanged, consistent with
% every other isolated test tonight. Writes to its own output/log files,
% named by PP_VALUE, only.

clc;

if ~exist('PP_VALUE', 'var')
    error('PP_VALUE must be set before running this script.');
end

S = load('out_pure2dmac_dtlarge_7steps.mat');
out = S.out;

fprintf('Run B step 1: t=%.6e s, dt used=%.6e s\n', out.t(1), out.dtHist(1));

state1 = out.stateHist{1};
par = out.par;
par.dt = 3e-4;
par.useUnsteadyStokes = true;
par.rho = 1000;
nFineSteps = 1;
par.tEnd = state1.t + nFineSteps * par.dt;

par.fluidPressurePenalty = PP_VALUE;
fprintf('PROBE: par.fluidPressurePenalty = %g\n', par.fluidPressurePenalty);

restart = struct();
restart.state = state1;
restart.par = par;
restart.meshE = out.meshE;
restart.interfaceE = out.interfaceE;
restart.meshL = out.meshL;
restart.interfaceL = out.interfaceL;
restart.z = out.z;
restart.PHist = out.PHist;
restart.urCHist = out.urCHist;
restart.uzCHist = out.uzCHist;
restart.speedCHist = out.speedCHist;
restart.RPHist = out.RPHist;
restart.ZPHist = out.ZPHist;

cfg = struct();
cfg.geometry.endotheliumPrestressFile = 'solid_endo_P300_wide.mat';
cfg.geometry.leukocytePrestressFile = 'solid_leu_P600.mat';

outFine = softlube_run_case_global_coupled(cfg, restart);

outFileName = sprintf('out_pressurepenalty_probe_%g.mat', PP_VALUE);
save(outFileName, 'outFine', '-v7.3');
fprintf('\nSaved full output to %s\n', outFileName);

fprintf('\n=== RESULT (pressurePenalty=%g probe) ===\n', PP_VALUE);
for k = 1:outFine.stopStep
    if k <= numel(outFine.fluidHist) && isstruct(outFine.fluidHist{k}) && isfield(outFine.fluidHist{k}, 'p')
        fprintf('Fine step %d: t=%.6e s, max|p|=%.4f Pa, min gap=%.6e m\n', ...
            k, outFine.t(k), max(abs(outFine.fluidHist{k}.p)), min(outFine.fluidHist{k}.gap(:)));
    else
        fprintf('Fine step %d: t=%.6e s, INVALID fluidHist entry (class=%s)\n', ...
            k, outFine.t(k), class(outFine.fluidHist{k}));
    end
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
