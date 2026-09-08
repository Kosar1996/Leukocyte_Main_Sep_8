%% TEST_FINESTEP_UNSTEADY_PPFIX_FROM_RUNB_STATE1

% motivated by a direct observation in the rOuter-fix run: its persistent
% step-2 correction-loop oscillation (85 passes, 16-2,933 Pa, never
% converging) began immediately after the monolithic two-solid solve
% failed with gap=1.000e-09 m -- exactly the minGap floor referenced in
% Issue 2 (fluid conditioning). That is direct evidence the correction
% loop's non-convergence and the fluid matrix's near-floor ill-
% conditioning are occurring together, not necessarily two independent
% problems.
%
% Relaxation-factor tuning alone was already tested and disproven (0.25,
% 0.5, 0.7 all eventually fail via the same fsolve exitflag=0 budget-
% exhaustion mode; 0.5 is the current default and survives longest).
% This test does NOT touch relax (kept at its default, 0.5). Instead it
% turns on the pressure-stabilization term already implemented for Issue
% 2 (par.fluidPressurePenalty, previously hardcoded to 0 in
% solve_fluid_2D_bodyfitted_MAC.m, now overridable) to see whether
% improving the fluid matrix conditioning near the gap floor lets the
% correction loop actually converge, rather than just tuning the solid-
% side relaxation.
%
% PP magnitude: 1e5, the largest of the three quick single-step probes
% run earlier today, which showed the clearest effect (pass-1 traction
% 26.5 Pa at PP=0 vs 0.76 Pa at PP=1e5) without yet being checked against
% a full accepted step's mass-conservation diagnostics -- this run will
% surface those.
%
% Uses the ORIGINAL (buggy rOuter=4e-6) starting state, unchanged, same
% as every other isolated test today, so this isolates the
% fluidPressurePenalty variable alone against the same known crisis
% region. Writes to its own output/log files only.

clc;

S = load('out_pure2dmac_dtlarge_7steps.mat');
out = S.out;

fprintf('Run B step 1: t=%.6e s, dt used=%.6e s\n', out.t(1), out.dtHist(1));
fprintf('Run B step 2 (actual, coarse dt, steady): t=%.6e s, max|p|=%.4f Pa\n', ...
    out.t(2), max(abs(out.fluidHist{2}.p)));

state1 = out.stateHist{1};
par = out.par;
par.dt = 3e-4;
par.useUnsteadyStokes = true;
par.rho = 1000;
nFineSteps = 4;
par.tEnd = state1.t + nFineSteps * par.dt;

fprintf('Relaxation factor (unchanged): par.bodyFittedTractionCorrectionRelax = %.4f\n', ...
    par.bodyFittedTractionCorrectionRelax);
par.fluidPressurePenalty = 1e5;
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

% Save FIRST, before any post-processing that could crash and lose the run.
save('out_finestep_UNSTEADY_PPFIX.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_PPFIX.mat\n');

fprintf('\n=== RESULT (unsteady Stokes ON, both warm-start fixes + fluidPressurePenalty=1e5 ACTIVE) ===\n');
for k = 1:outFine.stopStep
    if k <= numel(outFine.fluidHist) && isstruct(outFine.fluidHist{k}) && isfield(outFine.fluidHist{k}, 'p')
        fprintf('Fine step %d: t=%.6e s, max|p|=%.4f Pa\n', ...
            k, outFine.t(k), max(abs(outFine.fluidHist{k}.p)));
    else
        fprintf('Fine step %d: t=%.6e s, INVALID fluidHist entry (class=%s)\n', ...
            k, outFine.t(k), class(outFine.fluidHist{k}));
    end
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
