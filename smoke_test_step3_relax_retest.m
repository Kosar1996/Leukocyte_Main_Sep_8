%% SMOKE_TEST_STEP3_RELAX_RETEST
% Now that the frozen-correction bug is fixed (solve_finite_def_solid was
% exiting on iteration 1 due to an inherited absolute tolerance too loose
% for this context -- see apply_bodyfitted_MAC_traction_correction_feedback.m),
% retest with a larger relax to see if it actually helps convergence, or
% if the stall we saw at relax=0.05 (and confirmed even with 30 passes)
% persists regardless of relax.

clc; close all;
cd(fileparts(mfilename('fullpath')));

cfg = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', true);

testRelax = 0.3;
testMaxIter = 30;
cfg.fluid.bodyFittedTractionCorrectionRelax = testRelax;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = testRelax;
cfg.fluid.maxBodyFittedTractionCorrections = testMaxIter;
cfg.parOverrides.maxBodyFittedTractionCorrections = testMaxIter;

fprintf('Testing relax=%.2f (was 0.05), maxIterations=%d\n', testRelax, testMaxIter);
out = softlube_run_case_global_coupled(cfg);
fprintf('stopStep=%d\n', out.stopStep);
fprintf('\nDone. Watch the "EnE" column (endothelium normal, the dominant mismatch) --\n');
fprintf('does it trend clearly toward 10%% now, or still stall like before?\n');
