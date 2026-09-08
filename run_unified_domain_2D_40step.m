%% RUN_UNIFIED_DOMAIN_2D_40STEP
% The "unified domain" run for item #12: same 40-step production run as
% out_2D_t10_for_review.mat, but now on the wide endothelium domain
% (z=[-2,6], now the default in build_cfg_full2D_pressure2.m /
% run_input_full2D_pressure2.m), AND with the fixed, position-aware fluid
% interpolation patched into compute_interface_traction_mismatch.m (used
% internally by the traction-correction loop). This is the first real run
% where both fixes are actually active during the correction process
% itself, not just applied in post-hoc verification.

clc; close all;
cd(fileparts(mfilename('fullpath')));

cfg = build_cfg_full2D_pressure2('nSteps', 40, 'useHybridGap1DExterior2DFluid', false);
% endotheliumPrestressFile now defaults to solid_endo_P300_wide.mat -- not overridden here,
% confirming the "unified domain" default is actually in effect, not a one-off override.

fprintf('Running 2D, 40 steps, UNIFIED (wide, now-default) domain, WITH fixed interpolation ...\n');
fprintf('endotheliumPrestressFile = %s\n', cfg.geometry.endotheliumPrestressFile);
tic;
out = softlube_run_case_global_coupled(cfg);
fprintf('Elapsed: %.1f s\n', toc);

save('out_2D_unified_domain_40step.mat', 'out', '-v7.3');
fprintf('Saved out_2D_unified_domain_40step.mat, stopStep=%d\n', out.stopStep);
fprintf('SMOKE_TEST_STATUS: OK\n');
