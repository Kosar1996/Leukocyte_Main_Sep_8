%% DIAGNOSE_RADIAL_AXIAL_TRACTION_CORRECTION
% The real-time-step test (endothelium_normal_vs_tangential_over_time.m)
% ruled out "sigma_rr is intrinsically frozen by prestress" -- sigma_rr
% swings ~2x over real accepted steps, so it's clearly not physically inert.
% That means the freeze seen in EnE during the traction-CORRECTION sub-loop
% is something specific to the correction mechanism itself.
%
% This turns on a debug print (opts.debugRadialAxial, wired through
% par.debugRadialAxialTractionCorrection) inside
% apply_bodyfitted_MAC_traction_correction_feedback.m that reports, per
% correction pass, the norm of the endothelium's RADIAL vs AXIAL
% displacement-update DOFs (uEcorr-uEold, before under-relaxation), plus
% the recovered sigma_rr/sigma_rz right after that pass's solid solve.
%
% If |dRadial| goes to ~0 while |dAxial| keeps moving, the correction's
% radial-direction update itself is stuck (a real bug in how the
% correction applies to the endothelium's normal direction). If both
% keep moving but sigma_rr still doesn't move, the issue is downstream
% (stress recovery / how EnE's traction is measured), not the solid solve.

clc; close all;
cd(fileparts(mfilename('fullpath')));

cfg = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', true);

testRelax = 0.3;
testMaxIter = 15;
cfg.fluid.bodyFittedTractionCorrectionRelax = testRelax;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = testRelax;
cfg.fluid.maxBodyFittedTractionCorrections = testMaxIter;
cfg.parOverrides.maxBodyFittedTractionCorrections = testMaxIter;
cfg.fluid.debugRadialAxialTractionCorrection = true;
cfg.parOverrides.debugRadialAxialTractionCorrection = true;

fprintf('Running relax=%.2f, maxIterations=%d, with radial/axial debug on ...\n', testRelax, testMaxIter);
out = softlube_run_case_global_coupled(cfg);
fprintf('stopStep=%d\n', out.stopStep);
fprintf('SMOKE_TEST_STATUS: OK\n');
