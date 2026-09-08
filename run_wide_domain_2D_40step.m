%% RUN_WIDE_DOMAIN_2D_40STEP
% Item 6 domain-extension test: same 40-step 2D production run as
% out_2D_t10_for_review.mat, but using the wide-domain endothelium
% prestress (solid_endo_P300_wide.mat, z=[-2,6]um matching the
% leukocyte's own mesh span) instead of the original narrow one
% (solid_endo_P300.mat, z=[0,4]um). Same everything else -- same load,
% same solver settings, same time steps -- so any difference in the
% contact-zone stress/mismatch numbers isolates the effect of the free
% z-boundary's distance from the contact zone.

clc; close all;
cd(fileparts(mfilename('fullpath')));

cfg = build_cfg_full2D_pressure2('nSteps', 40, 'useHybridGap1DExterior2DFluid', false);
cfg.geometry.endotheliumPrestressFile = fullfile(fileparts(mfilename('fullpath')), 'solid_endo_P300_wide.mat');

fprintf('Running 2D, 40 steps, WIDE endothelium domain (z=[-2,6]um) ...\n');
out = softlube_run_case_global_coupled(cfg);

save('out_2D_wide_domain_40step.mat', 'out', '-v7.3');
fprintf('Saved out_2D_wide_domain_40step.mat, stopStep=%d\n', out.stopStep);
fprintf('SMOKE_TEST_STATUS: OK\n');
