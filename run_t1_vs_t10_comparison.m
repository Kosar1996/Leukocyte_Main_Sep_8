%% RUN_T1_VS_T10_COMPARISON

% wall-fix, no other patches -- just plot_select_native2d_stress.m as it
% stands) for the hybrid gap-1D/exterior-2D ("1D") version AND the full
% 2D ("2D") version, both for 10 time steps, then plot the 4-panel
% stress-component heat maps at t=1 (plotstep=1) and t=10 (plotstep=10)
% for each, so you can visually compare:
%   - whether both versions look physical
%   - whether the same issue is present at step 1 as at step 10 (i.e. is
%     it there from the very beginning, or does it grow/appear later)
%   - whether it's more pronounced in 1D or 2D
%
% Produces 4 figures total:
%   1D, t = step 1
%   1D, t = step 10
%   2D, t = step 1
%   2D, t = step 10
%
% Requires build_cfg_full2D_pressure2.m, add_native2d_field.m, and
% plot_select_native2d_stress.m on the path.

clc;
clearvars;
close all;

nSteps = 40;

%% --- 1D version (hybrid gap-1D / exterior-2D) ---
fprintf('\n================ Running 1D (hybrid) version, %d steps ================\n', nSteps);
cfg1D = build_cfg_full2D_pressure2('nSteps', nSteps, ...
    'useHybridGap1DExterior2DFluid', true);
out1D = softlube_run_case_global_coupled(cfg1D);
out1D = add_native2d_field(out1D);
save('out_1D_t10_for_review.mat', 'out1D');
fprintf('Saved 1D result (stopStep=%d) to out_1D_t10_for_review.mat\n', out1D.stopStep);

%% --- 2D version (full 2D everywhere, no hybrid gap shortcut) ---
fprintf('\n================ Running 2D (full) version, %d steps ================\n', nSteps);
cfg2D = build_cfg_full2D_pressure2('nSteps', nSteps, ...
    'useHybridGap1DExterior2DFluid', false);
out2D = softlube_run_case_global_coupled(cfg2D);
out2D = add_native2d_field(out2D);
save('out_2D_t10_for_review.mat', 'out2D');
fprintf('Saved 2D result (stopStep=%d) to out_2D_t10_for_review.mat\n', out2D.stopStep);

%% --- Four comparison figures ---
fprintf('\nPlotting 1D, t = step 1 ...\n');
plot_select_native2d_stress(out1D, 1);
sgtitle('1D (hybrid), step 1');

fprintf('Plotting 1D, t = step %d ...\n', out1D.stopStep);
plot_select_native2d_stress(out1D, out1D.stopStep);
sgtitle(sprintf('1D (hybrid), step %d', out1D.stopStep));

fprintf('Plotting 2D, t = step 1 ...\n');
plot_select_native2d_stress(out2D, 1);
sgtitle('2D (full, no hybrid), step 1');

fprintf('Plotting 2D, t = step %d ...\n', out2D.stopStep);
plot_select_native2d_stress(out2D, out2D.stopStep);
sgtitle(sprintf('2D (full, no hybrid), step %d', out2D.stopStep));

fprintf(['\nDone. Compare the four figures: same issue at step 1 as step %d? ', ...
    'Worse in 1D or 2D? Does either look physically implausible?\n'], out2D.stopStep);
