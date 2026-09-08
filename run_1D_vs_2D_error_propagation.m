%% RUN_1D_VS_2D_ERROR_PROPAGATION
% Runs the same leukocyte/endothelium case two ways:
%   "1D version" -- the gap is solved with the 1D lubrication p(z)
%                   formula, only the exterior reservoir is true 2D
%                   Stokes/MAC (this is what run_input_full2D_pressure2.m
%                   / build_cfg_full2D_pressure2.m does by default).
%   "2D version" -- the entire domain, including the gap, is solved as
%                   true 2D Stokes/MAC. No 1D shortcut, no gap/exterior
%                   stitching seam.
%
% For each version: runs nSteps timesteps, saves the full "out" struct
% to its own .mat file (so it can be reloaded later with
% read_and_plot.m-style scripts instead of rerunning the solve), and
% checks the pressure-vs-1/3-trace(sigma) error over time using
% check_pressure_stress_consistency_over_time.m.
%
% Requires build_cfg_full2D_pressure2.m,
% check_pressure_stress_consistency_over_time.m, and
% recover_fluid_nodes_pressure_stress_Q4.m (with the pRaw output added)
% on the MATLAB path.

clc;
clearvars;
close all;

nSteps = 40;   % number of timesteps for each version. Increase once you
              % know the 2D version's per-step runtime -- it may be
              % noticeably slower than the hybrid 1D version, since the
              % gap region is no longer shortcut by the analytic
              % lubrication profile.

results = struct();

%% --- 1D version (hybrid gap-1D / exterior-2D) ---
fprintf('\n================ Running 1D (hybrid) version ================\n');
cfg1D = build_cfg_full2D_pressure2('nSteps', nSteps, ...
    'useHybridGap1DExterior2DFluid', true);
out1D = softlube_run_case_global_coupled(cfg1D);
out1D.cfg = cfg1D;
save('out_1D_hybrid_version.mat', 'out1D');
fprintf('Saved 1D (hybrid) result to out_1D_hybrid_version.mat (stopStep=%d)\n', ...
    out1D.stopStep);

fprintf('\n--- Error check: 1D (hybrid) version ---\n');
summary1D = check_pressure_stress_consistency_over_time(out1D);
results.summary1D = summary1D;

%% --- 2D version (true 2D everywhere, no hybrid gap shortcut) ---
fprintf('\n================ Running 2D (full) version ================\n');
cfg2D = build_cfg_full2D_pressure2('nSteps', nSteps, ...
    'useHybridGap1DExterior2DFluid', false);
out2D = softlube_run_case_global_coupled(cfg2D);
out2D.cfg = cfg2D;
save('out_2D_full_version.mat', 'out2D');
fprintf('Saved 2D (full) result to out_2D_full_version.mat (stopStep=%d)\n', ...
    out2D.stopStep);

fprintf('\n--- Error check: 2D (full) version ---\n');
summary2D = check_pressure_stress_consistency_over_time(out2D);
results.summary2D = summary2D;

%% --- Side-by-side comparison ---
fprintf('\n============== 1D vs 2D error comparison ==============\n');
fprintf('%10s | %16s %16s | %16s %16s\n', ...
    't [s]', '1D diffP90[Pa]', '1D % of maxP', '2D diffP90[Pa]', '2D % of maxP');
nCompare = min(numel(summary1D.t), numel(summary2D.t));
for k = 1:nCompare
    fprintf('%10.4e | %16.6e %15.4f%% | %16.6e %15.4f%%\n', ...
        summary1D.t(k), summary1D.diffP90(k), 100*summary1D.diffRelToP(k), ...
        summary2D.diffP90(k), 100*summary2D.diffRelToP(k));
end

figure;
plot(summary1D.t*1e3, 100*summary1D.diffRelToP, 'o-', 'LineWidth', 2); hold on;
plot(summary2D.t*1e3, 100*summary2D.diffRelToP, 's--', 'LineWidth', 2);
xlabel('t [ms]');
ylabel('mismatch as % of max pressure');
legend({'1D (hybrid gap)', '2D (full, no hybrid)'}, 'Location', 'best');
title('Pressure/stress-trace mismatch: 1D vs 2D solver');
grid on;

save('results_1D_vs_2D_comparison.mat', 'results', 'summary1D', 'summary2D');
fprintf('\nSaved comparison summary to results_1D_vs_2D_comparison.mat\n');
fprintf(['\nTo re-analyze either run later without re-solving, load its .mat file ', ...
    '(out_1D_hybrid_version.mat or out_2D_full_version.mat) the same way ', ...
    'read_and_plot.m does, e.g.:\n', ...
    '  load(''out_2D_full_version.mat''); out = out2D;\n', ...
    '  plot_select_native2d_stress(out, 1);\n']);
