%% PLOT_FARFIELD_BOTH
% Plots the full-R-range (0 to the endothelium's true outer radius, not
% truncated at 4 um) stress components for both the hybrid (1D) and full
% 2D saved comparison runs, at their final step, so the far-field

% request. Uses the already-saved data from run_t1_vs_t10_comparison.m,
% no rerun needed.

clc; close all;
cd(fileparts(mfilename('fullpath')));

files = {'out_1D_t10_for_review.mat', 'out_2D_t10_for_review.mat'};
labels = {'1D (hybrid)', '2D (full)'};

for i = 1:2
    S = load(files{i});
    fn = fieldnames(S);
    out = S.(fn{1});

    plotstep = out.stopStep;
    fprintf('%s: stopStep=%d, t=%.4e s\n', files{i}, out.stopStep, out.t(plotstep));

    % Same fluid-node/pressure-stress recovery step read_and_plot.m does
    out.fluidHist{plotstep}.meshF = add_fluid_nodes(out.fluidHist{plotstep}.meshF);
    [pCell, sigmaCell, center] = recover_fluid_nodes_pressure_stress_Q4( ...
        out.fluidHist{plotstep}.meshF, out.fluidHist{plotstep}.ur2D, ...
        out.fluidHist{plotstep}.uz2D, out.par.mu, out.fluidHist{plotstep}.pCell);
    out.fluidHist{plotstep}.pCellNode = pCell;
    out.fluidHist{plotstep}.centerNode = center;
    out.fluidHist{plotstep}.sigmaCellNode = sigmaCell;

    plot_select_native2d_stress(out, plotstep);
    sgtitle(sprintf('%s, step %d (full R range)', labels{i}, plotstep));
end

fprintf('\nDone -- compare the two figures: does stress decay to zero by the true far edge in either one?\n');
