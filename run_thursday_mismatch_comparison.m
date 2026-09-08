%% RUN_THURSDAY_MISMATCH_COMPARISON

% between the fluid and solid phases, for the first timestep and for a
% longer simulation time, using the established reporting format
% (check_interface_traction_mismatch_report.m).

clc;

S = load('out_pure2dmac_dtlarge_7steps.mat');
out = S.out;

opts = struct('makePlot', false);

opts.label = 'first timestep';
report1 = check_interface_traction_mismatch_report(out, 1, opts);

opts.label = sprintf('longer simulation (step %d)', out.stopStep);
reportLast = check_interface_traction_mismatch_report(out, out.stopStep, opts);

fprintf('\nSMOKE_TEST_STATUS: OK\n');
