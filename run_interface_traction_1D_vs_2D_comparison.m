%% RUN_INTERFACE_TRACTION_1D_VS_2D_COMPARISON

% disk (no re-solving needed):
%   - item 2: does the interface-traction mismatch appear in BOTH the 1D
%     hybrid and the full 2D solver, and is it worse in one of them?
%   - item 4: is the mismatch already present at the very first accepted
%     time step (t = step 1), or does it only grow in over many steps?

%     interface/BC bug, not something that only accumulates.)
%
% This calls check_interface_traction_mismatch_report.m (item 3: en(s),
% et(s) with max/median/90th-percentile stats) on four cases:
%   1D hybrid, step 1        1D hybrid, final step
%   2D full,   step 1        2D full,   final step
%
% By default it loads the .mat files already produced by
% run_t1_vs_t10_comparison.m (out_1D_t10_for_review.mat /
% out_2D_t10_for_review.mat). Edit the file names below if you want to
% point this at a different pair of saved runs (e.g. the 40-step files
% from run_1D_vs_2D_error_propagation.m), or set out1D/out2D yourself in
% the workspace before running this script with runFromWorkspace = true.

clc;
close all;

runFromWorkspace = false;   % set true and define out1D/out2D beforehand to skip loading

file1D = 'out_1D_t10_for_review.mat';
file2D = 'out_2D_t10_for_review.mat';

if ~runFromWorkspace
    fprintf('Loading %s ...\n', file1D);
    S1 = load(file1D);
    out1D = S1.out1D;
    fprintf('Loading %s ...\n', file2D);
    S2 = load(file2D);
    out2D = S2.out2D;
end

cases = struct('out', {}, 'step', {}, 'label', {});
cases(end+1) = struct('out', out1D, 'step', 1,             'label', '1D hybrid, step 1');
cases(end+1) = struct('out', out1D, 'step', out1D.stopStep, 'label', sprintf('1D hybrid, step %d', out1D.stopStep));
cases(end+1) = struct('out', out2D, 'step', 1,             'label', '2D full, step 1');
cases(end+1) = struct('out', out2D, 'step', out2D.stopStep, 'label', sprintf('2D full, step %d', out2D.stopStep));

reports = cell(numel(cases),1);
opts = struct('makePlot', true);
for k = 1:numel(cases)
    opts.label = cases(k).label;
    reports{k} = check_interface_traction_mismatch_report( ...
        cases(k).out, cases(k).step, opts);
end

%% Side-by-side summary table
fprintf('\n================ Interface traction mismatch: side-by-side ================\n');
fprintf('%-22s | %10s %10s %10s | %10s %10s %10s\n', ...
    'case', 'Lmax|en|', 'Lmed|en|', 'Lp90|en|', 'Lmax|et|', 'Lmed|et|', 'Lp90|et|');
for k = 1:numel(cases)
    r = reports{k}.leuko.stats;
    fprintf('%-22s | %10.4f %10.4f %10.4f | %10.4f %10.4f %10.4f\n', ...
        cases(k).label, r.maxAbsEn, r.medianAbsEn, r.p90AbsEn, ...
        r.maxAbsEt, r.medianAbsEt, r.p90AbsEt);
end
fprintf('\n%-22s | %10s %10s %10s | %10s %10s %10s\n', ...
    'case', 'Emax|en|', 'Emed|en|', 'Ep90|en|', 'Emax|et|', 'Emed|et|', 'Ep90|et|');
for k = 1:numel(cases)
    r = reports{k}.endo.stats;
    fprintf('%-22s | %10.4f %10.4f %10.4f | %10.4f %10.4f %10.4f\n', ...
        cases(k).label, r.maxAbsEn, r.medianAbsEn, r.p90AbsEn, ...
        r.maxAbsEt, r.medianAbsEt, r.p90AbsEt);
end



save('interface_traction_1D_vs_2D_comparison.mat', 'reports', 'cases');
fprintf('\nSaved reports to interface_traction_1D_vs_2D_comparison.mat\n');
