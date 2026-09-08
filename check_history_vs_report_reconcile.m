S = load('/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat');
fn = fieldnames(S); out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

fprintf('Top-level fields:\n');
disp(fieldnames(out));

k = out.stopStep;
fprintf('\nstopStep = %d\n', k);

fl = out.fluidHist{k};
fprintf('\nfluidHist{k} fields:\n');
disp(fieldnames(fl));

if isfield(fl, 'history')
    h = fl.history;
    fprintf('\nfl.history found, %d entries, fields:\n', numel(h));
    disp(fieldnames(h));
    if isfield(h, 'pctEnL')
        n = numel(h);
        fprintf('\nLast few correction-pass entries (this run''s own internal tracking):\n');
        fprintf('%6s %8s %8s %8s %8s\n', 'iter', 'EnL', 'EtL', 'EnE', 'EtE');
        for i = max(1,n-5):n
            fprintf('%6d %8.2f %8.2f %8.2f %8.2f\n', i, h(i).pctEnL, h(i).pctEtL, h(i).pctEnE, h(i).pctEtE);
        end
        fprintf('\nMax EnL over all passes: %.2f%%\n', max([h.pctEnL]));
        fprintf('Max EnE over all passes: %.2f%%\n', max([h.pctEnE]));
        fprintf('Max EtL over all passes: %.2f%%\n', max([h.pctEtL]));
        fprintf('Max EtE over all passes: %.2f%%\n', max([h.pctEtE]));
    end
else
    fprintf('\nNo fl.history field. Checking st (stateHist) for similar tracking...\n');
    st = out.stateHist{k};
    disp(fieldnames(st));
end
