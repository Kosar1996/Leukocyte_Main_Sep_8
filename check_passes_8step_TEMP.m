S = load('out_short8_after_dtfix.mat');
out = S.out;
fprintf('%-6s %-12s %-12s %-14s\n', 'step', 'passesUsed', 'converged', 'min_gap[m]');
for k = 1:out.stopStep
    st = out.stateHist{k};
    gapMin = min(st.deltaE(:) - st.deltaL(:));
    if isfield(st, 'tractionCorrectionPassesUsed')
        fprintf('%-6d %-12d %-12d %-14.4e\n', k, st.tractionCorrectionPassesUsed, st.tractionCorrectionConverged, gapMin);
    else
        fprintf('%-6d %-12s %-12s %-14.4e\n', k, 'N/A', 'N/A', gapMin);
    end
end
