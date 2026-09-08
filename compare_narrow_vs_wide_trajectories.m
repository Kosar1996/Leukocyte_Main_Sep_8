clc;
Sorig = load('/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat');
fnO = fieldnames(Sorig); outO = Sorig.(fnO{1});
if isfield(outO,'cfg'), outO = rmfield(outO,'cfg'); end

Swide = load('out_2D_wide_domain_40step.mat');
outW = Swide.out;

nO = outO.stopStep; nW = outW.stopStep;
fprintf('orig stopStep=%d, wide stopStep=%d\n\n', nO, nW);

fprintf('%6s | %12s %12s | %10s %10s | %10s %10s\n', 'step', 'gap_orig[nm]', 'gap_wide[nm]', 'EnL_o', 'EnL_w', 'EnE_o', 'EnE_w');
for k = 1:min(nO,nW)
    dO = outO.diagHist{k}; dW = outW.diagHist{k};
    gapO = NaN; gapW = NaN;
    if ~isempty(dO) && isfield(dO,'gapMin'), gapO = dO.gapMin; end
    if ~isempty(dW) && isfield(dW,'gapMin'), gapW = dW.gapMin; end

    cO = outO.tractionCorrectionHistory{k}; cW = outW.tractionCorrectionHistory{k};
    EnLo=NaN; EnEo=NaN; EnLw=NaN; EnEw=NaN;
    if ~isempty(cO) && ~isempty(cO.history)
        EnLo = cO.history(end).pctEnL; EnEo = cO.history(end).pctEnE;
    end
    if ~isempty(cW) && ~isempty(cW.history)
        EnLw = cW.history(end).pctEnL; EnEw = cW.history(end).pctEnE;
    end
    fprintf('%6d | %12.2f %12.2f | %10.2f %10.2f | %10.2f %10.2f\n', ...
        k, gapO*1e9, gapW*1e9, EnLo, EnLw, EnEo, EnEw);
end
fprintf('\nSMOKE_TEST_STATUS: OK\n');
