S = load('/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat');
fn = fieldnames(S); out = S.(fn{1});
fprintf('meshE z range: %.4f to %.4f um\n', min(out.meshE.nodes(:,2))*1e6, max(out.meshE.nodes(:,2))*1e6);
fprintf('meshL z range: %.4f to %.4f um\n', min(out.meshL.nodes(:,2))*1e6, max(out.meshL.nodes(:,2))*1e6);
if isfield(out,'cfg')
    try
        disp(out.cfg.geometry.endotheliumPrestressFile);
    catch
        fprintf('cfg present but no endotheliumPrestressFile field directly\n');
    end
end
if isfield(out,'par') && isfield(out.par,'endotheliumPrestressFile')
    disp(out.par.endotheliumPrestressFile)
end
fprintf('stopStep = %d, total steps in dtHist = %d\n', out.stopStep, numel(out.dtHist));
