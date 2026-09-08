clc;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg')
    cfg = out.cfg;
    if isfield(cfg, 'geometry') && isfield(cfg.geometry, 'endotheliumPrestressFile')
        fprintf('endotheliumPrestressFile: %s\n', cfg.geometry.endotheliumPrestressFile);
    end
    fprintf('cfg saved: yes\n');
else
    fprintf('cfg field NOT present in this saved file.\n');
end
if isfield(out, 'par')
    par = out.par;
    fields_to_check = {'fineDz','coarseDz','NrExterior2D','Nr','NrFluid2D','useHybridGap1DExterior2DFluid','useBodyFittedMACFluid'};
    for i = 1:numel(fields_to_check)
        fn2 = fields_to_check{i};
        if isfield(par, fn2)
            v = par.(fn2);
            if isnumeric(v) && isscalar(v)
                fprintf('par.%s = %g\n', fn2, v);
            elseif islogical(v)
                fprintf('par.%s = %d\n', fn2, v);
            end
        end
    end
end
fprintf('stopStep=%d\n', out.stopStep);
if isfield(out, 'stateHist')
    fprintf('numel(stateHist)=%d\n', numel(out.stateHist));
end
fprintf('\nSMOKE_TEST_STATUS: OK\n');
