clc;
files = {'/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat', ...
         '/Users/kosarsafari/Desktop/Project_1/code/leukocyte-main/out_2D_unified_domain_40step.mat'};
labels = {'OLD (out_2D_t10_for_review)', 'NEW (out_2D_unified_domain_40step)'};

for fi = 1:2
    S = load(files{fi});
    fn = fieldnames(S);
    out = S.(fn{1});
    if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
    k = out.stopStep;
    fl = out.fluidHist{k};
    mraw = fl.meshF;
    Zcol = mraw.Zp(1,:);
    dZ = diff(Zcol)*1e6;
    fprintf('%s: Nz=%d, z-range=[%.4f,%.4f] um\n', labels{fi}, mraw.Nz, min(Zcol)*1e6, max(Zcol)*1e6);
    fprintf('  dz min=%.4f max=%.4f (uniform? %d)\n', min(dZ), max(dZ), max(dZ)-min(dZ) < 1e-8);
    fprintf('  dz near z=3.4-3.6um:\n');
    idx = find(Zcol*1e6 > 3.0 & Zcol*1e6 < 4.0);
    for ii = idx
        fprintf('    j=%d z=%.4f\n', ii, Zcol(ii)*1e6);
    end
    fprintf('\n');
end
fprintf('SMOKE_TEST_STATUS: OK\n');
