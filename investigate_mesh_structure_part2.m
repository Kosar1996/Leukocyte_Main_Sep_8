clc;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
k = out.stopStep;
fl = out.fluidHist{k};
mraw = fl.meshF;

Zcol = mraw.Zp(1,:);
fprintf('Full z-domain: [%.4f, %.4f] um, Nz=%d\n', min(Zcol)*1e6, max(Zcol)*1e6, mraw.Nz);
fprintf('\nSpacing (diff) across full range, sampled every 10:\n');
dZ = diff(Zcol)*1e6;
for j = 1:10:numel(dZ)
    fprintf('  j=%d->%.d: dz=%.4f um\n', j, j+1, dZ(j));
end
fprintf('\nMin dz=%.4f um, Max dz=%.4f um\n', min(dZ), max(dZ));

fprintf('\n--- Checking element (2,60), the one reported for L1/L2 ---\n');
i = 2; j = 60;
Nr = mraw.Nr;
e = (j-1)*(Nr-1) + i;
fprintf('Element index e = %d\n', e);
cn = fl.centerNode;
fprintf('centerNode(%d,:) = r=%.6f z=%.6f um\n', e, cn(e,1)*1e6, cn(e,2)*1e6);

fprintf('\nNeighboring elements centerNode z-values (same i=2, j=58..62):\n');
for jj = 58:62
    ee = (jj-1)*(Nr-1) + i;
    fprintf('  j=%d (elem %d): centerNode z=%.6f um  |  node Zp(i,j)=%.6f, Zp(i,j+1)=%.6f\n', ...
        jj, ee, cn(ee,2)*1e6, mraw.Zp(i,jj)*1e6, mraw.Zp(i,jj+1)*1e6);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
