clc;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
k = out.stopStep;
fl = out.fluidHist{k};
meshF2 = add_fluid_nodes(fl.meshF);
mraw = fl.meshF;

elems = struct('L1',[2,60],'L2',[2,60],'E1',[18,60],'E2',[18,60]);
labels = fieldnames(elems);
for li = 1:numel(labels)
    lab = labels{li};
    ij = elems.(lab);
    i = ij(1); j = ij(2);
    e = (j-1)*(mraw.Nr-1) + i;
    conn = meshF2.elems(e,:);
    xe = meshF2.nodes(conn,1)*1e6;
    ze = meshF2.nodes(conn,2)*1e6;
    fprintf('%s element (%d,%d): Corner1=(%.4f,%.4f) Corner2=(%.4f,%.4f) Corner3=(%.4f,%.4f) Corner4=(%.4f,%.4f)\n', ...
        lab, i, j, xe(1),ze(1), xe(2),ze(2), xe(3),ze(3), xe(4),ze(4));
end
fprintf('\nSMOKE_TEST_STATUS: OK\n');
