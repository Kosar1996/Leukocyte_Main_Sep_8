clc;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
k = out.stopStep;
fl = out.fluidHist{k};
par = out.par;

meshF2 = add_fluid_nodes(fl.meshF);
mraw = fl.meshF;

% L1 query point and its element (2,60), per the corrected Step 4 table
rq = 3.1051e-6; zq = 3.4500e-6;
i = 2; j = 60;
e = (j-1)*(mraw.Nr-1) + i;
conn = meshF2.elems(e,:);
xe = meshF2.nodes(conn,1);
ze = meshF2.nodes(conn,2);

fprintf('L1 query point: r=%.4f, z=%.4f um\n', rq*1e6, zq*1e6);
fprintf('Element (2,60) corners:\n');
for c = 1:4
    fprintf('  Corner %d: r=%.4f, z=%.4f um  (dist from query = %.4f um)\n', ...
        c, xe(c)*1e6, ze(c)*1e6, norm([xe(c)-rq, ze(c)-zq])*1e6);
end

% Solve for true local (xi,eta) -- same Newton iteration as locate_and_interp_fluid_stress.m
xi = 0; eta = 0;
for it = 1:50
    [N, dNdxi] = shape_Q4(xi, eta);
    r_cur = N.' * xe; z_cur = N.' * ze;
    resid = [rq - r_cur; zq - z_cur];
    if norm(resid) < 1e-13, break; end
    J = [xe ze].' * dNdxi;
    dxieta = J \ resid;
    xi = xi + dxieta(1); eta = eta + dxieta(2);
end
fprintf('\nTrue local coordinates: xi=%.4f, eta=%.4f\n', xi, eta);

[N, ~] = shape_Q4(xi, eta);
fprintf('\nShape function weights (proximity-based, NEW method):\n');
for c = 1:4
    fprintf('  Corner %d weight: %.4f\n', c, N(c));
end

[N_old, ~] = shape_Q4(0, 0);
fprintf('\nShape function weights (OLD method, always centroid xi=eta=0):\n');
for c = 1:4
    fprintf('  Corner %d weight: %.4f (equal, regardless of actual position)\n', c, N_old(c));
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
