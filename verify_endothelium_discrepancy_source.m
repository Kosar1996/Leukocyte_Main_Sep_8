clc;
%% Check 1: shape_Q4(0,0) vs q4_shape(0,0,1.0) -- are they identical?
[N1, dN1] = shape_Q4(0,0);
[N2, dN2, ~] = q4_shape(0,0,1.0);
fprintf('shape_Q4 N:  %s\n', mat2str(N1(:).',6));
fprintf('q4_shape N:  %s\n', mat2str(N2(:).',6));
fprintf('shape_Q4 dNdxi:\n'); disp(dN1);
fprintf('q4_shape dNdxi:\n'); disp(dN2);
fprintf('Max abs diff N: %.2e, dNdxi: %.2e\n\n', max(abs(N1(:)-N2(:))), max(abs(dN1(:)-dN2(:))));

%% Check 2: for endothelium element #585's node, what elements share that node,
%% and how much does each one's raw (single-element) zz/rz stress vary?
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
k = out.stopStep;
st = out.stateHist{k};
par = out.par;
dtStep = out.dtHist(k);

mesh = out.meshE;
elemIdx = 585;
targetNode = mesh.conn(elemIdx, 1); % pick corner-1 node of this element
touchingElems = find(any(mesh.conn == targetNode, 2));
fprintf('Node %d (corner-1 of element %d) is touched by elements: %s\n', targetNode, elemIdx, mat2str(touchingElems.'));

for ei = 1:numel(touchingElems)
    e = touchingElems(ei);
    conn = mesh.conn(e,:);
    Xe = mesh.nodes(conn,:);
    ue = zeros(8,1); ueOld = zeros(8,1);
    for a=1:4
        ue(2*a-1) = st.uE(2*conn(a)-1); ue(2*a) = st.uE(2*conn(a));
        ueOld(2*a-1) = st.uEPrev(2*conn(a)-1); ueOld(2*a) = st.uEPrev(2*conn(a));
    end
    [N, dNdxi] = shape_Q4(0,0);
    Rnod = Xe(:,1); Znod = Xe(:,2);
    rnod = Rnod + ue(1:2:end); znod = Znod + ue(2:2:end);
    rnodOld = Rnod + ueOld(1:2:end); znodOld = Znod + ueOld(2:2:end);
    Jmat = Xe.' * dNdxi;
    dNdX = dNdxi / Jmat;
    Rg = N.' * Rnod; rg = N.' * rnod;
    drdR = dNdX(:,1).'*rnod; drdZ = dNdX(:,2).'*rnod;
    dzdR = dNdX(:,1).'*znod; dzdZ = dNdX(:,2).'*znod;
    F = [drdR,0,drdZ; 0,rg/Rg,0; dzdR,0,dzdZ];
    rgOld = N.'*rnodOld;
    drdROld = dNdX(:,1).'*rnodOld; drdZOld = dNdX(:,2).'*rnodOld;
    dzdROld = dNdX(:,1).'*znodOld; dzdZOld = dNdX(:,2).'*znodOld;
    Fold = [drdROld,0,drdZOld; 0,rgOld/Rg,0; dzdROld,0,dzdZOld];
    Jdet = det(F); B = F*F.'; I3=eye(3); devB = B-(trace(B)/3)*I3;
    Te = par.Ge*Jdet^(-5/3)*devB + par.Ke*(Jdet-1)*I3;
    parVisc = par; parVisc.dt = dtStep;
    [~,kv] = objective_kelvin_voigt_piola(F,Fold,parVisc);
    T = Te + kv.Tvisc;
    fprintf('  elem %4d: zz=%.5f rz=%.5f\n', e, T(3,3), T(1,3));
end
fprintf('\nSMOKE_TEST_STATUS: OK\n');
