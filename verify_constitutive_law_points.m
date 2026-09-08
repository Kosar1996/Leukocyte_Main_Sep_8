%% VERIFY_CONSTITUTIVE_LAW_POINTS
% For 3 representative elements each in leukocyte and endothelium, show
% the full calculation explicitly: deformation gradient F, J=det(F),
% B=F*F', devB, and the constitutive-law stress prediction
% T = Ge*J^(-5/3)*devB + Ke*(J-1)*I -- so it can be checked by hand.

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

k = out.stopStep;
st = out.stateHist{k};
par = out.par;

parLmismatch = par;
if isfield(par, 'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par, 'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par, 'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

fprintf('Constitutive equation: T = Ge * J^(-5/3) * devB + Ke * (J-1) * I,  B = F*F^T,  devB = B - trace(B)/3 * I\n');
fprintf('Leukocyte: Ge=%.4f Pa, Ke=%.4f Pa | Endothelium: Ge=%.4f Pa, Ke=%.4f Pa\n', ...
    parLmismatch.Ge, parLmismatch.Ke, par.Ge, par.Ke);

function report_element(label, mesh, u, parUse, elemIdx)
    conn = mesh.conn(elemIdx,:);
    Xe = mesh.nodes(conn,:);
    ue = zeros(8,1);
    for a=1:4
        ue(2*a-1) = u(2*conn(a)-1);
        ue(2*a)   = u(2*conn(a));
    end

    [N, dNdxi] = shape_Q4(0,0);  % element centroid, matches internal recovery convention
    Rnod = Xe(:,1); Znod = Xe(:,2);
    rnod = Rnod + ue(1:2:end); znod = Znod + ue(2:2:end);
    Jmat = Xe.' * dNdxi;
    dNdX = dNdxi / Jmat;

    Rg = N.' * Rnod; rg = N.' * rnod;
    drdR = dNdX(:,1).' * rnod;
    drdZ = dNdX(:,2).' * rnod;
    dzdR = dNdX(:,1).' * znod;
    dzdZ = dNdX(:,2).' * znod;

    F = [drdR, 0, drdZ; 0, rg/Rg, 0; dzdR, 0, dzdZ];
    Jdet = det(F);
    B = F*F.';
    I3 = eye(3);
    devB = B - (trace(B)/3)*I3;
    T = parUse.Ge * Jdet^(-5/3) * devB + parUse.Ke * (Jdet-1) * I3;

    fprintf('\n--- %s element %d (centroid r=%.4f, z=%.4f um) ---\n', label, elemIdx, Rg*1e6, N.'*Znod*1e6);
    fprintf('  F = [%.5f %.5f %.5f; %.5f %.5f %.5f; %.5f %.5f %.5f]\n', F(1,1),F(1,2),F(1,3),F(2,1),F(2,2),F(2,3),F(3,1),F(3,2),F(3,3));
    fprintf('  J = det(F) = %.6f\n', Jdet);
    fprintf('  B = F*F^T diag = [%.5f, %.5f, %.5f], off-diag B13=%.5f\n', B(1,1),B(2,2),B(3,3),B(1,3));
    fprintf('  devB diag = [%.5f, %.5f, %.5f], off-diag devB13=%.5f\n', devB(1,1),devB(2,2),devB(3,3),devB(1,3));
    fprintf('  T (predicted from constitutive law): sigma_rr=%.4f sigma_tt=%.4f sigma_zz=%.4f sigma_rz=%.4f Pa\n', ...
        T(1,1), T(2,2), T(3,3), T(1,3));
end

fprintf('\n================ LEUKOCYTE ================\n');
nelL = out.meshL.nelem;
for e = [round(nelL*0.25), round(nelL*0.5), round(nelL*0.75)]
    report_element('Leukocyte', out.meshL, st.uL, parLmismatch, e);
end

fprintf('\n================ ENDOTHELIUM ================\n');
nelE = out.meshE.nelem;
for e = [round(nelE*0.25), round(nelE*0.5), round(nelE*0.75)]
    report_element('Endothelium', out.meshE, st.uE, par, e);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
