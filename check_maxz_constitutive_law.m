clc;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
k = out.stopStep;
st = out.stateHist{k};
par = out.par;
dtStep = out.dtHist(k);

parLmismatch = par;
if isfield(par, 'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par, 'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par, 'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, par);
stressL = recover_nodal_stress_axisym_viscoelastic(out.meshL, st.uL, st.uLPrev, dtStep, parLmismatch);

function [elemIdx, zmax] = find_maxz_element(mesh, u)
nel = mesh.nelem;
zmaxPerElem = zeros(nel,1);
for e = 1:nel
    conn = mesh.conn(e,:);
    znod = mesh.nodes(conn,2) + u(2*conn);
    zmaxPerElem(e) = mean(znod);
end
[zmax, elemIdx] = max(zmaxPerElem);
end

function report_element(label, mesh, u, uOld, dt, parUse, elemIdx, solverStress)
    conn = mesh.conn(elemIdx,:);
    Xe = mesh.nodes(conn,:);
    ue = zeros(8,1); ueOld = zeros(8,1);
    for a=1:4
        ue(2*a-1) = u(2*conn(a)-1);
        ue(2*a)   = u(2*conn(a));
        ueOld(2*a-1) = uOld(2*conn(a)-1);
        ueOld(2*a)   = uOld(2*conn(a));
    end

    [N, dNdxi] = shape_Q4(0,0);
    Rnod = Xe(:,1); Znod = Xe(:,2);
    rnod = Rnod + ue(1:2:end); znod = Znod + ue(2:2:end);
    rnodOld = Rnod + ueOld(1:2:end); znodOld = Znod + ueOld(2:2:end);
    Jmat = Xe.' * dNdxi;
    dNdX = dNdxi / Jmat;

    Rg = N.' * Rnod; rg = N.' * rnod;
    drdR = dNdX(:,1).' * rnod; drdZ = dNdX(:,2).' * rnod;
    dzdR = dNdX(:,1).' * znod; dzdZ = dNdX(:,2).' * znod;
    F = [drdR, 0, drdZ; 0, rg/Rg, 0; dzdR, 0, dzdZ];

    rgOld = N.' * rnodOld;
    drdROld = dNdX(:,1).' * rnodOld; drdZOld = dNdX(:,2).' * rnodOld;
    dzdROld = dNdX(:,1).' * znodOld; dzdZOld = dNdX(:,2).' * znodOld;
    Fold = [drdROld, 0, drdZOld; 0, rgOld/Rg, 0; dzdROld, 0, dzdZOld];

    Jdet = det(F);
    fprintf('%s element %d (r=%.4f, z=%.4f um): J=%.6f\n', label, elemIdx, Rg*1e6, N.'*Znod*1e6, Jdet);
    if Jdet <= 0
        fprintf('  *** NON-POSITIVE J -- element is inverted/degenerate! ***\n');
        return;
    end
    B = F*F.';
    I3 = eye(3);
    devB = B - (trace(B)/3)*I3;
    Telastic = parUse.Ge * Jdet^(-5/3) * devB + parUse.Ke * (Jdet-1) * I3;

    parVisc = parUse; parVisc.dt = dt;
    [~, kvdata] = objective_kelvin_voigt_piola(F, Fold, parVisc);
    Tfull = Telastic + kvdata.Tvisc;
    indepFull = [Tfull(1,1), Tfull(2,2), Tfull(3,3), Tfull(1,3)];

    solverRR = mean(solverStress.sigma_rr(conn));
    solverTT = mean(solverStress.sigma_tt(conn));
    solverZZ = mean(solverStress.sigma_zz(conn));
    solverRZ = mean(solverStress.sigma_rz(conn));
    solverVec = [solverRR, solverTT, solverZZ, solverRZ];

    relErr = abs(indepFull - solverVec) ./ max(abs(solverVec), 1e-6) * 100;

    fprintf('  Independent (elastic+visc): rr=%.4f tt=%.4f zz=%.4f rz=%.4f\n', indepFull);
    fprintf('  Solver-reported (nodal-avg): rr=%.4f tt=%.4f zz=%.4f rz=%.4f\n', solverVec);
    fprintf('  Relative error [%%]: rr=%.2f tt=%.2f zz=%.2f rz=%.2f\n', relErr);
end

[eL, zL] = find_maxz_element(out.meshL, st.uL);
fprintf('Leukocyte max-z element: %d, z=%.4f um\n', eL, zL*1e6);
report_element('Leukocyte', out.meshL, st.uL, st.uLPrev, dtStep, parLmismatch, eL, stressL);

fprintf('\n');
[eE, zE] = find_maxz_element(out.meshE, st.uE);
fprintf('Endothelium max-z element: %d, z=%.4f um\n', eE, zE*1e6);
report_element('Endothelium', out.meshE, st.uE, st.uEPrev, dtStep, par, eE, stressE);

fprintf('\nSMOKE_TEST_STATUS: OK\n');
