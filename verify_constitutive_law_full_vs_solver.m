%% VERIFY_CONSTITUTIVE_LAW_FULL_VS_SOLVER
% Corrected version: the independent reconstruction now includes BOTH the
% elastic term AND the Kelvin-Voigt viscous term (Tvisc, from
% objective_kelvin_voigt_piola.m), matching the FULL constitutive law
% recover_nodal_stress_axisym_viscoelastic.m actually applies -- not just
% the elastic part. Compares against the same nodal-averaged
% solver-recovery values as before.

clc; close all;
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
    drdR = dNdX(:,1).' * rnod;
    drdZ = dNdX(:,2).' * rnod;
    dzdR = dNdX(:,1).' * znod;
    dzdZ = dNdX(:,2).' * znod;
    F = [drdR, 0, drdZ; 0, rg/Rg, 0; dzdR, 0, dzdZ];

    rgOld = N.' * rnodOld;
    drdROld = dNdX(:,1).' * rnodOld;
    drdZOld = dNdX(:,2).' * rnodOld;
    dzdROld = dNdX(:,1).' * znodOld;
    dzdZOld = dNdX(:,2).' * znodOld;
    Fold = [drdROld, 0, drdZOld; 0, rgOld/Rg, 0; dzdROld, 0, dzdZOld];

    Jdet = det(F);
    B = F*F.';
    I3 = eye(3);
    devB = B - (trace(B)/3)*I3;
    Telastic = parUse.Ge * Jdet^(-5/3) * devB + parUse.Ke * (Jdet-1) * I3;

    parVisc = parUse; parVisc.dt = dt;
    [~, kvdata] = objective_kelvin_voigt_piola(F, Fold, parVisc);
    Tfull = Telastic + kvdata.Tvisc;

    indepElasticOnly = [Telastic(1,1), Telastic(2,2), Telastic(3,3), Telastic(1,3)];
    indepFull = [Tfull(1,1), Tfull(2,2), Tfull(3,3), Tfull(1,3)];

    solverRR = mean(solverStress.sigma_rr(conn));
    solverTT = mean(solverStress.sigma_tt(conn));
    solverZZ = mean(solverStress.sigma_zz(conn));
    solverRZ = mean(solverStress.sigma_rz(conn));
    solverVec = [solverRR, solverTT, solverZZ, solverRZ];

    relErrElasticOnly = abs(indepElasticOnly - solverVec) ./ max(abs(solverVec), 1e-6) * 100;
    relErrFull = abs(indepFull - solverVec) ./ max(abs(solverVec), 1e-6) * 100;

    fprintf('%s element %d (r=%.4f, z=%.4f um):\n', label, elemIdx, Rg*1e6, N.'*Znod*1e6);
    fprintf('  Elastic-only indep:  rr=%9.4f tt=%9.4f zz=%9.4f rz=%9.4f | relErr%% = %.2f %.2f %.2f %.2f\n', ...
        indepElasticOnly, relErrElasticOnly);
    fprintf('  Elastic+visc indep:  rr=%9.4f tt=%9.4f zz=%9.4f rz=%9.4f | relErr%% = %.4f %.4f %.4f %.4f\n', ...
        indepFull, relErrFull);
    fprintf('  Solver-reported:     rr=%9.4f tt=%9.4f zz=%9.4f rz=%9.4f\n\n', solverVec);
end

fprintf('================ LEUKOCYTE ================\n');
nelL = out.meshL.nelem;
for e = [round(nelL*0.25), round(nelL*0.5), round(nelL*0.75)]
    report_element('Leukocyte', out.meshL, st.uL, st.uLPrev, dtStep, parLmismatch, e, stressL);
end

fprintf('================ ENDOTHELIUM ================\n');
nelE = out.meshE.nelem;
for e = [round(nelE*0.25), round(nelE*0.5), round(nelE*0.75)]
    report_element('Endothelium', out.meshE, st.uE, st.uEPrev, dtStep, par, e, stressE);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
