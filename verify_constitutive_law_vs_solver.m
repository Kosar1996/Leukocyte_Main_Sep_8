%% VERIFY_CONSTITUTIVE_LAW_VS_SOLVER
% Extends verify_constitutive_law_points.m: for the same 6 sample
% elements, compares the independently-reconstructed stress (computed
% from scratch via F/J/B/devB/T, element-centroid, no averaging) against
% the solver's own stress-recovery pipeline
% (recover_nodal_stress_axisym_viscoelastic.m, nodal-averaged across all
% elements touching each node) evaluated at the same element's 4 corner
% nodes, averaged. This is a genuine cross-check: two independently
% implemented paths through the same constitutive law, compared at the
% same physical location.

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

stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, out.dtHist(k), par);
stressL = recover_nodal_stress_axisym_viscoelastic(out.meshL, st.uL, st.uLPrev, out.dtHist(k), parLmismatch);

function report_element(label, mesh, u, parUse, elemIdx, solverStress)
    conn = mesh.conn(elemIdx,:);
    Xe = mesh.nodes(conn,:);
    ue = zeros(8,1);
    for a=1:4
        ue(2*a-1) = u(2*conn(a)-1);
        ue(2*a)   = u(2*conn(a));
    end

    [N, dNdxi] = shape_Q4(0,0);
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
    indep = [T(1,1), T(2,2), T(3,3), T(1,3)]; % rr, tt, zz, rz

    % Solver-reported: average the 4 corner nodes' solver-recovered stress
    solverRR = mean(solverStress.sigma_rr(conn));
    solverTT = mean(solverStress.sigma_tt(conn));
    solverZZ = mean(solverStress.sigma_zz(conn));
    solverRZ = mean(solverStress.sigma_rz(conn));
    solverVec = [solverRR, solverTT, solverZZ, solverRZ];

    relErr = abs(indep - solverVec) ./ max(abs(solverVec), 1e-6) * 100;

    fprintf('%s element %d (r=%.4f, z=%.4f um):\n', label, elemIdx, Rg*1e6, N.'*Znod*1e6);
    fprintf('  Independent (from-scratch F/J/B/devB/T): rr=%.4f tt=%.4f zz=%.4f rz=%.4f\n', indep);
    fprintf('  Solver-reported (nodal-averaged recovery): rr=%.4f tt=%.4f zz=%.4f rz=%.4f\n', solverVec);
    fprintf('  Relative error [%%]: rr=%.4f tt=%.4f zz=%.4f rz=%.4f\n\n', relErr);
end

fprintf('================ LEUKOCYTE ================\n');
nelL = out.meshL.nelem;
for e = [round(nelL*0.25), round(nelL*0.5), round(nelL*0.75)]
    report_element('Leukocyte', out.meshL, st.uL, parLmismatch, e, stressL);
end

fprintf('================ ENDOTHELIUM ================\n');
nelE = out.meshE.nelem;
for e = [round(nelE*0.25), round(nelE*0.5), round(nelE*0.75)]
    report_element('Endothelium', out.meshE, st.uE, par, e, stressE);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
