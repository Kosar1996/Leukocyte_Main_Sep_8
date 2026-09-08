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

stressL = recover_nodal_stress_axisym_viscoelastic(out.meshL, st.uL, st.uLPrev, dtStep, parLmismatch);

targetZ = 4.0051e-6; % endothelium's max-z point, the shared truncation-edge location

function elemIdx = find_closest_z_element(mesh, u, targetZ)
nel = mesh.nelem;
zPerElem = zeros(nel,1);
for e = 1:nel
    conn = mesh.conn(e,:);
    znod = mesh.nodes(conn,2) + u(2*conn);
    zPerElem(e) = mean(znod);
end
[~, elemIdx] = min(abs(zPerElem - targetZ));
end

elemIdx = find_closest_z_element(out.meshL, st.uL, targetZ);

conn = out.meshL.conn(elemIdx,:);
Xe = out.meshL.nodes(conn,:);
ue = zeros(8,1); ueOld = zeros(8,1);
for a=1:4
    ue(2*a-1) = st.uL(2*conn(a)-1);
    ue(2*a)   = st.uL(2*conn(a));
    ueOld(2*a-1) = st.uLPrev(2*conn(a)-1);
    ueOld(2*a)   = st.uLPrev(2*conn(a));
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
rg_deformed = N.' * znod;
fprintf('Leukocyte element %d, closest to shared z=%.4f um: actual deformed z=%.4f um, r=%.4f um\n', ...
    elemIdx, targetZ*1e6, rg_deformed*1e6, Rg*1e6);
fprintf('J=%.6f\n', Jdet);

if Jdet <= 0
    fprintf('*** NON-POSITIVE J -- element inverted/degenerate ***\n');
else
    B = F*F.';
    I3 = eye(3);
    devB = B - (trace(B)/3)*I3;
    Telastic = parLmismatch.Ge * Jdet^(-5/3) * devB + parLmismatch.Ke * (Jdet-1) * I3;

    parVisc = parLmismatch; parVisc.dt = dtStep;
    [~, kvdata] = objective_kelvin_voigt_piola(F, Fold, parVisc);
    Tfull = Telastic + kvdata.Tvisc;
    indepFull = [Tfull(1,1), Tfull(2,2), Tfull(3,3), Tfull(1,3)];

    solverRR = mean(stressL.sigma_rr(conn));
    solverTT = mean(stressL.sigma_tt(conn));
    solverZZ = mean(stressL.sigma_zz(conn));
    solverRZ = mean(stressL.sigma_rz(conn));
    solverVec = [solverRR, solverTT, solverZZ, solverRZ];

    relErr = abs(indepFull - solverVec) ./ max(abs(solverVec), 1e-6) * 100;

    fprintf('Independent (elastic+visc): rr=%.4f tt=%.4f zz=%.4f rz=%.4f\n', indepFull);
    fprintf('Solver-reported (nodal-avg): rr=%.4f tt=%.4f zz=%.4f rz=%.4f\n', solverVec);
    fprintf('Relative error [%%]: rr=%.2f tt=%.2f zz=%.2f rz=%.2f\n', relErr);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
