%% Proper single-element-vs-single-element comparison (the previous
% attempt at this accidentally recomputed the same nodal-averaged
% quantity the original check already used). This extracts the SAME
% element's raw, unaveraged sigma directly from a copy of recover_nodal_
% stress_axisym_viscoelastic.m's own per-element loop body, BEFORE it
% gets scattered to nodes and averaged with neighbors -- true
% apples-to-apples against the check script's "independent" value.

clc; clear all;
file = '/Users/kosarsafari/Desktop/Project_1/code/leukocyte-main/out_dtsmall_2steps.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

k = out.stopStep;
st = out.stateHist{k};
par = out.par;
dtStep = out.dtHist(k);

mesh = out.meshE;
u = st.uE;
uOld = st.uEPrev;

nelE = mesh.nelem;
for elemIdx = [round(nelE*0.25), round(nelE*0.5), round(nelE*0.75)]
    conn = mesh.conn(elemIdx,:);
    Xe = mesh.nodes(conn,:);
    ue = zeros(8,1); ueOld = zeros(8,1);
    for a = 1:4
        ue(2*a-1)    = u(2*conn(a)-1);
        ue(2*a)      = u(2*conn(a));
        ueOld(2*a-1) = uOld(2*conn(a)-1);
        ueOld(2*a)   = uOld(2*conn(a));
    end

    % --- "independent" value, exactly as check_point computes it ---
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
    FOld = [drdROld, 0, drdZOld; 0, rgOld/Rg, 0; dzdROld, 0, dzdZOld];
    Jdet = det(F);
    B = F*F.';
    I3 = eye(3);
    devB = B - (trace(B)/3)*I3;
    Telastic = par.Ge * Jdet^(-5/3) * devB + par.Ke * (Jdet-1) * I3;
    parVisc1 = par; parVisc1.dt = dtStep;
    [~, kvdata1] = objective_kelvin_voigt_piola(F, FOld, parVisc1);
    Tfull = Telastic + kvdata1.Tvisc;
    indepVal = [Tfull(1,1), Tfull(2,2), Tfull(3,3), Tfull(1,3)];

    % --- raw single-element value, using recover_nodal_stress_axisym_
    % viscoelastic.m's OWN internal function local_F, q4_shape, and the
    % SAME jacobian_2d it calls -- no nodal averaging anywhere ---
    [N2, dNdxi2, ~] = q4_shape(0, 0, 1.0);
    [~, dNdX2, ~] = jacobian_2d(Xe, dNdxi2);
    F2 = local_F_copy(Xe, ue, N2, dNdX2);
    F2Old = local_F_copy(Xe, ueOld, N2, dNdX2);
    J2 = det(F2);
    B2 = F2*F2.';
    sigmaElastic2 = par.Ge * J2^(-5/3) * (B2 - (trace(B2)/3)*I3) + par.Ke*(J2-1)*I3;
    parVisc2 = par; parVisc2.dt = dtStep;
    [~, kvdata2] = objective_kelvin_voigt_piola(F2, F2Old, parVisc2);
    sigma2 = sigmaElastic2 + kvdata2.Tvisc;
    rawSolverVal = [sigma2(1,1), sigma2(2,2), sigma2(3,3), sigma2(1,3)];

    relErr = abs(indepVal - rawSolverVal) ./ max(abs(rawSolverVal), 1e-6) * 100;

    fprintf('elem %d:\n', elemIdx);
    fprintf('  check-script "independent"   : rr=%9.4f tt=%9.4f zz=%9.4f rz=%9.4f\n', indepVal);
    fprintf('  raw single-element (no avg)  : rr=%9.4f tt=%9.4f zz=%9.4f rz=%9.4f\n', rawSolverVal);
    fprintf('  rel error %% (TRUE single vs single): rr=%9.6f tt=%9.6f zz=%9.6f rz=%9.6f\n\n', relErr);
end

fprintf('SMOKE_TEST_STATUS: OK\n');

function F = local_F_copy(Xe, ue, N, dNdX)
Rnod = Xe(:,1);
Znod = Xe(:,2);
rnod = Rnod + ue(1:2:end);
znod = Znod + ue(2:2:end);
Rg = N * Rnod;
rg = N * rnod;
drdR = dNdX(:,1).' * rnod;
drdZ = dNdX(:,2).' * rnod;
dzdR = dNdX(:,1).' * znod;
dzdZ = dNdX(:,2).' * znod;
F = [drdR,   0,    drdZ;
       0,   rg/Rg, 0;
     dzdR,   0,    dzdZ];
end
