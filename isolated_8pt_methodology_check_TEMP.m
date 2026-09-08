%% Verify whether the 8-point check's "FAIL" results are a genuine formula
% bug or a nodal-averaging vs single-element methodology mismatch.
% recover_nodal_stress_axisym_viscoelastic.m reports NODAL-AVERAGED stress
% (each node blends every adjacent element's single-Gauss-point value).
% verify_solid_mechanics_8points_after_dtfix.m's "independent" value is
% ONE element's own unaveraged center value. Comparing those directly is
% not apples-to-apples wherever stress varies between neighbors. This
% recomputes the SAME single-element quantity recover_nodal_stress_axisym_
% viscoelastic.m computes internally (before scattering to nodes) and
% compares that directly against the check script's "independent" value.

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

    parVisc = par; parVisc.dt = dtStep;
    [~, kvdata] = objective_kelvin_voigt_piola(F, FOld, parVisc);
    Tfull = Telastic + kvdata.Tvisc;
    indepSingleElem = [Tfull(1,1), Tfull(2,2), Tfull(3,3), Tfull(1,3)];

    % Nodal-averaged value from the actual solver function, exactly as
    % the 8-point check does it: mean over the element's 4 corner nodes.
    stressFull = recover_nodal_stress_axisym_viscoelastic(mesh, u, uOld, dtStep, par);
    solverNodalAvg = [mean(stressFull.sigma_rr(conn)), mean(stressFull.sigma_tt(conn)), ...
                       mean(stressFull.sigma_zz(conn)), mean(stressFull.sigma_rz(conn))];

    relErrVsNodalAvg = abs(indepSingleElem - solverNodalAvg) ./ max(abs(solverNodalAvg), 1e-6) * 100;

    fprintf('elem %d:\n', elemIdx);
    fprintf('  single-element (independent) : rr=%9.4f tt=%9.4f zz=%9.4f rz=%9.4f\n', indepSingleElem);
    fprintf('  solver nodal-averaged         : rr=%9.4f tt=%9.4f zz=%9.4f rz=%9.4f\n', solverNodalAvg);
    fprintf('  rel error %% (single vs nodal-avg): rr=%7.3f tt=%7.3f zz=%7.3f rz=%7.3f\n\n', relErrVsNodalAvg);
end

fprintf('SMOKE_TEST_STATUS: OK\n');
