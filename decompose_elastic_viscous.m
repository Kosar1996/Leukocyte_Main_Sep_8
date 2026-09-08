function [stressElastic, stressVisc] = decompose_elastic_viscous(mesh, u, uOld, dt, par)
%DECOMPOSE_ELASTIC_VISCOUS
% Same per-element loop as recover_nodal_stress_axisym_viscoelastic.m, but
% returns the elastic and viscous nodal stress SEPARATELY instead of only
% their sum, so the two contributions can be compared directly (used to



nnode = size(mesh.nodes,1);
fields = {'sigma_rr','sigma_tt','sigma_zz','sigma_rz'};
sumsE = struct(); sumsV = struct();
for f = 1:numel(fields)
    sumsE.(fields{f}) = zeros(nnode,1);
    sumsV.(fields{f}) = zeros(nnode,1);
end
count = zeros(nnode,1);

xi = 0.0; eta = 0.0;
[N, dNdxi, ~] = q4_shape(xi, eta, 1.0);

useVisc = isfield(par,'etaE') && isfinite(par.etaE) && par.etaE > 0;
parVisc = par;
parVisc.dt = dt;

for e = 1:mesh.nelem
    conn = mesh.conn(e,:);
    Xe   = mesh.nodes(conn,:);
    dofs = reshape([2*conn-1; 2*conn], [], 1);
    ue    = u(dofs);
    ueOld = uOld(dofs);

    [~, dNdX, ~] = jacobian_2d(Xe, dNdxi);

    F = local_F(Xe, ue, N, dNdX);
    J = det(F);
    B = F * F.';
    I3 = eye(3);
    sigmaElastic = par.Ge * J^(-5/3) * ( B - (trace(B)/3)*I3 ) + par.Ke * (J - 1) * I3;

    if useVisc
        Fold = local_F(Xe, ueOld, N, dNdX);
        [~, kvdata] = objective_kelvin_voigt_piola(F, Fold, parVisc);
        Tvisc = kvdata.Tvisc;
    else
        Tvisc = zeros(3,3);
    end

    valsE = [sigmaElastic(1,1), sigmaElastic(2,2), sigmaElastic(3,3), sigmaElastic(1,3)];
    valsV = [Tvisc(1,1), Tvisc(2,2), Tvisc(3,3), Tvisc(1,3)];

    for a = 1:4
        node = conn(a);
        for f = 1:numel(fields)
            sumsE.(fields{f})(node) = sumsE.(fields{f})(node) + valsE(f);
            sumsV.(fields{f})(node) = sumsV.(fields{f})(node) + valsV(f);
        end
        count(node) = count(node) + 1;
    end
end

stressElastic = struct(); stressVisc = struct();
for f = 1:numel(fields)
    stressElastic.(fields{f}) = sumsE.(fields{f}) ./ max(count,1);
    stressVisc.(fields{f})    = sumsV.(fields{f}) ./ max(count,1);
end
end

function F = local_F(Xe, ue, N, dNdX)
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
