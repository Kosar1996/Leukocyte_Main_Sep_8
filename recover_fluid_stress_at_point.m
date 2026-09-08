function [sigma, p, rz_actual, xieta] = recover_fluid_stress_at_point(mraw, meshF, ur2D, uz2D, mu, plist, i, j, rq, zq)
%RECOVER_FLUID_STRESS_AT_POINT

% always evaluates shape_Q4(0,0) -- the element centroid -- regardless of
% where the actual query point (rq,zq) sits inside the element. This
% function instead solves for the TRUE local isoparametric coordinates
% (xi,eta) of the query point within its containing element (Newton
% iteration on the bilinear Q4 map), then evaluates the shape functions
% and their derivatives there, giving a genuinely position-dependent
% interpolated stress instead of a constant-per-element value.
%
% i,j: structured-grid element indices (as used elsewhere in this
%      codebase, e.g. full_redo_1234_complete.m) identifying which
%      element contains (rq,zq).

Nr = mraw.Nr;
e = (j-1)*(Nr-1) + i;
conn = meshF.elems(e,:);
xe = meshF.nodes(conn,1);
ze = meshF.nodes(conn,2);

% Newton iteration to invert the bilinear map: find (xi,eta) such that
% N(xi,eta)'*xe = rq, N(xi,eta)'*ze = zq
xi = 0; eta = 0;
for it = 1:50
    [N, dNdxi] = shape_Q4(xi, eta);
    r_cur = N.' * xe;
    z_cur = N.' * ze;
    resid = [rq - r_cur; zq - z_cur];
    if norm(resid) < 1e-14
        break;
    end
    J = [xe ze].' * dNdxi;   % 2x2 Jacobian d(r,z)/d(xi,eta)
    dxieta = J \ resid;
    xi = xi + dxieta(1);
    eta = eta + dxieta(2);
    xi = min(max(xi,-1.5),1.5);
    eta = min(max(eta,-1.5),1.5);
end
xieta = [xi, eta];

[N, dNdxi] = shape_Q4(xi, eta);
J = [xe ze].' * dNdxi;
dNdx = dNdxi / J;
dNdr = dNdx(:,1);
dNdz = dNdx(:,2);

r_actual = N.' * xe;
z_actual = N.' * ze;
rz_actual = [r_actual, z_actual];

pc = N.' * plist(conn);

ur = zeros(4,1); uz = zeros(4,1);
for a = 1:4
    ur(a) = ur2D(conn(a));
    uz(a) = uz2D(conn(a));
end

durdr = dNdr.' * ur;
duzdz = dNdz.' * uz;
durdz = dNdz.' * ur;
duzdr = dNdr.' * uz;
ur_over_r = (N.' * ur) / r_actual;

srr = -pc + 2*mu*durdr;
stt = -pc + 2*mu*ur_over_r;
szz = -pc + 2*mu*duzdz;
srz = mu*(durdz + duzdr);

sigma = [srr, stt, szz, srz];
p = pc;
end
