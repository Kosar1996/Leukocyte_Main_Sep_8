function [trL, trE] = compute_bodyfitted_wall_traction(mesh, fluid, par)
%COMPUTE_BODYFITTED_WALL_TRACTION
% Returns traction data in the same normal/tangent sign convention used by
% apply_interface_traction in the solid solver:
%   endothelium: normal = +pressure-like load, tangent = -tauE
%   leukocyte:   normal = -pressure-like load, tangent = +tauL
% The pressure/shear values are extracted from the body-fitted MAC solution.

    P = fluid.P;
    stateForShear = struct();
    stateForShear.deltaL = mesh.deltaL_c(:);
    stateForShear.deltaE = mesh.deltaE_c(:);
    if isfield(fluid,'bc') && isfield(fluid.bc,'uzL')
        stateForShear.UwL = fluid.bc.uzL(:);
    else
        stateForShear.UwL = zeros(mesh.Nz,1);
    end
    if isfield(fluid,'bc') && isfield(fluid.bc,'uzE')
        stateForShear.UwE = fluid.bc.uzE(:);
    else
        stateForShear.UwE = zeros(mesh.Nz,1);
    end
    [tauL, tauE] = estimate_wall_shear_bodyfitted(mesh, fluid, stateForShear, par);
    zc = mesh.zc(:);

    trE = struct();
    trE.z = zc;
    trE.normal = -P(end,:).';  % SIGNTEST: was missing minus sign (see locate_and_interp_fluid_stress.m: srr = -pc + visc)
    trE.tangent = -tauE(:);

    trL = struct();
    trL.z = zc;
    trL.normal = -P(1,:).';
    trL.tangent = tauL(:);
    trL = apply_leukocyte_traction_support(trL, par);
end
