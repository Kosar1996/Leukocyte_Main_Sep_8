%% VERIFY_TOLERANCE_BUG_FIX
% Permanent, re-runnable verification of the Aug 24 false-convergence fix.
%
% Uses REAL captured production correction-call data (endothelium and
% leukocyte, from an actual coupled run). Runs the solid solve twice on
% the identical input: once with the old, mis-scaled absolute tolerances
% (solidAbsTol=1e-10, solidFallbackAbsTol=1e-8, both tuned for the main
% monolithic solve's force scale, not this micro-scale correction loop)
% reintroduced, and once with today's fix (those tolerances stripped).
%
% Expected result: the OLD-tolerance run falsely reports "converged"
% after very few iterations while the true relative residual is still
% large; the FIXED run genuinely converges to near machine precision.
clc;

fprintf('=== ENDOTHELIUM ===\n');
S = load('/tmp/corr_call_dump_-6000_1.mat');
run_comparison(S.meshE, S.uOldE, S.tractionE, S.interfaceE, S.baseE, S.supportE, S.parCorrE, S.uInitE);

fprintf('\n=== LEUKOCYTE ===\n');
SL = load('/tmp/corr_L_call_dump.mat');
run_comparison(SL.meshL, SL.uOldL, SL.tractionL, SL.interfaceL, SL.baseL, SL.supportL, SL.parCorrL, SL.uInitL);

fprintf('\nSMOKE_TEST_STATUS: OK\n');

function run_comparison(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, parFixed, uInit)
    % Reconstruct the OLD (buggy) parameter set by adding back the
    % mis-scaled tolerances that were stripped by the Aug 24 fix.
    parOld = parFixed;
    parOld.solidAbsTol = 1e-10;
    parOld.solidFallbackAbsTol = 1e-8;

    fprintf('--- OLD (buggy) tolerances: solidAbsTol=1e-10, solidFallbackAbsTol=1e-8 ---\n');
    uOldTol = solve_finite_def_solid(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, parOld, uInit);
    report_quality(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, uOldTol, parOld);

    fprintf('\n--- FIXED (current production) tolerances: stripped, internal refNormFloor default 1e-6 ---\n');
    uFixed = solve_finite_def_solid(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, parFixed, uInit);
    report_quality(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, uFixed, parFixed);
end

function report_quality(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, u, par)
    % Compute the TRUE relative residual for this converged state, using
    % the identical formula the solver itself uses internally (Fint+Fvisc
    % minus Fext, normalized by the larger of the two force scales, over
    % the FREE degrees of freedom only -- fixed/boundary DOFs carry
    % reaction forces that are never supposed to go to zero and are
    % correctly excluded from the solver's own convergence check, so
    % they must be excluded here too for this to mean the same thing).
    ndof = size(mesh.nodes,1)*2;
    [fixDofs, ~] = solid_support_conditions(baseNodes, supportType);
    free = setdiff((1:ndof).', unique(fixDofs(:)));

    Fext = zeros(ndof,1);
    [Fext, ~] = apply_interface_traction(mesh, u, Fext, interfaceNodes, traction);
    [Fint, ~] = assemble_finite_def_axisym(mesh, u, par);
    [Fvisc, ~] = assemble_axisym_kelvin_voigt_viscous(mesh, u, uOld, par);
    Fint = Fint + Fvisc;
    refNorm = max([norm(Fext(free), inf), norm(Fint(free), inf), 1e-14]);
    relNorm = norm(Fint(free) - Fext(free), inf) / refNorm;
    fprintf('    true relative residual (free DOFs only) = %.4e (%.2f%%)\n', relNorm, relNorm*100);
end
