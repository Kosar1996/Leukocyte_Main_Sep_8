function uNew = solve_finite_def_solid(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial)
%SOLVE_FINITE_DEF_SOLID
% Finite-deformation (neo-Hookean + Kelvin-Voigt viscous) Newton solve for
% a single solid under a prescribed interface traction (normal + tangent).
%
% Extracted from apply_bodyfitted_MAC_traction_correction.m (was a private
% local function there) so it can be shared with
% apply_bodyfitted_MAC_traction_correction_feedback.m without duplicating
% the Newton/line-search/trust-region logic. No behavior change from the
% original -- same code, just made reusable across files.

    ndof = size(mesh.nodes,1)*2;
    if nargin >= 8 && ~isempty(uInitial)
        u = uInitial;
    else
        u = uOld;
    end

    [fixDofs, fixVals] = solid_support_conditions(baseNodes, supportType);
    free = setdiff((1:ndof).', unique(fixDofs(:)));
    u(fixDofs) = fixVals;

    absTol = 1e-13;
    if isfield(par, 'solidAbsTol')
        absTol = par.solidAbsTol;
    end

    fallbackAbsTol = [];
    if isfield(par, 'solidFallbackAbsTol')
        fallbackAbsTol = par.solidFallbackAbsTol;
    end
    useFallbackAbsTol = ~isempty(fallbackAbsTol) && ...
        isfinite(fallbackAbsTol) && fallbackAbsTol > 0;

    fallbackRelTol = inf;
    if isfield(par, 'solidFallbackRelTol')
        fallbackRelTol = par.solidFallbackRelTol;
    end
    fallbackMinIterations = 2;
    if isfield(par, 'solidFallbackMinIterations')
        fallbackMinIterations = max(1, round(par.solidFallbackMinIterations));
    end

    trustU = par.trustU0;
    if isfield(par, 'solidTrustU0')
        trustU = par.solidTrustU0;
    end

    trustUMin = par.trustUMin;
    if isfield(par, 'solidTrustUMin')
        trustUMin = par.solidTrustUMin;
    end

    trustUMax = par.trustUMax;
    if isfield(par, 'solidTrustUMax')
        trustUMax = par.solidTrustUMax;
    end

    bestRel = inf;
    bestAbs = inf;
    bestU = u;
    bestUpdated = false;

    for it = 1:par.newtonMaxItSolid

        Fext = zeros(ndof,1);
        [Fext, Kext] = apply_interface_traction(mesh, u, Fext, interfaceNodes, traction);

        [Fint, Ktan] = assemble_finite_def_axisym(mesh, u, par);
        [Fvisc, Kvisc] = assemble_axisym_kelvin_voigt_viscous(mesh, u, uOld, par);
        Fint = Fint + Fvisc;
        Ktan = Ktan + Kvisc;

        R = Fint - Fext;
        Ktot = Ktan - Kext;

        Rf  = R(free);
        Kff = Ktot(free, free);

        resNorm = norm(Rf, inf);
        refNorm = max([norm(Fext(free), inf), norm(Fint(free), inf), 1e-14]);
        relNorm = resNorm / refNorm;

        if isfield(par, 'debugVerbose') && par.debugVerbose
            fprintf('      solid Newton %d: rel=%.3e abs=%.3e trustU=%.3e\n', ...
                it, relNorm, resNorm, trustU);
        end

        if relNorm < bestRel
            bestRel = relNorm;
            bestAbs = resNorm;
            bestU = u;
        end

        if relNorm < par.newtonTolSolid || resNorm < absTol
            uNew = u;
            return;
        end

        if useFallbackAbsTol && it >= fallbackMinIterations && ...
                resNorm < fallbackAbsTol && relNorm < fallbackRelTol
            uNew = u;
            return;
        end

        % Diagonal (Jacobi/Ruiz) equilibration before the linear solve.
        % Root-caused Aug 21-24: this system spans ~6 orders of magnitude
        % between its length scale (m) and force scale (N) at this
        % micro-geometry, which leaves Kff badly scaled -- force, both
        % analytical tangent matrices (Kext, Ke), and the constitutive law
        % were all independently verified correct via finite-difference
        % checks, yet a rigorously-benchmarked (exact Lame thick-cylinder
        % solution) test still showed the solve landing 3-6x off from the
        % true stress state, with global force balance violated by 100%+
        % despite a tiny per-DOF residual -- the signature of an
        % ill-conditioned linear system, not a formula bug. This is a
        % mathematically exact change of variables (Kff*du = -Rf becomes
        % (D*Kff*D)*(D^-1*du) = -(D*Rf), solved then unscaled back) that
        % only improves the direct solver's achievable precision; it does
        % not alter the physics or the converged solution in exact
        % arithmetic.
        dscale = sqrt(abs(full(diag(Kff))));
        dscale(dscale < eps(class(dscale)) | ~isfinite(dscale)) = 1;
        Dinv = spdiags(1./dscale, 0, numel(dscale), numel(dscale));
        Kscaled = Dinv * Kff * Dinv;
        Rscaled = Dinv * Rf;
        du_scaled = -Kscaled \ Rscaled;
        du_free = dscale .\ du_scaled;
        if ~all(isfinite(du_free))
            error('Solid Newton produced non-finite displacement increments.');
        end

        duMax = max(abs(du_free));
        if duMax > trustU
            stepScale = trustU / duMax;
        else
            stepScale = 1.0;
        end

        alpha = 1.0;
        accepted = false;
        resAccepted = inf;

        for ls = 1:par.lineSearchMax
            uTrial = u;
            uTrial(free) = uTrial(free) + alpha * stepScale * du_free;
            uTrial(fixDofs) = fixVals;

            try
                FextTrial = zeros(ndof,1);
                [FextTrial, ~] = apply_interface_traction(mesh, uTrial, FextTrial, interfaceNodes, traction);

                FintTrial = assemble_finite_def_internal_force_only(mesh, uTrial, par) + ...
                    assemble_axisym_kelvin_voigt_viscous_force_only(mesh, uTrial, uOld, par);
                Rtrial = FintTrial - FextTrial;
                resTrial = norm(Rtrial(free), inf);
                refTrial = max([norm(FextTrial(free), inf), norm(FintTrial(free), inf), 1e-14]);
                relTrial = resTrial / refTrial;

                if relTrial < bestRel
                    bestRel = relTrial;
                    bestAbs = resTrial;
                    bestU = uTrial;
                    bestUpdated = true;
                end

                if resTrial < resNorm || resTrial < absTol
                    u = uTrial;
                    resAccepted = resTrial;
                    accepted = true;
                    break;
                end

            catch ME
                if contains(ME.message, 'Negative or zero J') || ...
                   contains(ME.message, 'Non-positive radius') || ...
                   contains(ME.message, 'Element inverted')
                    % Bad trial step: reject and reduce alpha.
                else
                    rethrow(ME);
                end
            end

            alpha = 0.5 * alpha;
        end

        if ~accepted
            trustU = 0.5 * trustU;
            if trustU < trustUMin
                error(['Solid Newton line search failed before equilibrium. ', ...
                       'best relative residual = %.3e, best absolute residual = %.3e.'], ...
                       bestRel, bestAbs);
            end
            continue;
        end

        if resAccepted < 0.25 * resNorm
            trustU = min(2.0 * trustU, trustUMax);
        elseif alpha < 0.25 || stepScale < 0.25
            trustU = max(0.5 * trustU, trustUMin);
        end
    end

    if useFallbackAbsTol && bestUpdated && ...
            bestAbs < fallbackAbsTol && bestRel < fallbackRelTol
        uNew = bestU;
        return;
    end

    error(['Finite-deformation solid solve hit max iterations before equilibrium. ', ...
           'best relative residual = %.3e, best absolute residual = %.3e.'], ...
           bestRel, bestAbs);
end
