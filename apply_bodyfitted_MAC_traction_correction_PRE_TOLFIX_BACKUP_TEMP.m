function [state, fluid, ok, stopReason] = apply_bodyfitted_MAC_traction_correction( ...
    z, old, state, fluid, meshE, interfaceE, baseE, meshL, interfaceL, baseL, parL, par)
%APPLY_BODYFITTED_MAC_TRACTION_CORRECTION
% Performs one or more partitioned solid corrections using the current
% body-fitted MAC fluid traction. This is the nonintrusive way to make the
% post-step solid state feel the 2D Stokes pressure/shear without rebuilding
% the whole monolithic residual around MAC unknowns.

    ok = true;
    stopReason = '';

    nCorr = 1;
    if isfield(par,'maxBodyFittedTractionCorrections') && isfinite(par.maxBodyFittedTractionCorrections)
        nCorr = max(0, round(par.maxBodyFittedTractionCorrections));
    end
    if nCorr == 0
        return;
    end

    relax = 1.0;
    if isfield(par,'bodyFittedTractionCorrectionRelax') && isfinite(par.bodyFittedTractionCorrectionRelax)
        relax = min(1.0, max(0.0, par.bodyFittedTractionCorrectionRelax));
    end

    % Convergence tolerance for the inner partitioned iteration: this loop
    % is a fixed-point iteration between the solid state and the fluid
    % traction, distinct from the outer adaptive time-step retries (which
    % shrink dt and redo the whole step). It should run until the solid
    % correction stops changing meaningfully, not always a fixed count of
    % passes -- some steps may need far fewer than nCorr, others may need
    % close to it. maxBodyFittedTractionCorrections is now the CAP on that
    % iteration, not a fixed count.
    corrTol = 1e-4;
    if isfield(par,'bodyFittedTractionCorrectionTol') && isfinite(par.bodyFittedTractionCorrectionTol)
        corrTol = par.bodyFittedTractionCorrectionTol;
    end

    useRLoutInner = use_RLout_fluid_interface_for_solid_leukocyte(par);
    hasL = ~useRLoutInner && ~isempty(meshL) && ~isempty(interfaceL) && ...
        isfield(state,'uL') && ~isempty(state.uL);

    passesUsed = 0;
    converged = false;

    for ic = 1:nCorr
        if ~isfield(fluid,'tractionE') || ~isfield(fluid,'tractionL')
            [fluid.tractionL, fluid.tractionE] = compute_bodyfitted_wall_traction(fluid.meshF, fluid, par);
        end

        stateBeforeCorr = state;
        fluidBeforeCorr = fluid;
        try
            uEold = state.uE;
            uEcorr = solve_finite_def_solid(meshE, old.uE, fluid.tractionE, ...
                interfaceE, baseE, par.supportE, par, state.uE);
            state.uE = uEold + relax * (uEcorr - uEold);
            [state.deltaE, state.UwE] = monolithic_interface_kinematics_value_only( ...
                meshE, state.uE, old.uE, interfaceE, z, par);

            % Normalize the pass-to-pass change against THIS TIME STEP's
            % own increment (state.uE - old.uE), not the full accumulated
            % displacement history (old.uE can be large after many prior
            % steps, which would make any correction look artificially
            % small and falsely signal convergence after a single pass).
            stepIncrE = norm(uEold - old.uE);
            relChangeE = norm(state.uE - uEold) / max(stepIncrE, 1e-30);

            relChangeL = 0;
            if hasL
                uLold = state.uL;
                uLcorr = solve_finite_def_solid(meshL, old.uL, fluid.tractionL, ...
                    interfaceL, baseL, par.supportL, parL, state.uL);
                state.uL = uLold + relax * (uLcorr - uLold);
                [state.deltaL, state.UwL] = monolithic_interface_kinematics_value_only( ...
                    meshL, state.uL, old.uL, interfaceL, z, par);
                stepIncrL = norm(uLold - old.uL);
                relChangeL = norm(state.uL - uLold) / max(stepIncrL, 1e-30);
            end

            if useRLoutInner
                state.deltaL = par.RLout * ones(size(z));
                state.UwL = zeros(size(z));
            end

            state = attach_physical_solid_interface_fields( ...
                state, old, meshE, interfaceE, meshL, interfaceL, z, par);

            if any(state.deltaE(:) - state.deltaL(:) <= par.minGap)
                error('Body-fitted traction correction produced a gap below minGap.');
            end

            if ~isfield(par,'resolveFluidAfterTractionCorrection') || par.resolveFluidAfterTractionCorrection
                [fluid, ok, stopReason] = solve_selected_poststep_fluid(z, old, state, par);
                if ~ok
                    return;
                end
            end

            passesUsed = ic;
            if max(relChangeE, relChangeL) < corrTol
                converged = true;
                break;
            end
        catch ME
            state = stateBeforeCorr;
            fluid = fluidBeforeCorr;
            stopReason = ME.message;

            failMode = "error";
            if isfield(par, 'bodyFittedTractionCorrectionFailMode') && ...
                    ~isempty(par.bodyFittedTractionCorrectionFailMode)
                failMode = lower(string(par.bodyFittedTractionCorrectionFailMode));
            end

            if failMode == "warn" || failMode == "skip"
                ok = true;
                if failMode == "warn"
                    warning(['Skipping body-fitted MAC traction correction ', ...
                        'after correction %d failed: %s'], ic, ME.message);
                end
                if ~isfield(par,'resolveFluidAfterTractionCorrection') || par.resolveFluidAfterTractionCorrection
                    [fluid, ok, stopReason] = solve_selected_poststep_fluid(z, old, state, par);
                end
                return;
            end

            ok = false;
            return;
        end
    end

    state.tractionCorrectionPassesUsed = passesUsed;
    state.tractionCorrectionConverged = converged;
    if ~converged
        warning(['Body-fitted traction correction did not converge within %d passes ', ...
            '(tol=%.3g). Consider raising maxBodyFittedTractionCorrections.'], nCorr, corrTol);
    end
end

% solve_finite_def_solid was moved to its own file, solve_finite_def_solid.m,
% so apply_bodyfitted_MAC_traction_correction_feedback.m can reuse it without
% duplicating the Newton/line-search/trust-region logic. No behavior change.