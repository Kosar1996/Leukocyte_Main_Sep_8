%% TEST_AITKEN_TRACTION_CORRECTION
% Every fixed relax value from 0.02 to 0.2 made the traction correction
% monotonically WORSE than doing nothing (iteration 0 stayed the best
% result in sweep_relax_after_metric_fix.m). That's the classic
% loosely-coupled-partitioned-FSI failure mode for a lubrication problem:
% the fluid re-solve after each correction pass reacts strongly to even a
% tiny gap change (pressure ~ 1/h^3), so the correction is always chasing
% a target that has already moved by the time it's measured again.
%
% Aitken's Delta^2 relaxation (opts.useAitkenRelax, wired through
% par.useAitkenTractionCorrectionRelax) replaces the fixed step fraction
% with one derived each pass from how the correction's own residual
% actually changed -- the standard fix for exactly this failure mode.
% This compares it directly against the fixed relax=0.1 baseline.

clc; close all;
cd(fileparts(mfilename('fullpath')));

maxIter = 20;

fprintf('\n============ Fixed relax = 0.10 (baseline) ============\n');
cfgFixed = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', true);
cfgFixed.fluid.bodyFittedTractionCorrectionRelax = 0.1;
cfgFixed.parOverrides.bodyFittedTractionCorrectionRelax = 0.1;
cfgFixed.fluid.maxBodyFittedTractionCorrections = maxIter;
cfgFixed.parOverrides.maxBodyFittedTractionCorrections = maxIter;
outFixed = softlube_run_case_global_coupled(cfgFixed);
ciFixed = outFixed.tractionCorrectionHistory{1};

fprintf('\n============ Aitken relaxation ============\n');
cfgAitken = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', true);
cfgAitken.fluid.bodyFittedTractionCorrectionRelax = 0.1; % initial guess for pass 1 only
cfgAitken.parOverrides.bodyFittedTractionCorrectionRelax = 0.1;
cfgAitken.fluid.maxBodyFittedTractionCorrections = maxIter;
cfgAitken.parOverrides.maxBodyFittedTractionCorrections = maxIter;
cfgAitken.fluid.useAitkenTractionCorrectionRelax = true;
cfgAitken.parOverrides.useAitkenTractionCorrectionRelax = true;
outAitken = softlube_run_case_global_coupled(cfgAitken);
ciAitken = outAitken.tractionCorrectionHistory{1};

fprintf('\n\n=================== SUMMARY ===================\n');
fprintf('Fixed  relax=0.10: converged=%d  finalWorstPct=%.2f%%  iters=%d\n', ...
    ciFixed.converged, ciFixed.finalWorstPct, ciFixed.iterations);
fprintf('Aitken relaxation: converged=%d  finalWorstPct=%.2f%%  iters=%d\n', ...
    ciAitken.converged, ciAitken.finalWorstPct, ciAitken.iterations);
fprintf('\nSMOKE_TEST_STATUS: OK\n');
