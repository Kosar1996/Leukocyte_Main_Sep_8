%% SWEEP_RELAX_AFTER_METRIC_FIX
% Now that the endothelium query-range bug (compute_interface_traction_
% mismatch.m dropping 65% of its sample points) is fixed, the starting
% point is already close to target (worst=10.84% at iter 0, vs the old
% buggy 78.27%). At relax=0.30 the correction OVERSHOOTS on pass 1
% (10.84% -> 12.44%) and only partially recovers to a ~11.62% plateau,
% never beating the uncorrected starting point. Try smaller relax values
% to see if a gentler step avoids the overshoot and actually gets under
% the 10% target.

clc; close all;
cd(fileparts(mfilename('fullpath')));

relaxValues = [0.02, 0.05, 0.1, 0.15, 0.2];
maxIter = 20;

results = struct('relax', {}, 'converged', {}, 'finalWorstPct', {}, 'iterUsed', {});

for i = 1:numel(relaxValues)
    r = relaxValues(i);
    cfg = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', true);
    cfg.fluid.bodyFittedTractionCorrectionRelax = r;
    cfg.parOverrides.bodyFittedTractionCorrectionRelax = r;
    cfg.fluid.maxBodyFittedTractionCorrections = maxIter;
    cfg.parOverrides.maxBodyFittedTractionCorrections = maxIter;

    fprintf('\n============================================================\n');
    fprintf('relax = %.3f\n', r);
    fprintf('============================================================\n');
    out = softlube_run_case_global_coupled(cfg);

    ci = out.tractionCorrectionHistory{1};
    results(i).relax = r;
    results(i).converged = ci.converged;
    results(i).finalWorstPct = ci.finalWorstPct;
    results(i).iterUsed = ci.iterations;
end

fprintf('\n\n=================== SUMMARY ===================\n');
fprintf('%8s | %10s | %10s | %6s\n', 'relax', 'converged', 'worstPct', 'iters');
for i = 1:numel(results)
    fprintf('%8.3f | %10d | %10.2f | %6d\n', results(i).relax, results(i).converged, ...
        results(i).finalWorstPct, results(i).iterUsed);
end
fprintf('\nSMOKE_TEST_STATUS: OK\n');
