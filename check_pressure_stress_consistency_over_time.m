function summary = check_pressure_stress_consistency_over_time(out)
%CHECK_PRESSURE_STRESS_CONSISTENCY_OVER_TIME
% Repeats the pressure vs. 1/3-trace(sigma) consistency check (see
% check_pressure_stress_consistency.m) at every accepted timestep in
% out, and reports how the mismatch compares to the pressure scale and
% to the viscous stress scale as the simulation evolves in time.
%
% Uses the 90th percentile of the mismatch (not the domain-wide max) as
% "the" mismatch at each step, since the earlier single-step
% investigation found the domain-wide max is dominated by an isolated
% artifact at the outer-wall/grid-stitching seam, while the 90th
% percentile is the value that actually matched the 0.37 Pa reported
% for step 1.
%
% NOTE: if out.stopStep == 1 (a single-timestep run, e.g. from
% run_input_full2D_pressure2.m / build_cfg_full2D_pressure2.m with the
% default nSteps=1), this will only have one point to show. Rerun with
% more steps, e.g.:
%   cfg = build_cfg_full2D_pressure2('nSteps', 10);
%   out = softlube_run_case_global_coupled(cfg);
% to see genuine time evolution.
%
%   summary = check_pressure_stress_consistency_over_time(out)
%
% Fields of summary (each a vector of length out.stopStep):
%   t              simulation time [s]
%   diffMax        max_e |pRaw - (-(1/3)trace(sigma))| [Pa], domain-wide
%   diffP90        90th percentile of the same mismatch [Pa] ("typical")
%   maxP           max |pRaw| [Pa]
%   diffRelToP     diffP90 / maxP
%   viscousScale   representative interior viscous normal-stress scale [Pa]
%   diffRelToVisc  diffP90 / viscousScale

mu = out.par.mu;
nSteps = out.stopStep;

t = nan(nSteps,1);
diffMax = nan(nSteps,1);
diffP90 = nan(nSteps,1);
maxP = nan(nSteps,1);
viscousScale = nan(nSteps,1);

for k = 1:nSteps
    fluid = out.fluidHist{k};
    if isempty(fluid) || ~isfield(fluid,'meshF') || isempty(fluid.meshF)
        continue;
    end

    meshF = add_fluid_nodes(fluid.meshF);
    [~, sigmaCell, ~, pRaw] = recover_fluid_nodes_pressure_stress_Q4( ...
        meshF, fluid.ur2D, fluid.uz2D, mu, fluid.pCell);

    traceThird = -(sigmaCell(:,1) + sigmaCell(:,2) + sigmaCell(:,3)) / 3;
    diffAll = abs(pRaw(:) - traceThird(:));
    valid = isfinite(diffAll);

    t(k) = out.t(k);
    diffMax(k) = max(diffAll(valid));
    diffP90(k) = prctile(diffAll(valid), 90);
    maxP(k) = max(abs(pRaw(valid)));

    viscNormal = abs(sigmaCell(:,1:3) + pRaw(:));
    viscNormalFlat = viscNormal(:);
    viscNormalFinite = viscNormalFlat(isfinite(viscNormalFlat) & viscNormalFlat > 0);
    if ~isempty(viscNormalFinite)
        viscousScale(k) = median(viscNormalFinite);
    end
end

diffRelToP = diffP90 ./ max(maxP, eps);
diffRelToVisc = diffP90 ./ max(viscousScale, eps);

summary = struct('t', t, 'diffMax', diffMax, 'diffP90', diffP90, ...
    'maxP', maxP, 'viscousScale', viscousScale, ...
    'diffRelToP', diffRelToP, 'diffRelToVisc', diffRelToVisc);

fprintf('\n=== Pressure/stress-trace consistency over time ===\n');
fprintf('%10s %14s %14s %14s %12s %12s\n', ...
    't [s]', 'diffMax[Pa]', 'diffP90[Pa]', 'maxP[Pa]', '% of maxP', '% of visc');
for k = 1:nSteps
    if isnan(t(k)), continue; end
    fprintf('%10.4e %14.6e %14.6e %14.6e %12.4f %12.4f\n', ...
        t(k), diffMax(k), diffP90(k), maxP(k), ...
        100*diffRelToP(k), 100*diffRelToVisc(k));
end

if nnz(~isnan(t)) >= 2
    figure;
    yyaxis left;
    plot(t*1e3, diffP90, 'o-', 'LineWidth', 2); hold on;
    plot(t*1e3, maxP, 's--', 'LineWidth', 1.5);
    ylabel('Pa');
    yyaxis right;
    plot(t*1e3, 100*diffRelToP, '^:', 'LineWidth', 1.5);
    ylabel('mismatch as % of max pressure');
    xlabel('t [ms]');
    legend({'90th-pctile mismatch', 'max pressure', '% of max pressure'}, ...
        'Location', 'best');
    title('Pressure/stress-trace consistency vs. simulation time');
    grid on;
else
    fprintf(['\nOnly one accepted timestep is available (out.stopStep=%d), ', ...
        'so there is nothing to plot over time yet. Rerun with nSteps > 1 ', ...
        '(e.g. build_cfg_full2D_pressure2(''nSteps'', 10)) to see genuine ', ...
        'time evolution.\n'], nSteps);
end

end
