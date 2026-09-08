%% TEST_FINESTEP_UNSTEADY_PPLOW_FROM_RUNB_STATE1
% Sanity-check companion to test_finestep_UNSTEADY_PPFIX_from_runB_state1.m
% (Issue 3, pressurePenalty=1e5). That run avoided the persistent
% correction-loop oscillation and converged (exitflag=3) cleanly for 5
% consecutive steps -- but the resulting pressures were suspiciously tiny
% (0.04-0.29 Pa, vs 25-14,000 Pa in every other test this session). That
% raises a real concern: pressurePenalty=1e5 may be so large it
% artificially decouples pressure from the actual incompressibility
% constraint (making the fluid effectively "too compressible"),
% suppressing the physical pressure field rather than genuinely
% stabilizing the solve -- a correction loop with near-zero corrective
% traction "converges" trivially regardless of whether that traction is
% physically real.
%
% This test uses a much smaller penalty (100, matching the more modest of
% the two quick single-step probes run earlier today) to see whether
% convergence still holds without suppressing the pressure field this
% much. If PP=100 both avoids the oscillation AND produces pressures in
% the physically-expected range (tens to low thousands of Pa, consistent
% with the known baseline), that is much stronger evidence of a genuine
% fix than PP=1e5 alone. If PP=100 falls back into the same oscillation
% baseline, that would suggest the effect is a threshold/magnitude-
% sensitive artifact rather than a robust fix.
%
% Uses the ORIGINAL (buggy rOuter=4e-6) starting state, unchanged, same
% as every other isolated test today. Relax kept at its default (0.5),
% matching the PPFIX run so only the penalty magnitude differs between
% the two. Writes to its own output/log files only.

clc;

S = load('out_pure2dmac_dtlarge_7steps.mat');
out = S.out;

fprintf('Run B step 1: t=%.6e s, dt used=%.6e s\n', out.t(1), out.dtHist(1));
fprintf('Run B step 2 (actual, coarse dt, steady): t=%.6e s, max|p|=%.4f Pa\n', ...
    out.t(2), max(abs(out.fluidHist{2}.p)));

state1 = out.stateHist{1};
par = out.par;
par.dt = 3e-4;
par.useUnsteadyStokes = true;
par.rho = 1000;
nFineSteps = 4;
par.tEnd = state1.t + nFineSteps * par.dt;

fprintf('Relaxation factor (unchanged): par.bodyFittedTractionCorrectionRelax = %.4f\n', ...
    par.bodyFittedTractionCorrectionRelax);
par.fluidPressurePenalty = 100;
fprintf('PROBE: par.fluidPressurePenalty = %g\n', par.fluidPressurePenalty);

restart = struct();
restart.state = state1;
restart.par = par;
restart.meshE = out.meshE;
restart.interfaceE = out.interfaceE;
restart.meshL = out.meshL;
restart.interfaceL = out.interfaceL;
restart.z = out.z;
restart.PHist = out.PHist;
restart.urCHist = out.urCHist;
restart.uzCHist = out.uzCHist;
restart.speedCHist = out.speedCHist;
restart.RPHist = out.RPHist;
restart.ZPHist = out.ZPHist;

cfg = struct();
cfg.geometry.endotheliumPrestressFile = 'solid_endo_P300_wide.mat';
cfg.geometry.leukocytePrestressFile = 'solid_leu_P600.mat';

outFine = softlube_run_case_global_coupled(cfg, restart);

% Save FIRST, before any post-processing that could crash and lose the run.
save('out_finestep_UNSTEADY_PPLOW.mat', 'outFine', '-v7.3');
fprintf('\nSaved full output to out_finestep_UNSTEADY_PPLOW.mat\n');

fprintf('\n=== RESULT (unsteady Stokes ON, both warm-start fixes + fluidPressurePenalty=100 ACTIVE) ===\n');
for k = 1:outFine.stopStep
    if k <= numel(outFine.fluidHist) && isstruct(outFine.fluidHist{k}) && isfield(outFine.fluidHist{k}, 'p')
        fprintf('Fine step %d: t=%.6e s, max|p|=%.4f Pa\n', ...
            k, outFine.t(k), max(abs(outFine.fluidHist{k}.p)));
    else
        fprintf('Fine step %d: t=%.6e s, INVALID fluidHist entry (class=%s)\n', ...
            k, outFine.t(k), class(outFine.fluidHist{k}));
    end
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
