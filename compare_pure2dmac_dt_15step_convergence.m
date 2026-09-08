%% COMPARE_PURE2DMAC_DT_15STEP_CONVERGENCE
% Real dt-sensitivity check for the standalone Pure 2D MAC fluid solve,

% a larger dt over multiple steps (not just step 1), compare ALL
% fluid-related quantities, and report how much they differ.
%
% Run A: dt = 3e-4 s nominal, adaptive retries -> 17 accepted steps
%        (out_pure2dmac_dtsmall_14steps.mat)
% Run B: dt = 6e-4 s nominal, 7 accepted steps, no retries
%        (out_pure2dmac_dtlarge_7steps.mat)
%
% IMPORTANT: adaptive retries in Run A (starting around its step 12, where
% the gap tightens toward near-contact) mean step-index 2*kB does NOT
% reliably land on the same physical time as Run B's step kB beyond that
% point -- dt was repeatedly halved, so later steps advance less physical
% time per accepted step than the nominal 3e-4 s would suggest. This
% script therefore matches by ACTUAL physical time (nearest available
% step in Run A to each Run B time), not by a fixed index relationship,
% and reports both times explicitly so any mismatch is visible rather
% than silently assumed away.

clc;

outA = load_out('out_pure2dmac_dtsmall_14steps.mat');
outB = load_out('out_pure2dmac_dtlarge_7steps.mat');

nB = outB.stopStep;
fprintf('Run A (dt=3e-4 nominal, small): stopStep=%d, final t=%.6e s\n', outA.stopStep, outA.t(outA.stopStep));
fprintf('Run B (dt=6e-4 nominal, large): stopStep=%d, final t=%.6e s\n\n', nB, outB.t(nB));

fprintf('%-6s %-12s %-12s %-12s %-14s %-14s %-10s\n', ...
    'kB', 'tB [s]', 'tA [s]', '|dt| [s]', 'A (small dt)', 'B (large dt)', '%diff');

for kB = 1:nB
    tB = outB.t(kB);

    % Find the Run A step whose time is closest to tB (nearest match, not
    % an assumed index relationship).
    [dtMismatch, kA] = min(abs(outA.t(1:outA.stopStep) - tB));
    tA = outA.t(kA);

    flA = outA.fluidHist{kA};
    flB = outB.fluidHist{kB};

    pA = max(abs(flA.p(:)));
    pB = max(abs(flB.p(:)));
    pctP = abs(pA-pB) / max(abs(pB),1e-6) * 100;

    tauLA = max(abs(flA.tauL(:))); tauLB = max(abs(flB.tauL(:)));
    pctTauL = abs(tauLA-tauLB) / max(abs(tauLB),1e-6) * 100;

    tauEA = max(abs(flA.tauE(:))); tauEB = max(abs(flB.tauE(:)));
    pctTauE = abs(tauEA-tauEB) / max(abs(tauEB),1e-6) * 100;

    fprintf('%-6d %-12.4e %-12.4e %-12.2e (stepA=%d, stepB=%d)\n', kB, tB, tA, dtMismatch, kA, kB);
    fprintf('  max|p|    %-14.4f %-14.4f %-10.3f\n', pA, pB, pctP);
    fprintf('  max|tauL| %-14.4f %-14.4f %-10.3f\n', tauLA, tauLB, pctTauL);
    fprintf('  max|tauE| %-14.4f %-14.4f %-10.3f\n', tauEA, tauEB, pctTauE);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');

function out = load_out(file)
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
end
