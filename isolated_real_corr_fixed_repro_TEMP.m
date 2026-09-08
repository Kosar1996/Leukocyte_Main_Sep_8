%% Reproduce the real production traction-correction call with BOTH
% false-convergence tolerances stripped (solidFallbackAbsTol, already
% absent from this dump; solidAbsTol=1e-10, just fixed) -- to see whether
% the Newton solve genuinely converges on this REAL, warm-started,
% modest-magnitude (median ~166 Pa, not the synthetic -500 Pa cold-start
% test) production traction, or whether it still stalls once the false
% "convergence" exit is removed.

clc; clear all;
cd('/Users/kosarsafari/Desktop/Project_1/code/leukocyte-main');

S = load('/tmp/corr_call_dump_-6000_1.mat');

parFixed = S.parCorrE;
if isfield(parFixed,'solidAbsTol'), parFixed = rmfield(parFixed,'solidAbsTol'); end
if isfield(parFixed,'solidFallbackAbsTol'), parFixed = rmfield(parFixed,'solidFallbackAbsTol'); end
parFixed.debugVerbose = true;

fprintf('=== reproducing with solidAbsTol AND solidFallbackAbsTol stripped ===\n');
try
    uEcorr = solve_finite_def_solid(S.meshE, S.uOldE, S.tractionE, S.interfaceE, S.baseE, S.supportE, parFixed, S.uInitE);
    fprintf('\nCONVERGED genuinely.\n');
    fprintf('norm(uEcorr - uInitE) = %.6e\n', norm(uEcorr - S.uInitE));
    fprintf('max|uEcorr - uInitE| = %.6e m\n', max(abs(uEcorr - S.uInitE)));
catch ME
    fprintf('\nFAILED: %s\n', ME.message);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
