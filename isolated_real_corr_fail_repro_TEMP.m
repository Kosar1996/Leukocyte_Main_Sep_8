%% Load the REAL, production-captured inputs from the first traction-
% correction failure (real coupled run, step 1, correction pass 1) and
% reproduce/instrument the Newton solve directly -- no synthetic test
% assumptions, exact production state.

clc; clear all;
cd('/Users/kosarsafari/Desktop/Project_1/code/leukocyte-main');

S = load('/tmp/corr_call_dump_-6000_1.mat');

fprintf('=== traction magnitude summary ===\n');
fprintf('normal traction: min=%.4f max=%.4f median=%.4f Pa\n', ...
    min(S.tractionE.normal), max(S.tractionE.normal), median(S.tractionE.normal));
fprintf('tangent traction: min=%.4f max=%.4f median=%.4f Pa\n', ...
    min(S.tractionE.tangent), max(S.tractionE.tangent), median(S.tractionE.tangent));

fprintf('\n=== initial guess vs old state ===\n');
duInit = S.uInitE - S.uOldE;
fprintf('norm(uInitE - uOldE) = %.6e   (nonzero = warm start, as expected)\n', norm(duInit));
fprintf('max|uInitE - uOldE| = %.6e m\n', max(abs(duInit)));

parTest = S.parCorrE;
parTest.debugVerbose = true;
parTest.solidLoadSteps = 1;  % see the raw single-shot behavior first, exactly as production calls it via the wrapper default... actually check what parCorrE.solidLoadSteps is
fprintf('\nparCorrE.solidLoadSteps (as passed from production) = %s\n', mat2str(getfield_default(S.parCorrE,'solidLoadSteps',NaN)));
fprintf('parCorrE.newtonMaxItSolid = %d, newtonTolSolid = %.3e\n', S.parCorrE.newtonMaxItSolid, S.parCorrE.newtonTolSolid);
fprintf('parCorrE.trustU0 = %.3e, trustUMin = %.3e, trustUMax = %.3e\n', S.parCorrE.trustU0, S.parCorrE.trustUMin, S.parCorrE.trustUMax);

fprintf('\n=== reproducing the EXACT production call (with debugVerbose) ===\n');
parTest2 = S.parCorrE;
parTest2.debugVerbose = true;
try
    uEcorr = solve_finite_def_solid(S.meshE, S.uOldE, S.tractionE, S.interfaceE, S.baseE, S.supportE, parTest2, S.uInitE);
    fprintf('CONVERGED (unexpected given production log)\n');
catch ME
    fprintf('FAILED as expected: %s\n', ME.message);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');

function v = getfield_default(s, f, d)
if isfield(s,f), v = s.(f); else, v = d; end
end
