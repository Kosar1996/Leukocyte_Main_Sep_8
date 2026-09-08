%% DEMO_TRACTION_CORRECTION_FEEDBACK
% Demonstrates apply_bodyfitted_MAC_traction_correction_feedback.m against
% a saved run, and plots the %mismatch-vs-iteration convergence history.
%
% Runs TWO scenarios so both behaviors are visible in one place:
%   Scenario A: the correction loop applied to an already-accepted,
%     converged step. Expected: it should recognize the interface is
%     already at (or very near) the correction map's fixed point --
%     mismatch flat, no further displacement.
%   Scenario B: the same fluid traction, but paired with the PREVIOUS
%     step's (not-yet-corrected-for-this-traction) solid shape, so there
%     is genuine room to improve. Expected: mismatch should actually move
%     iteration to iteration, converging (or reporting that it did not,
%     honestly) instead of silently stopping after a fixed pass count.
%
% Edit savedRunFile / stepToUse below to point at a different run/step.

clc; close all;
cd(fileparts(mfilename('fullpath')));

savedRunFile = 'out_1D_t10_for_review.mat';
S = load(savedRunFile);
runFieldName = fieldnames(S);
out = S.(runFieldName{1});   % out1D or out2D, whichever variable the file stores

SE = load(out.par.endotheliumPrestressFile);
SL = load(out.par.leukocytePrestressFile);
baseE = SE.baseE;
baseL = SL.baseL;

% Leukocyte-remapped material parameters, same remap used throughout this
% codebase (leukocyte_solid_parameters.m is a private local function
% inside softlube_run_case_global_coupled.m, so it's inlined here).
parL = out.par;
if isfield(out.par, 'EL'), parL.Ee = out.par.EL; end
if isfield(out.par, 'nuL'), parL.nuE = out.par.nuL; end
if isfield(out.par, 'GL')
    parL.Ge = out.par.GL;
elseif isfield(parL, 'Ee') && isfield(parL, 'nuE')
    parL.Ge = parL.Ee/(2*(1+parL.nuE));
end
if isfield(out.par, 'KL')
    parL.Ke = out.par.KL;
elseif isfield(parL, 'Ee') && isfield(parL, 'nuE')
    parL.Ke = parL.Ee/(3*(1-2*parL.nuE));
end
if isfield(out.par, 'etaL'), parL.etaE = out.par.etaL; end
if isfield(out.par, 'useViscoelasticLeukocyte')
    parL.useViscoelasticEndothelium = out.par.useViscoelasticLeukocyte;
end

k = out.stopStep;
z = out.z;
par = out.par;
par.dt = out.dtHist(k);

opts = struct();
opts.mismatchThresholdPct = 10;
opts.statForConvergence = 'p90';
opts.maxIterations = 8;
opts.verbose = true;

%% Scenario A: already-accepted, self-consistent state
fprintf('\n===================== SCENARIO A: already-converged step =====================\n');
old_A = out.stateHist{k-1};
state_A = out.stateHist{k};
fluid_A = out.fluidHist{k};

[~, ~, okA, ~, convA] = apply_bodyfitted_MAC_traction_correction_feedback( ...
    z, old_A, state_A, fluid_A, out.meshE, out.interfaceE, baseE, ...
    out.meshL, out.interfaceL, baseL, parL, par, opts);

%% Scenario B: previous step's solid shape under the later step's fluid traction
% (deliberately mismatched pairing, so there's real room for the loop to work)
fprintf('\n===================== SCENARIO B: deliberately non-equilibrated =====================\n');
old_B = out.stateHist{k-1};
state_B = out.stateHist{k-1};   % start from the OLDER solid shape...
fluid_B = out.fluidHist{k};      % ...under the NEWER step's fluid traction

[~, ~, okB, ~, convB] = apply_bodyfitted_MAC_traction_correction_feedback( ...
    z, old_B, state_B, fluid_B, out.meshE, out.interfaceE, baseE, ...
    out.meshL, out.interfaceL, baseL, parL, par, opts);

%% Plot both convergence histories
figure('Name', 'Traction correction feedback: convergence history', 'Position', [100 100 900 700]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot_history(convA.history, opts.mismatchThresholdPct);
title(sprintf('Scenario A (already converged): converged=%d, final worst=%.2f%%', ...
    convA.converged, convA.finalWorstPct));

nexttile;
plot_history(convB.history, opts.mismatchThresholdPct);
title(sprintf('Scenario B (room to improve): converged=%d, final worst=%.2f%%', ...
    convB.converged, convB.finalWorstPct));

sgtitle('Feedback-loop traction correction: %mismatch vs iteration (worst of Leuko en/et, Endo en/et)');

function plot_history(history, thresholdPct)
iters = [history.iter];
worst = [history.worstPct];
enL = [history.pctEnL]; etL = [history.pctEtL];
enE = [history.pctEnE]; etE = [history.pctEtE];

semilogy(iters, worst, 'k-o', 'LineWidth', 2, 'DisplayName', 'worst of all 4'); hold on;
semilogy(iters, enL, '--', 'DisplayName', 'Leuko e_n');
semilogy(iters, etL, '--', 'DisplayName', 'Leuko e_t');
semilogy(iters, enE, ':', 'DisplayName', 'Endo e_n');
semilogy(iters, etE, ':', 'DisplayName', 'Endo e_t');
yline(thresholdPct, 'r-', 'LineWidth', 1.5, 'DisplayName', 'threshold');
xlabel('correction pass'); ylabel('%mismatch (log scale)');
legend('Location','best'); grid on;
end
