%% CHECK_SHAPE_DEGRADATION_OVER_TIME

% and then start to deviate at later time points, using the ALREADY-SAVED
% 40-step 2D run -- no new simulation needed.
%
% Pulls per-step diagnostics already logged during the run (out.diagHist)
% across all accepted steps: mesh quality (minJE/minJL -- Jacobian minima,
% a measure of element distortion; approaching 0 means a element is
% collapsing/inverting, i.e. genuinely nonphysical), mass conservation
% (volumeChangeRel, fluxJump), and the minimum gap each step.

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

nSteps = out.stopStep;
t = out.t(1:nSteps);

minJE = nan(nSteps,1);
minJL = nan(nSteps,1);
gapMin = nan(nSteps,1);
volChangeRel = nan(nSteps,1);
fluxJump = nan(nSteps,1);
fluidResInf = nan(nSteps,1);
retries = nan(nSteps,1);

for k = 1:nSteps
    d = out.diagHist{k};
    if isempty(d), continue; end
    if isfield(d,'minJE'),           minJE(k) = d.minJE; end
    if isfield(d,'minJL'),           minJL(k) = d.minJL; end
    if isfield(d,'gapMin'),          gapMin(k) = d.gapMin; end
    if isfield(d,'volumeChangeRel'), volChangeRel(k) = d.volumeChangeRel; end
    if isfield(d,'fluxJump'),        fluxJump(k) = d.fluxJump; end
    if isfield(d,'fluidResidualInf'),fluidResInf(k) = d.fluidResidualInf; end
    if isfield(d,'retries'),         retries(k) = d.retries; end
end

fprintf('%6s %10s | %10s %10s | %10s | %10s %10s | %8s\n', ...
    'step', 't[s]', 'minJE', 'minJL', 'gapMin[m]', 'volChgRel', 'fluxJump', 'retries');
for k = 1:nSteps
    fprintf('%6d %10.4e | %10.4f %10.4f | %10.4e | %10.4e %10.4e | %8d\n', ...
        k, t(k), minJE(k), minJL(k), gapMin(k), volChangeRel(k), fluxJump(k), retries(k));
end

fprintf('\n--- Trend check ---\n');
fprintf('minJE:  first=%.4f  last=%.4f  min-over-run=%.4f (at step %d)\n', ...
    minJE(1), minJE(end), min(minJE), find(minJE==min(minJE),1));
fprintf('minJL:  first=%.4f  last=%.4f  min-over-run=%.4f (at step %d)\n', ...
    minJL(1), minJL(end), min(minJL), find(minJL==min(minJL),1));
fprintf('gapMin: first=%.4e  last=%.4e  min-over-run=%.4e (at step %d)\n', ...
    gapMin(1), gapMin(end), min(gapMin), find(gapMin==min(gapMin),1));
fprintf('volChangeRel: max-over-run=%.4e (at step %d)\n', ...
    max(abs(volChangeRel)), find(abs(volChangeRel)==max(abs(volChangeRel)),1));
fprintf('total retries across run: %d\n', sum(retries));
fprintf(['\nA Jacobian trending toward 0 (element inversion), gapMin trending toward 0 ' ...
    '(solids about to touch/overlap), or volChangeRel growing over time, would all indicate ' ...
    'genuine degradation as the run progresses -- not just noise.\n']);

fig = figure('Position',[100 100 900 700],'Color','w');
tl = tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');

ax1 = nexttile(tl);
plot(ax1, 1:nSteps, minJE, 'o-', 'DisplayName','minJE'); hold(ax1,'on');
plot(ax1, 1:nSteps, minJL, 's-', 'DisplayName','minJL');
xlabel(ax1,'step'); ylabel(ax1,'Jacobian min'); legend(ax1,'Location','best'); grid(ax1,'on');
title(ax1,'Mesh quality (lower = more distorted)');

ax2 = nexttile(tl);
plot(ax2, 1:nSteps, gapMin*1e6, 'o-');
xlabel(ax2,'step'); ylabel(ax2,'min gap [\mum]'); grid(ax2,'on');
title(ax2,'Minimum leukocyte-endothelium gap');

ax3 = nexttile(tl);
semilogy(ax3, 1:nSteps, abs(volChangeRel), 'o-');
xlabel(ax3,'step'); ylabel(ax3,'|volume change rel|'); grid(ax3,'on');
title(ax3,'Mass conservation drift');

ax4 = nexttile(tl);
plot(ax4, 1:nSteps, retries, 'o-');
xlabel(ax4,'step'); ylabel(ax4,'adaptive retries'); grid(ax4,'on');
title(ax4,'Retries per step (spikes = solver struggling)');

sgtitle(fig, '2D solver, 40 steps -- shape/quality trend over time');
outPng = fullfile(fileparts(mfilename('fullpath')), 'shape_degradation_2D.png');
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('\nSaved: %s\nSMOKE_TEST_STATUS: OK\n', outPng);
