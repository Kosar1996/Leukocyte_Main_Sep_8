%% CHECK_GAP_MISMATCH_IS_REAL
% Tests whether the mismatch blowup at steps 37-40 (as the gap shrinks
% toward ~0.05 um) is a genuine physical effect, or a measurement artifact
% from the epsFrac*gap sampling offset shrinking along with the gap.
%
% Method: recompute the mismatch at those same steps using several
% DIFFERENT epsFrac values. If the readings are robust (similar) across
% choices of epsFrac, that's evidence of a real effect. If the readings
% change dramatically depending on epsFrac, that points to a measurement
% sensitivity artifact rather than genuine physics.

clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

par = out.par;
parLmismatch = par;
if isfield(par, 'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par, 'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par, 'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

stepsToCheck = [30 34 36 37 38 39 40]; % includes some pre-blowup steps for context
epsFracValues = [0.05, 0.10, 0.20, 0.30];

fprintf('%6s | ', 'step');
for e = epsFracValues
    fprintf('epsFrac=%.2f(EnL) ', e);
end
fprintf('\n');

results = nan(numel(stepsToCheck), numel(epsFracValues));

for si = 1:numel(stepsToCheck)
    k = stepsToCheck(si);
    st = out.stateHist{k};
    fl = out.fluidHist{k};
    dtStep = out.dtHist(k);

    fprintf('%6d | ', k);
    for ei = 1:numel(epsFracValues)
        mismatchOpts = struct('nQuery', 41, 'epsFrac', epsFracValues(ei), 'trimFrac', 0.05);
        try
            [leuko, ~] = compute_interface_traction_mismatch( ...
                out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);
            results(si,ei) = leuko.stats.pctMaxEn;
            fprintf('%14.2f%% ', results(si,ei));
        catch ME
            fprintf('%14s ', 'FAILED');
        end
    end
    fprintf('\n');
end

fprintf('\n--- Robustness check: how much does EnL vary with epsFrac, at each step? ---\n');
for si = 1:numel(stepsToCheck)
    row = results(si,:);
    row = row(isfinite(row));
    if numel(row) < 2, continue; end
    spreadPct = 100*(max(row)-min(row))/max(mean(row),1e-9);
    fprintf('  step %2d: EnL ranges %.2f%% to %.2f%% across epsFrac choices -- spread = %.1f%% of the mean\n', ...
        stepsToCheck(si), min(row), max(row), spreadPct);
end
fprintf(['\nIf the spread stays small and consistent (similar %% at every step, including the\n' ...
    'pre-blowup steps like 30/34), the blowup is likely real physics, not a sampling artifact.\n' ...
    'If the spread grows sharply specifically at the small-gap steps (37-40), that is evidence\n' ...
    'the measurement itself becomes unstable as epsFrac*gap gets tiny.\n']);
fprintf('\nSMOKE_TEST_STATUS: OK\n');
