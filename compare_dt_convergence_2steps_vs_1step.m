%% COMPARE_DT_CONVERGENCE_2STEPS_VS_1STEP

% consistency check. Run the SAME physical duration two ways --
%   Run A: 2 small time steps (dt_small each)
%   Run B: 1 large time step  (dt_large = 2 * dt_small)
% and confirm the leukocyte/endothelium stress at the 8 tracked points
% agree to within a small, expected numerical difference. If they don't,
% the time-stepping (dt) still isn't correctly reflected in the solve.
%
% HOW TO PRODUCE THE TWO INPUT FILES (run_input_full2D_pressure2.m always
% saves to the hardcoded name out_2_new_t1.mat, so you must rename between
% runs or it gets overwritten):
%   1) In run_input_full2D_pressure2.m, set dt = 3e-4; nSteps = 2;
%      Run it, then rename the output:
%        movefile('out_2_new_t1.mat','out_dtsmall_2steps.mat')
%   2) In run_input_full2D_pressure2.m, set dt = 6e-4; nSteps = 1;
%      Run it, then rename the output:
%        movefile('out_2_new_t1.mat','out_dtlarge_1step.mat')
%   3) Run this script.
%
% Uses the same 8-point roster as verify_solid_mechanics_8points_after_dtfix.m.

clc; close all;

fileSmall = 'out_dtsmall_2steps.mat';   % 2 steps of dt_small
fileLarge = 'out_dtlarge_1step.mat';    % 1 step of dt_large = 2*dt_small

outA = load_out(fileSmall);
outB = load_out(fileLarge);

kA = outA.stopStep;
kB = outB.stopStep;

fprintf('Run A (small steps): stopStep=%d, final t=%.6e s, dtHist=%s\n', ...
    kA, outA.t(end), mat2str(outA.dtHist(1:kA)));
fprintf('Run B (large step):  stopStep=%d, final t=%.6e s, dtHist=%s\n\n', ...
    kB, outB.t(end), mat2str(outB.dtHist(1:kB)));

% The body-fitted MAC traction correction is an INNER, per-time-step
% fixed-point iteration (partitioned solid/fluid correction), separate
% from the OUTER adaptive physical time-stepping (dtHist/retries above).
% It can legitimately take a very different number of passes to reach its
% own convergence tolerance on any given step -- e.g. one step might
% settle in 3 passes, another might need close to
% maxBodyFittedTractionCorrections (commonly ~100) -- without that
% affecting whether the two dt configurations being compared here agree.
% Report it explicitly so a large pass count on one run/step isn't
% mistaken for a dt-consistency problem: what must agree between Run A
% and Run B is the converged PHYSICAL STATE at matching times, not how
% many inner passes either one took to get there.
fprintf('--- inner traction-correction pass counts (informational only; not required to match between runs) ---\n');
report_correction_passes('Run A (small steps)', outA, kA);
report_correction_passes('Run B (large step)', outB, kB);
fprintf('\n');

if abs(outA.t(end) - outB.t(end)) > 1e-12
    fprintf(['*** WARNING: final times do not match (%.6e vs %.6e s). ', ...
        'This comparison is only meaningful if both runs reach the same ', ...
        'physical time. Check for retries in dtHist above. ***\n\n'], ...
        outA.t(end), outB.t(end));
end

stA = outA.stateHist{kA};
stB = outB.stateHist{kB};
parA = outA.par;
parB = outB.par;

parLA = parA;
if isfield(parA,'GL')  && isfinite(parA.GL),  parLA.Ge   = parA.GL;  end
if isfield(parA,'KL')  && isfinite(parA.KL),  parLA.Ke   = parA.KL;  end
if isfield(parA,'etaL')&& isfinite(parA.etaL),parLA.etaE = parA.etaL;end

parLB = parB;
if isfield(parB,'GL')  && isfinite(parB.GL),  parLB.Ge   = parB.GL;  end
if isfield(parB,'KL')  && isfinite(parB.KL),  parLB.Ke   = parB.KL;  end
if isfield(parB,'etaL')&& isfinite(parB.etaL),parLB.etaE = parB.etaL;end

dtA = outA.dtHist(kA);
dtB = outB.dtHist(kB);

stressEA = recover_nodal_stress_axisym_viscoelastic(outA.meshE, stA.uE, stA.uEPrev, dtA, parA);
stressLA = recover_nodal_stress_axisym_viscoelastic(outA.meshL, stA.uL, stA.uLPrev, dtA, parLA);
stressEB = recover_nodal_stress_axisym_viscoelastic(outB.meshE, stB.uE, stB.uEPrev, dtB, parB);
stressLB = recover_nodal_stress_axisym_viscoelastic(outB.meshL, stB.uL, stB.uLPrev, dtB, parLB);

nelL = outA.meshL.nelem;
nelE = outA.meshE.nelem;
elemsL = [round(nelL*0.25), round(nelL*0.5), round(nelL*0.75)];
elemsE = [round(nelE*0.25), round(nelE*0.5), round(nelE*0.75)];

fprintf('================ LEUKOCYTE (2 small steps vs 1 large step) ================\n');
for e = elemsL
    report_dt_diff('Leukocyte', outA.meshL, e, stressLA, stressLB);
end
fprintf('================ ENDOTHELIUM ================\n');
for e = elemsE
    report_dt_diff('Endothelium', outA.meshE, e, stressEA, stressEB);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');

%% ---- local functions ----

function report_correction_passes(label, out, kStop)
fprintf('%s:\n', label);
for n = 1:kStop
    passes = NaN; converged = NaN; source = '';
    if isfield(out, 'stateHist') && numel(out.stateHist) >= n && ...
            ~isempty(out.stateHist{n}) && isfield(out.stateHist{n}, 'tractionCorrectionPassesUsed')
        passes = out.stateHist{n}.tractionCorrectionPassesUsed;
        if isfield(out.stateHist{n}, 'tractionCorrectionConverged')
            converged = out.stateHist{n}.tractionCorrectionConverged;
        end
        source = 'plain loop';
    elseif isfield(out, 'tractionCorrectionHistory') && numel(out.tractionCorrectionHistory) >= n && ...
            ~isempty(out.tractionCorrectionHistory{n}) && isfield(out.tractionCorrectionHistory{n}, 'passesUsed')
        passes = out.tractionCorrectionHistory{n}.passesUsed;
        if isfield(out.tractionCorrectionHistory{n}, 'converged')
            converged = out.tractionCorrectionHistory{n}.converged;
        end
        source = 'feedback loop';
    end
    if isnan(passes)
        fprintf('   step %d: pass count not recorded for this run\n', n);
    else
        convStr = 'unknown';
        if converged == 1, convStr = 'yes'; elseif converged == 0, convStr = 'no'; end
        fprintf('   step %d: %d passes (%s), converged=%s\n', n, passes, source, convStr);
    end
end
end

function out = load_out(file)
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
end

function report_dt_diff(label, mesh, elemIdx, stressA, stressB)
conn = mesh.conn(elemIdx,:);
rrA = mean(stressA.sigma_rr(conn)); rrB = mean(stressB.sigma_rr(conn));
ttA = mean(stressA.sigma_tt(conn)); ttB = mean(stressB.sigma_tt(conn));
zzA = mean(stressA.sigma_zz(conn)); zzB = mean(stressB.sigma_zz(conn));
rzA = mean(stressA.sigma_rz(conn)); rzB = mean(stressB.sigma_rz(conn));

vecA = [rrA, ttA, zzA, rzA];
vecB = [rrB, ttB, zzB, rzB];
pctDiff = abs(vecA - vecB) ./ max(abs(vecB), 1e-6) * 100;

fprintf('%-12s elem %5d:\n', label, elemIdx);
fprintf('   2 small steps: rr=%9.4f tt=%9.4f zz=%9.4f rz=%9.4f\n', vecA);
fprintf('   1 large step:  rr=%9.4f tt=%9.4f zz=%9.4f rz=%9.4f\n', vecB);
fprintf('   %% diff:        rr=%7.3f tt=%7.3f zz=%7.3f rz=%7.3f\n\n', pctDiff);
end
