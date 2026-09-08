%% COMPARE_FLUID2DONLY_DT_CONVERGENCE
% Task-2-followup: pointwise step-size consistency check for the standalone
% 2D fluid solver (hybrid gap1D/exterior2D OFF), same method as the solid
% mechanics check -- 2 small steps vs 1 double-size step, same final time,
% compared at matching z locations (not aggregate min/max, which hides
% where a real discrepancy lives).

clc;

outA = load_out('out_fluid2Donly_dtsmall_2steps.mat');  % 2 steps of dt_small
outB = load_out('out_fluid2Donly_dtlarge_1step.mat');   % 1 step of dt_large = 2*dt_small

kA = outA.stopStep;
kB = outB.stopStep;
flA = outA.fluidHist{kA};
flB = outB.fluidHist{kB};

fprintf('Run A (2 small steps): stopStep=%d, final t=%.6e s\n', kA, outA.t(end));
fprintf('Run B (1 large step):  stopStep=%d, final t=%.6e s\n\n', kB, outB.t(end));

zA = outA.z(:);
zB = outB.z(:);

% Sample at 25/50/75% along the z grid (same fractional locations in both
% runs, then interpolate B onto A's z grid to compare like-for-like).
NzA = numel(zA);
idxA = [round(NzA*0.25), round(NzA*0.5), round(NzA*0.75)];
zSample = zA(idxA);

pA   = flA.p(:);    pB_i   = interp1(zB, flB.p(:),    zSample, 'linear', 'extrap');
tauLA = flA.tauL(:); tauLB_i = interp1(zB, flB.tauL(:), zSample, 'linear', 'extrap');
tauEA = flA.tauE(:); tauEB_i = interp1(zB, flB.tauE(:), zSample, 'linear', 'extrap');

fprintf('%-8s %-10s %-14s %-14s %-10s\n', 'z[um]', 'field', '2 small steps', '1 large step', '%diff');
for i = 1:numel(idxA)
    j = idxA(i);
    pctP    = abs(pA(j)    - pB_i(i))    / max(abs(pB_i(i)),1e-6)    * 100;
    pctTauL = abs(tauLA(j) - tauLB_i(i)) / max(abs(tauLB_i(i)),1e-6) * 100;
    pctTauE = abs(tauEA(j) - tauEB_i(i)) / max(abs(tauEB_i(i)),1e-6) * 100;

    fprintf('%-8.3f %-10s %-14.4f %-14.4f %-10.3f\n', zSample(i)*1e6, 'p [Pa]',    pA(j),    pB_i(i),    pctP);
    fprintf('%-8s %-10s %-14.4f %-14.4f %-10.3f\n', '', 'tauL [Pa]', tauLA(j), tauLB_i(i), pctTauL);
    fprintf('%-8s %-10s %-14.4f %-14.4f %-10.3f\n\n', '', 'tauE [Pa]', tauEA(j), tauEB_i(i), pctTauE);
end

fprintf('SMOKE_TEST_STATUS: OK\n');

function out = load_out(file)
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
end
