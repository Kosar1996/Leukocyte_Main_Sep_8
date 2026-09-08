%% COMPARE_EXTERIOR_ZPROFILE
% The actual intended test: compare pressure specifically AWAY from the
% near-contact gap (the "exterior" region the timescale analysis flagged),
% not just the domain-wide max|p| which is dominated by the gap spike.
%
% Uses the now-valid, bug-fixed saved outputs from both fine-step runs.
% Compares at matching physical time (nearest available step in each),
% across the full z-profile, and reports separately near the gap vs away
% from it based on each step's own gap(z) profile.

clc;

S1 = load('out_finestep_STEADY.mat');
S2 = load('out_finestep_UNSTEADY.mat');
outS = S1.outFine;
outU = S2.outFine;

fprintf('Steady:   stopStep=%d, t final=%.6e s\n', outS.stopStep, outS.t(outS.stopStep));
fprintf('Unsteady: stopStep=%d, t final=%.6e s\n\n', outU.stopStep, outU.t(outU.stopStep));

z = outS.z(:);

for kU = 1:outU.stopStep
    tU = outU.t(kU);
    [dtMismatch, kS] = min(abs(outS.t(1:outS.stopStep) - tU));
    tS = outS.t(kS);

    stU = outU.stateHist{kU};
    gapU = stU.deltaE(:) - stU.deltaL(:);
    [gapMin, iGapMin] = min(gapU);

    pU = outU.fluidHist{kU}.p(:);
    pS = outS.fluidHist{kS}.p(:);

    % "Near gap": within the closest 20% of z-points to the minimum-gap
    % location. "Exterior": the rest.
    nZ = numel(z);
    distToMin = abs((1:nZ)' - iGapMin);
    nearMask = distToMin <= round(0.1*nZ);
    farMask = ~nearMask;

    pctDiff = abs(pU - pS) ./ max(abs(pS), 1e-6) * 100;

    fprintf('kU=%d (t=%.4e s, matched kS=%d t=%.4e s, mismatch=%.2e s):\n', ...
        kU, tU, kS, tS, dtMismatch);
    fprintf('  gapMin=%.4e m at z=%.4e um\n', gapMin, z(iGapMin)*1e6);
    fprintf('  near-gap  (%d pts): max|p| steady=%.4f unsteady=%.4f, median %%diff=%.3f, max %%diff=%.3f\n', ...
        nnz(nearMask), max(abs(pS(nearMask))), max(abs(pU(nearMask))), ...
        median(pctDiff(nearMask)), max(pctDiff(nearMask)));
    fprintf('  exterior  (%d pts): max|p| steady=%.4f unsteady=%.4f, median %%diff=%.3f, max %%diff=%.3f\n', ...
        nnz(farMask), max(abs(pS(farMask))), max(abs(pU(farMask))), ...
        median(pctDiff(farMask)), max(pctDiff(farMask)));
    fprintf('\n');
end

fprintf('SMOKE_TEST_STATUS: OK\n');
