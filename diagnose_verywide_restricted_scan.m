clc;
Sw = load('solid_endo_P300_wide.mat');
Sv = load('solid_endo_P300_verywide.mat');
SL = load('solid_leu_P600.mat');

% restrict to the leukocyte's OWN true physical z-extent (-2,6)um, where
% deltaL is real data, not extrapolated
zScan = linspace(-2e-6, 6e-6, 41).';

dEw = interp1(Sw.z(:), Sw.deltaE_pre(:), zScan, 'linear', 'extrap');
dEv = interp1(Sv.z(:), Sv.deltaE_pre(:), zScan, 'linear', 'extrap');
dL  = interp1(SL.z(:), SL.deltaL_pre(:), zScan, 'linear', 'extrap');

fprintf('%8s | %10s %10s | %10s %10s\n', 'z[um]', 'deltaE_w', 'deltaE_v', 'gap_w[nm]', 'gap_v[nm]');
for i = 1:numel(zScan)
    gw = (dEw(i)-dL(i))*1e9;
    gv = (dEv(i)-dL(i))*1e9;
    fprintf('%8.3f | %10.4f %10.4f | %10.2f %10.2f\n', zScan(i)*1e6, dEw(i)*1e6, dEv(i)*1e6, gw, gv);
end
[minGw, iw] = min(dEw-dL);
[minGv, iv] = min(dEv-dL);
fprintf('\nmin gap (wide, restricted [-2,6]um)     = %.2f nm at z=%.3f um\n', minGw*1e9, zScan(iw)*1e6);
fprintf('min gap (verywide, restricted [-2,6]um)  = %.2f nm at z=%.3f um\n', minGv*1e9, zScan(iv)*1e6);
fprintf('SMOKE_TEST_STATUS: OK\n');
