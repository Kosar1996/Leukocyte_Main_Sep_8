clc;
Sorig = load('solid_endo_P300.mat');
Swide = load('solid_endo_P300_wide.mat');

zCommon = linspace(0, 4, 30)' * 1e-6;

zOrig = Sorig.z(:);
zWide = Swide.z(:);
dOrig = interp1(zOrig, Sorig.deltaE_pre(:), zCommon, 'linear', 'extrap');
dWide = interp1(zWide, Swide.deltaE_pre(:), zCommon, 'linear', 'extrap');

fprintf('%8s | %10s %10s %10s\n', 'z[um]', 'orig[um]', 'wide[um]', 'diff[nm]');
for i = 1:numel(zCommon)
    fprintf('%8.3f | %10.4f %10.4f %10.2f\n', zCommon(i)*1e6, dOrig(i)*1e6, dWide(i)*1e6, (dWide(i)-dOrig(i))*1e9);
end
fprintf('\nmax abs diff over [0,4]um: %.3f nm\n', max(abs(dWide-dOrig))*1e9);
fprintf('max abs diff over [3.5,3.9]um (near original contact zone): %.3f nm\n', ...
    max(abs(dWide(zCommon>=3.5e-6 & zCommon<=3.9e-6) - dOrig(zCommon>=3.5e-6 & zCommon<=3.9e-6)))*1e9);
fprintf('SMOKE_TEST_STATUS: OK\n');
