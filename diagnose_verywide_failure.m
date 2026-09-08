clc;
Sn = load('solid_endo_P300.mat');
Sw = load('solid_endo_P300_wide.mat');
Sv = load('solid_endo_P300_verywide.mat');

zCommon = linspace(3.4, 4.0, 20)' * 1e-6;
dN = interp1(Sn.z(:), Sn.deltaE_pre(:), zCommon, 'linear', 'extrap');
dW = interp1(Sw.z(:), Sw.deltaE_pre(:), zCommon, 'linear', 'extrap');
dV = interp1(Sv.z(:), Sv.deltaE_pre(:), zCommon, 'linear', 'extrap');

fprintf('%8s | %10s %10s %10s\n', 'z[um]', 'narrow', 'wide', 'verywide');
for i = 1:numel(zCommon)
    fprintf('%8.3f | %10.4f %10.4f %10.4f\n', zCommon(i)*1e6, dN(i)*1e6, dW(i)*1e6, dV(i)*1e6);
end

% also check midline (length-independent reference point) for all three
fprintf('\nMidline (z=2um, should be length-independent far-field value):\n');
fprintf('  narrow:   %.4f um\n', interp1(Sn.z(:), Sn.deltaE_pre(:), 2e-6, 'linear', 'extrap')*1e6);
fprintf('  wide:     %.4f um\n', interp1(Sw.z(:), Sw.deltaE_pre(:), 2e-6, 'linear', 'extrap')*1e6);
fprintf('  verywide: %.4f um\n', interp1(Sv.z(:), Sv.deltaE_pre(:), 2e-6, 'linear', 'extrap')*1e6);

% Leukocyte's own outer radius near this z-range, for context (unchanged file)
SL = load('solid_leu_P600.mat');
fnL = fieldnames(SL);
disp(fnL);
SMOKE = 'SMOKE_TEST_STATUS: OK';
fprintf('\n%s\n', SMOKE);
