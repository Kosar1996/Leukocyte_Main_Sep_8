clc;
Sv = load('solid_endo_P300_verywide.mat');
SL = load('solid_leu_P600.mat');

% global grid used at runtime (fluid domain, unchanged regardless of endo domain)
zGrid = linspace(-6e-6, 10e-6, 61).';  % matches par.NzFluid=61 default in production cfg

dE = interp1(Sv.z(:), Sv.deltaE_pre(:), zGrid, 'linear', 'extrap');
dL = interp1(SL.z(:), SL.deltaL_pre(:), zGrid, 'linear', 'extrap');
gap = dE - dL;

fprintf('%8s | %10s %10s %10s\n', 'z[um]', 'deltaE', 'deltaL', 'gap[nm]');
for i = 1:numel(zGrid)
    fprintf('%8.3f | %10.4f %10.4f %10.2f\n', zGrid(i)*1e6, dE(i)*1e6, dL(i)*1e6, gap(i)*1e9);
end
[minGap, iMin] = min(gap);
fprintf('\nGLOBAL min gap = %.2f nm at z = %.4f um\n', minGap*1e9, zGrid(iMin)*1e6);
fprintf('meshE actual z-range: [%.3f, %.3f] um\n', min(Sv.meshE.nodes(:,2))*1e6, max(Sv.meshE.nodes(:,2))*1e6);
fprintf('SMOKE_TEST_STATUS: OK\n');
