clc; close all;
cd(fileparts(mfilename('fullpath')));

generate_prestress_endothelium_wide(-4e-6, 8e-6, 'solid_endo_P300_verywide.mat');

S = load('solid_endo_P300_verywide.mat');
if ~isreal(S.meshE.nodes) || ~isreal(S.uE_pre)
    error('Complex values detected in very-wide prestress -- do not proceed.');
end
fprintf('Very-wide prestress verified real.\n');

cfg = build_cfg_full2D_pressure2('nSteps', 40, 'useHybridGap1DExterior2DFluid', false);
cfg.geometry.endotheliumPrestressFile = fullfile(fileparts(mfilename('fullpath')), 'solid_endo_P300_verywide.mat');

fprintf('Running 2D, 40 steps, VERY WIDE endothelium domain (z=[-4,8]um) ...\n');
out = softlube_run_case_global_coupled(cfg);

save('out_2D_verywide_domain_40step.mat', 'out', '-v7.3');
fprintf('Saved out_2D_verywide_domain_40step.mat, stopStep=%d\n', out.stopStep);
fprintf('SMOKE_TEST_STATUS: OK\n');
