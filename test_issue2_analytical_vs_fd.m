%% TEST_ISSUE2_ANALYTICAL_VS_FD
% Before proposing anything as a "fix," verify whether the exploding
% d(tauE)/d(re) is a genuine mathematical property of the exact
% no-slip annular Couette-Poiseuille solution, or an artifact of the
% finite-difference sensitivity computation itself (e.g. a poorly-scaled
% perturbation epsR near a highly nonlinear point).
%
% Hand-derived analytical formula (no-slip case, ll=le=0):
%   C1(b) = N(b)/D(b), N(b) = dUw - (G/4)*(a^2-b^2), D(b) = log(a)-log(b)
%   dC1/db = [N'(b)*D(b) - N(b)*D'(b)] / D(b)^2
%          = [(G/2)*b*D(b) + N(b)/b] / D(b)^2
%   tauE = mu*(G*b/2 + C1/b)
%   d(tauE)/db = mu*(G/2 + (dC1/db)/b - C1/b^2)
%
% Compares this closed-form derivative directly against the existing,
% unmodified local_flux_shear_sensitivities_fd_full.m's FD result, at
% the same range of gaps. Purely diagnostic; does not modify any code.

clc;

S = load('out_pure2dmac_dtlarge_7steps.mat');
par = S.out.par;
par.NzFluid = 2;
mu = par.mu;

r0 = 3.26e-6;
G = 1e3;
UwL0 = 1e-7;
UwE0 = -1e-7;
dUw = UwL0 - UwE0;

gaps = [1e-6, 1e-7, 1e-8, 1e-9, 1.01e-9];

fprintf('%12s %16s %16s %14s\n', 'gap (m)', 'analytical', 'FD (code)', 'relative diff');
for k = 1:numel(gaps)
    h = gaps(k);
    a = r0 - h/2;
    b = r0 + h/2;

    % Analytical
    N = dUw - (G/4)*(a^2 - b^2);
    D = log(a) - log(b);
    Np = (G/2)*b;
    Dp = -1/b;
    dC1db = (Np*D - N*Dp) / D^2;
    C1 = N/D;
    dtauEdb_analytical = mu*(G/2 + dC1db/b - C1/b^2);

    % FD, via the existing unmodified code
    z = [0; 1];
    p = [0; -G*mu*(z(2)-z(1))];
    rl = [a; a];
    re = [b; b];
    UwL = [UwL0; UwL0];
    UwE = [UwE0; UwE0];
    fs = local_flux_shear_sensitivities_fd_full(z, p, rl, re, UwL, UwE, par);
    dtauEdb_fd = full(fs.dtauEdre(1,1));

    relDiff = abs(dtauEdb_analytical - dtauEdb_fd) / max(abs(dtauEdb_analytical), 1e-300);

    fprintf('%12.4e %16.6e %16.6e %14.4e\n', h, dtauEdb_analytical, dtauEdb_fd, relDiff);
end

fprintf('\nIf relative diff stays small (<1e-3) across all gaps, the FD result matches\n');
fprintf('the exact analytical derivative -- confirming this is genuine math, not an FD bug.\n');
fprintf('\nSMOKE_TEST_STATUS: OK\n');
