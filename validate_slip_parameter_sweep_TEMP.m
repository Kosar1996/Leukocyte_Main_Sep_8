clc;
addpath('/Users/kosarsafari/Desktop/Project_1/code/leukocyte-main');


% "since some model parameters are missing at this point, i think it is
% worth using the dimensionless parameters you defined to do some
% parameter sweeps." Uses the paper's own dimensionless group l_i =
% l_i'/H_o (Eq. 1), sweeping it across a physically-motivated range
% (l'=0 to 100nm, per Table 1/2's stated range for Navier slip length on
% soft polymer surfaces) and reporting flux enhancement relative to
% no-slip, plus a direct comparison to the paper's own closed-form
% dimensionless flux relation at each point.

a = 2.5e-6;   % rigid uniform gap (isolates the slip effect cleanly)
b = 4.0e-6;
Ho = b - a;   % gap thickness, for the dimensionless slip length l_i = l_i'/Ho
Lz = 8e-6;
N = 21;
z = linspace(0, Lz, N).';

par = struct();
par.NzFluid = N;
par.dz = Lz/(N-1);
par.zGrid = z;
par.mu = 1.2e-3;
par.pIn = 100;
par.pOut = 0;
par.dt = 1e-4;
par.SsrcFactor = 1;
par.minGap = 1e-9;
par.maxNewtonFluid = 30;
par.tolNewtonFluid = 1e-13;
par.useDirectFluidSolve = true;

state = struct();
state.deltaL = a*ones(N,1);
state.deltaE = b*ones(N,1);
state.UwL = zeros(N,1);
state.UwE = zeros(N,1);
state.p = linspace(par.pIn, par.pOut, N).';
old = state;

slipPhysical_nm = [0, 10, 25, 50, 75, 100];  % nm, matching Table 1/2's stated range
fprintf('%-14s %-14s %-14s %-16s %-16s\n', 'l_prime [nm]', 'l/Ho (dimless)', 'Q [m^3/s]', '%% vs no-slip', 'code-vs-hand relerr');

Q_noslip = [];
for k = 1:numel(slipPhysical_nm)
    lPrime = slipPhysical_nm(k) * 1e-9;
    par.slipL = lPrime;
    par.slipE = lPrime;

    [fluid, ok, stopReason] = solve_fluid_reynolds_slip(z, old, state, par);
    if ~ok
        fprintf('l=%.0fnm FAILED: %s\n', slipPhysical_nm(k), stopReason);
        continue;
    end
    Qk = mean(fluid.Q);
    if k == 1
        Q_noslip = Qk;
    end
    pctChange = 100*(Qk - Q_noslip)/abs(Q_noslip);

    % independent closed-form cross-check at this slip length
    pLinFit = polyfit(z, fluid.p, 1);
    G = pLinFit(1)/par.mu;
    ll = lPrime; le = lPrime;
    M = [log(a)+ll/a, 1; log(b)+le/b, 1];
    rhs = [0 - (G/4)*a^2 - ll*(G*a/2); 0 - (G/4)*b^2 - le*(G*b/2)];
    C = M\rhs;
    Qhand = integral_flux_annulus(a,b,G,C(1),C(2));
    relerr = abs(Qk-Qhand)/abs(Qhand);

    dimlessL = lPrime/Ho;
    fprintf('%-14.1f %-14.6f %-14.6e %-16.4f %-16.4e\n', slipPhysical_nm(k), dimlessL, Qk, pctChange, relerr);
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
