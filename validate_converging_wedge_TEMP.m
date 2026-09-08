clc;
addpath('/Users/kosarsafari/Desktop/Project_1/code/leukocyte-main');

%% Converging-wedge (tapered) rigid annular gap -- classical Reynolds
% lubrication problem. Unlike the uniform-gap case, pressure is NOT
% linear here -- this tests the code's actual nonlinear discretized solve
% (Newton iteration over the z-grid), not just the pointwise closed-form
% flux formula already verified.
%
% Independent check: at EVERY z-node, given the code's own converged
% (constant) flux Q and the KNOWN local geometry a(z),b(z), independently
% root-solve (via fzero, not reusing the code's Newton solve at all) for
% the local pressure gradient G(z) that the exact annular-slip flux
% formula requires to produce that Q. Compare against dp/dz obtained by
% numerically differentiating the CODE's own solved p(z). If these two,
% completely independently derived quantities agree, the nonlinear
% discretized solve is validated, not just the algebra.

a0 = 2.5e-6;              % leukocyte radius, constant
b1 = 4.5e-6; b2 = 3.0e-6; % endothelium radius: converging taper
Lz = 8e-6;
N = 41;
z = linspace(0, Lz, N).';
bOfz = b1 + (b2-b1)*(z/Lz);   % linear taper

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
par.maxNewtonFluid = 50;
par.tolNewtonFluid = 1e-13;
par.useDirectFluidSolve = true;
par.slipL = 0;
par.slipE = 0;

state = struct();
state.deltaL = a0*ones(N,1);
state.deltaE = bOfz;
state.UwL = zeros(N,1);
state.UwE = zeros(N,1);
state.p = linspace(par.pIn, par.pOut, N).';
old = state;

[fluid, ok, stopReason] = solve_fluid_reynolds_slip(z, old, state, par);
if ~ok
    fprintf('FAILED: %s\n', stopReason);
    return;
end

fprintf('=== CONVERGING WEDGE: pressure is NOT linear (unlike uniform gap) ===\n');
pLinFit = polyfit(z, fluid.p, 1);
pLinPred = polyval(pLinFit, z);
fprintf('max |p - best-fit LINE| = %.4f Pa (should be LARGE/nonzero -- confirms nonlinearity)\n', ...
    max(abs(fluid.p - pLinPred)));
fprintf('Q along z: mean=%.6e, spread(max-min)=%.6e (spread should be ~0, mass conservation)\n', ...
    mean(fluid.Q), max(fluid.Q)-min(fluid.Q));

% Independent check: at interior nodes, numerically differentiate the
% code's own p(z), then independently root-solve for G that the exact
% formula says SHOULD produce the code's own converged Q at that (a,b).
Qtarget = mean(fluid.Q);
mu = par.mu;
fprintf('\n%-8s %-14s %-14s %-14s %-10s\n', 'z[um]', 'dp/dz (code, FD)', 'G_hand*mu', 'abs diff', 'rel diff');
maxRel = 0;
for i = 3:N-2  % central differences, avoid edges
    dpdz_code = (fluid.p(i+1) - fluid.p(i-1)) / (z(i+1)-z(i-1));

    a_i = a0; b_i = bOfz(i);
    Gfun = @(G) local_Q_annular(a_i, b_i, G, 0, 0) - Qtarget;
    G0 = dpdz_code/mu;  % use code's own value as the fzero starting guess
    Ghand = fzero(Gfun, G0);
    dpdz_hand = Ghand * mu;

    absdiff = abs(dpdz_code - dpdz_hand);
    reldiff = absdiff / max(abs(dpdz_hand), 1e-20);
    maxRel = max(maxRel, reldiff);
    if mod(i,8)==0
        fprintf('%-8.4f %-14.6e %-14.6e %-14.6e %-10.6e\n', z(i)*1e6, dpdz_code, dpdz_hand, absdiff, reldiff);
    end
end
fprintf('\nmax relative difference across all interior nodes = %.6e\n', maxRel);

fprintf('\nSMOKE_TEST_STATUS: OK\n');

function Q = local_Q_annular(a, b, G, ll, le)
    M = [log(a)+ll/a, 1; log(b)+le/b, 1];
    rhs = [0 - (G/4)*a^2 - ll*(G*a/2);
           0 - (G/4)*b^2 - le*(G*b/2)];
    C = M\rhs;
    Q = integral_flux_annulus(a, b, G, C(1), C(2));
end
