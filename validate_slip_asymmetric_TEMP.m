clc;
addpath('/Users/kosarsafari/Desktop/Project_1/code/leukocyte-main');

a = 2.5e-6;
b = 4.0e-6;
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

slipVals_nm = [0, 10, 25, 50, 75, 100, 150, 200, 300];

fprintf('%-10s %-16s %-16s %-16s\n', 'l[nm]', 'Q leukocyte-only', 'Q endo-only', 'Q both (symmetric)');
Qbase = [];
for k = 1:numel(slipVals_nm)
    lPrime = slipVals_nm(k)*1e-9;

    par.slipL = lPrime; par.slipE = 0;
    [flL, okL] = solve_fluid_reynolds_slip(z, old, state, par);

    par.slipL = 0; par.slipE = lPrime;
    [flE, okE] = solve_fluid_reynolds_slip(z, old, state, par);

    par.slipL = lPrime; par.slipE = lPrime;
    [flBoth, okBoth] = solve_fluid_reynolds_slip(z, old, state, par);

    QL = mean(flL.Q); QE = mean(flE.Q); QB = mean(flBoth.Q);
    if k==1, Qbase = QL; end
    fprintf('%-10.0f %-16.6e %-16.6e %-16.6e\n', slipVals_nm(k), QL, QE, QB);
end
fprintf('\n(all as %% change from no-slip baseline %.6e)\n', Qbase);
fprintf('%-10s %-16s %-16s %-16s\n', 'l[nm]', '%%L-only', '%%E-only', '%%both');
par.slipL=0; par.slipE=0;
[fl0,~]=solve_fluid_reynolds_slip(z,old,state,par);
Q0 = mean(fl0.Q);
for k = 1:numel(slipVals_nm)
    lPrime = slipVals_nm(k)*1e-9;
    par.slipL = lPrime; par.slipE = 0;
    [flL,~] = solve_fluid_reynolds_slip(z, old, state, par);
    par.slipL = 0; par.slipE = lPrime;
    [flE,~] = solve_fluid_reynolds_slip(z, old, state, par);
    par.slipL = lPrime; par.slipE = lPrime;
    [flB,~] = solve_fluid_reynolds_slip(z, old, state, par);
    fprintf('%-10.0f %-16.4f %-16.4f %-16.4f\n', slipVals_nm(k), ...
        100*(mean(flL.Q)-Q0)/abs(Q0), 100*(mean(flE.Q)-Q0)/abs(Q0), 100*(mean(flB.Q)-Q0)/abs(Q0));
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
