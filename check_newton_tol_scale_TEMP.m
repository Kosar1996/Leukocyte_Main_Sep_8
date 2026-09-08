clc;
S = load('out_debug_radialaxial.mat');
out = S.out;
k = out.stopStep;
st = out.stateHist{k};
par = out.par;
fl = out.fluidHist{k};

if ~isfield(fl,'tractionE') || ~isfield(fl,'tractionL')
    [fl.tractionL, fl.tractionE] = compute_bodyfitted_wall_traction(fl.meshF, fl, par);
end

meshE = out.meshE;
interfaceE = out.interfaceE;
baseE = out.baseE;

parCorrE = par;
if isfield(parCorrE, 'solidAbsTol'), parCorrE = rmfield(parCorrE, 'solidAbsTol'); end
parCorrE.debugVerbose = true;

fprintf('par.newtonTolSolid = %.3e\n', par.newtonTolSolid);
fprintf('par.newtonMaxItSolid = %d\n', par.newtonMaxItSolid);

ndof = size(meshE.nodes,1)*2;
u = st.uE;
[fixDofs, fixVals] = solid_support_conditions(baseE, par.supportE);
free = setdiff((1:ndof).', unique(fixDofs(:)));
u(fixDofs) = fixVals;

Fext = zeros(ndof,1);
[Fext, ~] = apply_interface_traction(meshE, u, Fext, interfaceE, fl.tractionE);
[Fint, ~] = assemble_finite_def_axisym(meshE, u, par);
[Fvisc, ~] = assemble_axisym_kelvin_voigt_viscous(meshE, u, st.uEPrev, par);
Fint = Fint + Fvisc;

R = Fint - Fext;
Rf = R(free);
resNorm = norm(Rf, inf);
refNorm = max([norm(Fext(free), inf), norm(Fint(free), inf), 1e-14]);
relNorm = resNorm / refNorm;

fprintf('\n||Fext(free)||_inf = %.6e N\n', norm(Fext(free), inf));
fprintf('||Fint(free)||_inf = %.6e N\n', norm(Fint(free), inf));
fprintf('resNorm (||Fint-Fext||_inf) = %.6e N\n', resNorm);
fprintf('refNorm = %.6e N\n', refNorm);
fprintf('relNorm = %.6e\n', relNorm);
fprintf('relNorm < newtonTolSolid? %d  (would trigger immediate return)\n', relNorm < par.newtonTolSolid);

fprintf('\nSMOKE_TEST_STATUS: OK\n');
