clc;
S = load('out_pure2DMAC_fixedpoint.mat');
out = S.out;
k = out.stopStep;
st = out.stateHist{k};
par = out.par;
dtStep = out.dtHist(k);
fl = out.fluidHist{k};

parLmismatch = par;
if isfield(par,'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par,'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par,'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

mismatchOpts = struct('nQuery', 31, 'epsFrac', 0.1, 'trimFrac', 0.05);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, dtStep, fl, mismatchOpts);

normFromTractionE = interp1(fl.tractionE.z(:), fl.tractionE.normal(:), endo.z(:), 'linear', 'extrap');

% Independent estimate of the viscous normal-stress term (2*mu*durdr) at
% the endothelium wall, mirroring estimate_wall_shear_bodyfitted.m's
% local-polyfit approach but for radial velocity's radial gradient
% instead of axial velocity's radial gradient. No-penetration wall BC:
% approximate u_r at the wall as 0 (rigid/near-rigid interface over one
% step; the correction loop's own displacement updates are ~1e-7 m, tiny
% vs the flow's radial velocity scale).
mesh = fl.meshF;
Nz = mesh.Nz; Nr = mesh.Nr;
nFit = min(4, Nr);
viscNormalE = nan(Nz,1);
for j = 1:Nz
    rE = st.deltaE(j);
    idxE = (Nr-nFit+1):Nr;
    rFitE = [mesh.Rp(idxE,j); rE];
    uFitE = [fl.urC(idxE,j); 0];
    cE = polyfit(rFitE, uFitE, 1);
    viscNormalE(j) = 2*par.mu*cE(1);
end
mesh_zc = mesh.zc(:);
viscNormalE_atEndoZ = interp1(mesh_zc, viscNormalE, endo.z(:), 'linear', 'extrap');

fprintf('%-8s %-14s %-14s %-12s %-14s %-14s\n', ...
    'z[um]', 'tnFluid(meas)', 'normal(loaded)', 'gap', '2*mu*durdr', 'gap-2*mu*durdr');
for i = 1:5:numel(endo.z)
    gapVal = endo.tnFluid(i) - normFromTractionE(i);
    fprintf('%-8.3f %-14.4f %-14.4f %-12.4f %-14.4f %-14.4f\n', ...
        endo.z(i)*1e6, endo.tnFluid(i), normFromTractionE(i), gapVal, ...
        viscNormalE_atEndoZ(i), gapVal - viscNormalE_atEndoZ(i));
end

fprintf('\nSMOKE_TEST_STATUS: OK\n');
