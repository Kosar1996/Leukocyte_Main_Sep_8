S = load('out_2D_unified_domain_40step.mat');
fn = fieldnames(S); out = S.(fn{1});
fl = out.fluidHist{out.stopStep};
mraw = fl.meshF;
fprintf('pCell size: %s, any NaN: %d, any Inf: %d\n', mat2str(size(fl.pCell)), any(isnan(fl.pCell)), any(isinf(fl.pCell)));
fprintf('pCell range: [%.4f, %.4f]\n', min(fl.pCell), max(fl.pCell));
fprintf('Nr=%d Nz=%d\n', mraw.Nr, mraw.Nz);
Zcol = mraw.Zp(1,:);
fprintf('Zp(1,:) range: [%.4f, %.4f] um (should be full z domain)\n', min(Zcol)*1e6, max(Zcol)*1e6);
for jtest = round(linspace(1, mraw.Nz, 8))
    fprintf('j=%d z=%.4f um: Rp column = %s\n', jtest, mraw.Zp(1,jtest)*1e6, mat2str(mraw.Rp(:,jtest)'*1e6, 4));
end
pgrid = reshape(fl.pCell, mraw.Nr, mraw.Nz);
fprintf('pgrid col 50 (z=%.3f um): %s\n', mraw.Zp(1,50)*1e6, mat2str(pgrid(:,50)', 4));
