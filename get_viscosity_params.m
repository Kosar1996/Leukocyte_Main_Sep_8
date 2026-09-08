clc;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
par = out.par;

parLmismatch = par;
if isfield(par, 'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par, 'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par, 'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

etaBulkE = par.etaE * par.Ke / par.Ge;
etaBulkL = parLmismatch.etaE * parLmismatch.Ke / parLmismatch.Ge;

fprintf('Endothelium: Ge=%.4f Ke=%.4f eta=%.4f etaBulk=%.4f\n', par.Ge, par.Ke, par.etaE, etaBulkE);
fprintf('Leukocyte:   Ge=%.4f Ke=%.4f eta=%.4f etaBulk=%.4f\n', parLmismatch.Ge, parLmismatch.Ke, parLmismatch.etaE, etaBulkL);
