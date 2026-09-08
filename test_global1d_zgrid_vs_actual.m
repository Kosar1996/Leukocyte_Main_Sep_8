clc;
file = 'out_2D_unified_domain_40step.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end
par = out.par;

fields_to_check = {'useGlobal1DPressure','useGlobal1DCoarseEdgeMesh','zMin','zMax','NzFluid', ...
    'global1DFineWindow','global1DCoarseDz','global1DFineDz'};
for i = 1:numel(fields_to_check)
    fn2 = fields_to_check{i};
    if isfield(par, fn2)
        v = par.(fn2);
        if isnumeric(v)
            fprintf('par.%s = %s\n', fn2, mat2str(v));
        elseif islogical(v)
            fprintf('par.%s = %d\n', fn2, v);
        end
    else
        fprintf('par.%s = <not present>\n', fn2);
    end
end

fprintf('\nCalling make_global_1d_z_grid(par) directly:\n');
z = make_global_1d_z_grid(par);
fprintf('numel(z) = %d\n', numel(z));
dZ = diff(z)*1e6;
fprintf('dz range: [%.4f, %.4f] um (uniform if equal)\n', min(dZ), max(dZ));
fprintf('First 10 z values (um): %s\n', mat2str(z(1:10)'*1e6, 4));
idx = find(z*1e6 > 3 & z*1e6 < 4);
fprintf('z values near 3.4-3.5um: %s\n', mat2str(z(idx)'*1e6, 6));

fprintf('\nCompare to actual meshF.Zp(1,:) in this file:\n');
k = out.stopStep;
fl = out.fluidHist{k};
actualZ = fl.meshF.Zp(1,:);
fprintf('numel = %d, dz = %.4f um (uniform)\n', numel(actualZ), (actualZ(2)-actualZ(1))*1e6);

fprintf('\nSMOKE_TEST_STATUS: OK\n');
