S = load('solid_endo_P300.mat');
disp(fieldnames(S));
disp('par fields:');
disp(fieldnames(S.par));
fprintf('par.zMin=%g par.zMax=%g\n', S.par.zMin, S.par.zMax);
if isfield(S.par,'NzSolid'); fprintf('par.NzSolid=%g\n', S.par.NzSolid); end
if isfield(S,'meshE')
    fprintf('meshE nodes z range: [%g, %g], nNodes=%d\n', min(S.meshE.nodes(:,2)), max(S.meshE.nodes(:,2)), size(S.meshE.nodes,1));
else
    fprintf('No meshE field in this file.\n');
end
if isfield(S,'uE_pre'); fprintf('uE_pre length=%d\n', numel(S.uE_pre)); end
fprintf('SMOKE_TEST_STATUS: OK\n');
