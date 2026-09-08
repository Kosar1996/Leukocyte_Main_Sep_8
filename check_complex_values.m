S = load('solid_endo_P300_wide.mat');
fprintf('meshE.nodes isreal: %d\n', isreal(S.meshE.nodes));
fprintf('uE_pre isreal: %d\n', isreal(S.uE_pre));
fprintf('any complex nodes: %d\n', nnz(imag(S.meshE.nodes(:))~=0));
fprintf('any complex uE_pre: %d\n', nnz(imag(S.uE_pre(:))~=0));
if ~isreal(S.meshE.nodes)
    idx = find(imag(S.meshE.nodes(:,1))~=0);
    disp(idx(1:min(5,numel(idx))));
    disp(S.meshE.nodes(idx(1:min(5,numel(idx))),:));
end
fprintf('SMOKE_TEST_STATUS: OK\n');
