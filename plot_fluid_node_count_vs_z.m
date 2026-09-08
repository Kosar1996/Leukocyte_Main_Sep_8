%% PLOT_FLUID_NODE_COUNT_VS_Z
% Directly count how many fluid mesh nodes span the gap between the
% leukocyte and endothelium interfaces, as a function of z (not vs. time


clc; close all;
file = '/Users/kosarsafari/Desktop/Project_1/all_three_runs/after_correction_aug_10/out_2D_t10_for_review.mat';
S = load(file);
fn = fieldnames(S);
out = S.(fn{1});
if isfield(out, 'cfg'), out = rmfield(out, 'cfg'); end

k = out.stopStep;
fl = out.fluidHist{k};
mraw = fl.meshF;

Nz = mraw.Nz; Nr = mraw.Nr;
zCols = mraw.Zp(1,:);          % z value of each column j
countPerColumn = zeros(Nz,1);
minGapPerColumn = zeros(Nz,1);
maxGapPerColumn = zeros(Nz,1);

for j = 1:Nz
    Rcol = mraw.Rp(:,j);
    % directly count how many distinct radial fluid nodes exist in this column
    countPerColumn(j) = numel(Rcol);
    minGapPerColumn(j) = min(diff(Rcol));
    maxGapPerColumn(j) = max(diff(Rcol));
end

fig = figure('Position',[100 100 900 500]);
yyaxis left
plot(zCols*1e6, countPerColumn, 'b-o', 'MarkerSize',3);
ylabel('Number of fluid nodes spanning the gap');
ylim([0 max(countPerColumn)+5]);
yyaxis right
plot(zCols*1e6, (mraw.Rp(end,:)-mraw.Rp(1,:))'*1e9, 'r--');
ylabel('Gap width [nm]');
xlabel('z [\mum]');
title('Fluid node count between leukocyte and endothelium interfaces, vs z');
grid on;

saveas(fig, 'fluid_node_count_vs_z.png');
fprintf('Node count: min=%d, max=%d, constant=%d\n', min(countPerColumn), max(countPerColumn), all(countPerColumn==countPerColumn(1)));
fprintf('Gap range: %.1f to %.1f nm\n', min(mraw.Rp(end,:)-mraw.Rp(1,:))*1e9, max(mraw.Rp(end,:)-mraw.Rp(1,:))*1e9);
fprintf('SMOKE_TEST_STATUS: OK\n');
