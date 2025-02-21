clear all
close all

%% LOAD RESULTS

res = load('../results/matfiles/productivity_gains_results_symmetric.mat');

% Add functions
addpath(genpath('../functions'));

%% Compute HHI

res.HHI_pre = res.out_none.HHIwn_j;

% Compute post-merger HHI    
res.out_none.swn_ij = bsxfun(@rdivide,res.out_none.w_ij.*res.out_none.n_ij,sum(res.out_none.w_ij.*res.out_none.n_ij));

% Compute delta HHI
res.delta_HHI = zeros(1,size(res.out_none.swn_ij,2));    
for jj=1:size(res.out_none.swn_ij,2)
    idx = logical(res.merge_ij(:,jj));
    res.delta_HHI(1,jj) = sum(res.out_none.swn_ij(idx, jj))^2 - sum(res.out_none.swn_ij(idx, jj).^2);
end

%% FIGURE 4: AVERAGE WAGES

res.out_none.ave_wages = 0.5*sum(reshape(res.out_none.w_ij(logical(res.merge_ij)), [2, sum(res.Mj>=res.Mjcut)]));
res.out_merge.ave_wages = 0.5*sum(reshape(res.out_merge.w_ij(logical(res.merge_ij)), [2, sum(res.Mj>=res.Mjcut)]));
res.pct_chg_ave_wages = 100*(log(res.out_merge.ave_wages) - log(res.out_none.ave_wages));

%% FIGURE 5: TOTAL EMPLOYMENT

res.out_none.total_emp = sum(reshape(res.out_none.n_ij(logical(res.merge_ij)), [2, sum(res.Mj>=res.Mjcut)]));
res.out_merge.total_emp = sum(reshape(res.out_merge.n_ij(logical(res.merge_ij)), [2, sum(res.Mj>=res.Mjcut)]));
res.pct_chg_total_emp = 100*(log(res.out_merge.total_emp) - log(res.out_none.total_emp));

%% FIGURE 6: DIVERSION RATIO

% Compute pre-merger shares
res.out_none.swn_ij = bsxfun(@rdivide,res.out_none.w_ij.*res.out_none.n_ij,sum(res.out_none.w_ij.*res.out_none.n_ij));

% Compute diversion rations
res.out_none.D = zeros(1,res.J);

for jj=1:size(res.out_none.swn_ij,2) 
    res.out_none.D(1,jj) = compute_diversion(res.out_none.swn_ij(1,jj), res.out_none.swn_ij(1,jj), res.out_none.swn_ij(2,jj), res.param);
end

%% COMBINE FIGURES

F1 = figure(1);
set(F1,'Pos',[0 0 1600 1000]);
font = 24;
lw = 3;

subplot(2, 3, 1);
plot(res.Mj, 10000*res.HHI_pre, 'LineWidth', lw, 'color', rgb('DarkBlue'));
ylabel('$HHI$', 'FontSize', font, 'interpreter', 'latex');
xlabel('Number of firms: $M$', 'FontSize', font, 'interpreter', 'latex');
title('A. Pre-merger concentration', 'FontSize', font, 'interpreter', 'latex');
set(gca,'ticklabelinterpreter','latex');
ax = gca;
ax.FontSize = font;
grid on;

subplot(2, 3, 2);
plot(res.Mj, 10000*res.delta_HHI, 'LineWidth', lw, 'color', rgb('DarkBlue'));
ylabel('Change in $HHI$', 'FontSize', font, 'interpreter', 'latex');
xlabel('Number of firms: $M$', 'FontSize', font, 'interpreter', 'latex');
title('B. Concentration increase','FontSize', font, 'interpreter', 'latex');
set(gca,'ticklabelinterpreter','latex');
ax = gca;
ax.FontSize = font;
grid on;

subplot(2, 3, 3);
plot(res.Mj, res.pct_chg_total_emp, 'LineWidth', lw, 'color', rgb('DarkBlue'));
ylabel('Change in total employment ($\%$)', 'FontSize', font, 'interpreter', 'latex');
xlabel('Number of firms: $M$', 'FontSize', font, 'interpreter', 'latex');
title('C. Employment at merging firms','FontSize', font, 'interpreter', 'latex');
set(gca,'ticklabelinterpreter','latex');
ax = gca;
ax.FontSize = font;
grid on;

subplot(2, 3, 4);
plot(res.Mj, res.pct_chg_ave_wages, 'LineWidth', lw, 'color', rgb('DarkBlue'));
ylabel('Change in average wage ($\%$)', 'FontSize', font, 'interpreter', 'latex');
xlabel('Number of firms: $M$', 'FontSize', font, 'interpreter', 'latex');
title('D. Wages of merging firms','FontSize', font, 'interpreter', 'latex');
set(gca,'ticklabelinterpreter','latex');
ax = gca;
ax.FontSize = font;
grid on;

subplot(2, 3, 5);
plot(res.Mj, res.out_none.D, 'LineWidth', lw, 'color', rgb('DarkBlue'));
ylabel('Downward wage pressure', 'FontSize', font, 'interpreter', 'latex');
xlabel('Number of firms: $M$', 'FontSize', font, 'interpreter', 'latex');
title('E. Pre-merger downward wage pressure','FontSize', font, 'interpreter', 'latex');
set(gca,'ticklabelinterpreter','latex');
ax = gca;
ax.FontSize = font;
grid on;

subplot(2, 3, 6);
plot(res.Mj, res.out.Delta_j, 'LineWidth', lw, 'color', rgb('DarkBlue'));
ylabel('WS-neutral $\Delta^{*}$ ($\%$)', 'FontSize', font, 'interpreter', 'latex');
xlabel('Number of firms: $M$', 'FontSize', font, 'interpreter', 'latex');
title('F. WS-neutral productivity gain','FontSize', font, 'interpreter', 'latex');
set(gca,'ticklabelinterpreter','latex');
ax = gca;
ax.FontSize = font;
grid on;

saveas(gcf, '../results/figures/figb1_symmetric_mergers.png')

close all;



