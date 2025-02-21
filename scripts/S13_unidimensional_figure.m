clear all;
close all;

load_results = true;

%% SOLVE MODEL AFTER INCREASING MERGING FIRMS PRODUCTIVITY
if load_results == false
    
    load('../results/matfiles/productivity_gains_results.mat');
    HHI = out_none.HHIwn_j;

    delta_HHI = zeros(1,size(out_none.s_ij,2));    
    for jj=1:size(out_none.s_ij,2)
        idx = logical(merge_ij(:,jj));
        delta_HHI(1,jj) = sum(out_none.s_ij(idx, jj))^2 - sum(out_none.s_ij(idx, jj).^2);
    end

    HHI_cutoffs = [0.15, 0.2, 0.25, 0.3, 0.5, 1];
    delta_HHI_cutoffs = [0.005, 0.01, 0.02, 0.03, 0.1, 0.5];

    % Grid for productivity increases
    grid_size = 51;
    Delta_grid = linspace(0, 0.15, grid_size);

    wj_HHI = zeros(size(Delta_grid,2), size(HHI_cutoffs,2), J);
    wj_deltaHHI = zeros(size(Delta_grid,2), size(delta_HHI_cutoffs,2), J);

    for i=1:size(Delta_grid,2)

        z_ij1 = z_ij;
        for jj=1:J
            if Mj(jj)>=Mjcut
                z_ij1(logical(merge_ij(:,jj)), jj) = exp(Delta_grid(i))*z_ij(logical(merge_ij(:,jj)), jj);
            end
        end

        out_merge = solve_model('select', z_ij1, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none);

        for k = 1:size(HHI_cutoffs,2)
            idx = logical(HHI <= HHI_cutoffs(k));
            wj_HHI(i,k,idx) = out_merge.w_j(1,idx);
            wj_HHI(i,k,~idx) = out_none.w_j(1,~idx);
        end

        for k = 1:size(delta_HHI_cutoffs,2)
            idx = logical(delta_HHI <= delta_HHI_cutoffs(k));
            wj_deltaHHI(i,k,idx) = out_merge.w_j(1,idx);
            wj_deltaHHI(i,k,~idx) = out_none.w_j(1,~idx);
        end

    end

    wj_HHI_plot = zeros(grid_size, size(HHI_cutoffs,2));
    wj_deltaHHI_plot = zeros(grid_size, size(delta_HHI_cutoffs,2));

    mean_none = mean(out_none.w_j);

    for i = 1:size(HHI_cutoffs,2)
        for j = 1:size(Delta_grid,2)
            wj_HHI_plot(j,i) = 100*(log(mean(squeeze(wj_HHI(j,i,:)))) - log(mean_none));       
        end
    end

    for i = 1:size(delta_HHI_cutoffs,2)
        for j = 1:size(Delta_grid,2)
            wj_deltaHHI_plot(j,i) = 100*(log(mean(squeeze(wj_deltaHHI(j,i,:)))) - log(mean_none));             
        end
    end

    save('../results/matfiles/cutoff_rules', '-v7.3');

else
    
    load('../results/matfiles/cutoff_rules.mat');
    
end

%% FIND ARNOLD CUTOFF

n_ij_merge = reshape(out_none.n_ij(logical(merge_ij)), [2, sum(Mj>=Mjcut)]);
ntilde_j = (1/2)*sum(n_ij_merge);

ncut_max        = 200;
ncut            = [1:1:ncut_max];
for i=1:length(ncut)
    n_ij_temp       = n_ij_merge(:,ntilde_j>ncut(i));
    n_ij_temp       = n_ij_temp(:);
    median_size(i)  = median(n_ij_temp);
end

median_size_target  = 116;
[~,icut] = min(abs(median_size - median_size_target));

jstar = zeros(1, J);
jstar(Mj>=Mjcut) = (ntilde_j>ncut(icut));
jstar = logical(jstar);

%%
wj_HHI_plot = zeros(grid_size, size(HHI_cutoffs,2));
wj_deltaHHI_plot = zeros(grid_size, size(delta_HHI_cutoffs,2));

mean_none = mean(out_none.w_j(jstar));

for i = 1:size(HHI_cutoffs,2)
    for j = 1:size(Delta_grid,2)
        wj_HHI_plot(j,i) = 100*(log(mean(squeeze(wj_HHI(j,i,jstar)))) - log(mean_none));       
    end
end

for i = 1:size(delta_HHI_cutoffs,2)
    for j = 1:size(Delta_grid,2)
        wj_deltaHHI_plot(j,i) = 100*(log(mean(squeeze(wj_deltaHHI(j,i,jstar)))) - log(mean_none));             
    end
end

%% PLOT

F1 = figure(1);
set(F1,'Pos', [0 0 1600 700]);
font = 28;
lw = 3;
ms = 15;

subplot(1, 2, 1);
plot(100*Delta_grid, wj_HHI_plot(:,1),'-o', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid), 'LineWidth', lw, 'Color', [0, 48, 143]/255)
hold on;
plot(100*Delta_grid, wj_HHI_plot(:,2), '-+', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid), 'LineWidth', lw, 'Color', [6, 67, 165]/255)
hold on
plot(100*Delta_grid, wj_HHI_plot(:,3), '-v', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid),'LineWidth', lw, 'Color', [12, 86, 188]/255)
hold on;
plot(100*Delta_grid, wj_HHI_plot(:,4), '-square', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid),'LineWidth', lw, 'Color', [18, 106, 210]/255)
hold on;
plot(100*Delta_grid, wj_HHI_plot(:,5), '-diamond', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid),'LineWidth', lw, 'Color', [24, 125, 233]/255)
hold on;
plot(100*Delta_grid, wj_HHI_plot(:,6), '-^', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid),'LineWidth', lw, 'Color', [30, 144, 255]/255)
hold on;
xlabel('Efficiency gain ($\%$)', 'FontSize', font, 'interpreter', 'latex');
ylabel('Expected change in $\mathbf{W}_j$ due to merger ($\%$)', 'FontSize', font, 'interpreter', 'latex');
legend('Block if HHI$>$1500', 'Block if HHI$>$2000', 'Block if HHI$>$2500', 'Block if HHI$>$3000', 'Block if HHI$>$5000', 'Never block', 'FontSize', 24, 'interpreter', 'latex', 'Location', 'SouthEast')
title('A. HHI','FontSize', font, 'interpreter', 'latex');
set(gca,'ticklabelinterpreter','latex');
yline(0,'k--','LineWidth', lw, 'HandleVisibility', 'off')
xtickformat('%.d')
ax = gca;
ax.FontSize = font;
grid on;

subplot(1, 2, 2);
plot(100*Delta_grid, wj_deltaHHI_plot(:,1), '-o', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid), 'LineWidth', lw, 'Color', [0, 48, 143]/255)
hold on;
plot(100*Delta_grid, wj_deltaHHI_plot(:,2), '-+', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid), 'LineWidth', lw, 'Color', [6, 67, 165]/255)
hold on
plot(100*Delta_grid, wj_deltaHHI_plot(:,3), '-v', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid), 'LineWidth', lw, 'Color', [12, 86, 188]/255)
hold on;
plot(100*Delta_grid, wj_deltaHHI_plot(:,4), '-square', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid), 'LineWidth', lw, 'Color', [18, 106, 210]/255)
hold on;
plot(100*Delta_grid, wj_deltaHHI_plot(:,5), '-diamond', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid), 'LineWidth', lw, 'Color', [24, 125, 233]/255)
hold on;
plot(100*Delta_grid, wj_deltaHHI_plot(:,6), '-^', 'MarkerSize', ms, 'MarkerIndices', 1:4:length(Delta_grid), 'LineWidth', lw, 'Color', [30, 144, 255]/255)
hold on;
xlabel('Efficiency gain ($\%$)', 'FontSize', font, 'interpreter', 'latex');
ylabel('Expected change in $\mathbf{W}_j$ due to merger ($\%$)', 'FontSize', font, 'interpreter', 'latex');
legend('Block if $\Delta$HHI$>$50', 'Block if $\Delta$HHI$>$100', 'Block if $\Delta$HHI$>$200', 'Block if $\Delta$HHI$>$300', 'Block if $\Delta$HHI$>$1000', 'Never block', 'FontSize', 24, 'interpreter', 'latex', 'Location', 'SouthEast')
title('B. $\Delta$HHI','FontSize', font, 'interpreter', 'latex');
set(gca,'ticklabelinterpreter','latex');
yline(0,'k--','LineWidth', lw, 'HandleVisibility', 'off')
xtickformat('%.d')
ax = gca;
ax.FontSize = font;
grid on;ytickformat('%2.1f');

saveas(gcf, '../results/figures/figd1_unidimensional.png')

close all;


