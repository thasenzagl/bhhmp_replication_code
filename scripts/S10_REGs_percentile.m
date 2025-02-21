clear all
close all

%% LOAD RESULTS

load('../results/matfiles/productivity_gains_results');

% Add functions
addpath(genpath('../functions'))

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

out_none.s_ij = out_none.s_ij(:,jstar);
out.Delta_j = out.Delta_j(jstar);
merge_ij = merge_ij(:,jstar);

%% GENERATE DATA FOR FIGURE

HHI_bins = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 1];
delta_HHI_bins = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.05, 0.1, 0.5];

[hm20, delta_HHI20] = make_chart_hhi_delta(out, out_none, merge_ij, HHI_bins, delta_HHI_bins, 20);
[hm50, delta_HHI50] = make_chart_hhi_delta(out, out_none, merge_ij, HHI_bins, delta_HHI_bins, 50);

hm20 = flipud(hm20);
hm50 = flipud(hm50);

%% CUTOFF FIGURES

CustomXLabels = {'[0, 50)', '[50, 100)', '[100, 150)', '[150, 200)', '[200, 250)', '[250, 500)', '[500, 1000)', '[1000, 5000)'};
CustomYLabels = {'[0, 500)','[500, 1000)','[1000, 1500)', '[1500, 2000)', '[2000, 2500)', '[2500, 10000)'};
CustomYLabels = fliplr(CustomYLabels);

F1 = figure(1);
set(F1,'Pos',[0 0 1600 1000]);
font = 24;

color = generatecolormapthreshold([0 5 10 20 30 40],[0.9, 0.9447, 0.9741; 0.5, 0.7235, 0.8705; 0.3, 0.6129, 0.8187; 0.1, 0.5023, 0.7669; 0, 0.4470, 0.7410]);

subplot(1, 2, 1);
h = heatmap(hm20, 'Colormap',  color, 'ColorLimits', [0 40], 'FontSize', font);
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomYLabels;

h.FontSize = font;
h.MissingDataColor = 'w';
h.ColorbarVisible = 'off';
h.Title = 'A. 20th percentile of $\Delta^{*}$';
h.YLabel = 'Post-merger concentration ($HHI$)';
h.XLabel = 'Change in concentration ($\Delta HHI$)';
h.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).XLabel.Interpreter = 'latex';
h.NodeChildren(3).YLabel.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.CellLabelFormat = '%.2f';

subplot(1, 2, 2);
h = heatmap(hm50, 'Colormap',  color, 'ColorLimits', [0 40], 'FontSize', font);
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomYLabels;

h.XLabel = 'Small firm share';
h.YLabel = 'Large firm share';
h.FontSize = font;
h.MissingDataColor = 'w';
h.ColorbarVisible = 'off';
h.Title = 'B. Median value of $\Delta^{*}$';
h.YLabel = 'Post-merger concentration ($HHI$)';
h.XLabel = 'Change in concentration ($\Delta HHI$)';
h.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).XLabel.Interpreter = 'latex';
h.NodeChildren(3).YLabel.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.CellLabelFormat = '%.2f';

saveas(h,'../results/figures/figC1_REGs_percentiles_heatmap.png')

close all;
