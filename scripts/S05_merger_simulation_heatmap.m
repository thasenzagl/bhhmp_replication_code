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
Delta_j = out.Delta_j(jstar);
merge_ij = merge_ij(:,jstar);
Mj = Mj(jstar);

%% GENERATE DATA FOR FIGURE

bins_small = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
bins_large = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
hm_share = make_chart_share_fraction(bins_small, bins_large, merge_ij, Delta_j, Mj, out_none, 5);
hm_share = 100.*hm_share(:,1:end-5);
hm_share = flipud(hm_share);

HHI_bins = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 1];
delta_HHI_bins = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.05, 0.1, 0.5];

hm_hhi = make_chart_hhi_fraction(HHI_bins, delta_HHI_bins, merge_ij, Delta_j, out_none, 5);
hm_hhi = 100*hm_hhi;
hm_hhi = flipud(hm_hhi);

%% FIGURES

CustomXLabels_hhi = {'[0, 50)', '[50, 100)', '[100, 150)', '[150, 200)', '[200, 250)', '[250, 500)', '[500, 1000)', '[1000, 5000)'};
CustomYLabels_hhi = {'[0, 500)','[500, 1000)','[1000, 1500)', '[1500, 2000)', '[2000, 2500)', '[2500, 10000)'};
CustomYLabels_hhi = fliplr(CustomYLabels_hhi);

F1 = figure(1);
set(F1,'Pos',[0 0 1600 1000]);
font = 24;

color = generatecolormapthreshold([0 20 40 60 80 100],[0.6, 0, 0; 0.7, 0.1, 0.1; 0.8, 0.3, 0.3; 0.9, 0.5, 0.5; 1, 0.7, 0.7]);

subplot(1, 2, 1);
h = heatmap(hm_hhi, 'Colormap',  color, 'ColorLimits', [0 100], 'FontSize', font);
h.XDisplayLabels = CustomXLabels_hhi;
h.YDisplayLabels = CustomYLabels_hhi;

h.FontSize = font;
h.MissingDataColor = 'w';
h.ColorbarVisible = 'off';
h.Title = 'A. 5 percent efficiency increase';
h.XLabel = 'Change in concentration ($\Delta HHI$)';
h.YLabel = 'Post-merger concentration ($HHI$)';
h.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).XLabel.Interpreter = 'latex';
h.NodeChildren(3).YLabel.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.CellLabelFormat = '%.1f';


CustomXLabels_share = {'[0, 0.05)', '[0.05, 0.1)', '[0.1, 0.2)', '[0.2, 0.3)', '[0.3, 0.4)', '[0.4, 0.5)'};
CustomYLabels_share = {'[0, 0.05)', '[0.05, 0.1)', '[0.1, 0.2)', '[0.2, 0.3)', '[0.3, 0.4)', '[0.4, 0.5)', '[0.5, 0.6)', '[0.6, 0.7)', '[0.7, 0.8)', '[0.8, 0.9)', '[0.9, 1.0)'};
CustomYLabels_share = fliplr(CustomYLabels_share);

subplot(1, 2, 2);
h = heatmap(hm_share, 'Colormap',  color, 'ColorLimits', [0 100], 'FontSize', font);
h.XDisplayLabels = CustomXLabels_share;
h.YDisplayLabels = CustomYLabels_share;

h.FontSize = font;
h.MissingDataColor = 'w';
h.ColorbarVisible = 'off';
h.Title = 'B. 5 percent efficiency increase';
h.XLabel = 'Initial small firm payroll share';
h.YLabel = 'Initial large firm payroll share';
h.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).XLabel.Interpreter = 'latex';
h.NodeChildren(3).YLabel.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.CellLabelFormat = '%.1f';

saveas(h,'../results/figures/fig1_merger_simulation_heatmap.png')


close all;
