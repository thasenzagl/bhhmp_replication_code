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

bins = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
table20 = make_chart_share_delta(bins, merge_ij, Delta_j, Mj, Mjcut, J, out_none, 20);
table50 = make_chart_share_delta(bins, merge_ij, Delta_j, Mj, Mjcut, J, out_none, 50);

hm20 = table20(:,1:end-5);
hm50 = table50(:,1:end-5);

hm20 = flipud(hm20);
hm50 = flipud(hm50);

%% CUTOFF FIGURES

CustomXLabels = {'[0, 0.05)', '[0.05, 0.1)', '[0.1, 0.2)', '[0.2, 0.3)', '[0.3, 0.4)', '[0.4, 0.5)'};
CustomYLabels = {'[0, 0.05)', '[0.05, 0.1)', '[0.1, 0.2)', '[0.2, 0.3)', '[0.3, 0.4)', '[0.4, 0.5)', '[0.5, 0.6)', '[0.6, 0.7)', '[0.7, 0.8)', '[0.8, 0.9)', '[0.9, 1.0)'};
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
h.Title = 'A: 20th percentile of $\Delta^{*}$';
h.XLabel = 'Initial small firm payroll share';
h.YLabel = 'Initial large firm payroll share';
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
h.XLabel = 'Initial small firm payroll share';
h.YLabel = 'Initial large firm payroll share';
h.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).XLabel.Interpreter = 'latex';
h.NodeChildren(3).YLabel.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.CellLabelFormat = '%.2f';

saveas(h,'../results/figures/figC2_REGs_shares_heatmap.png')

close all;
