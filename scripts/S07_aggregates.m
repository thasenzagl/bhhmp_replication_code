clear all;
close all;

% Add functions
addpath(genpath('../functions'));

arnold_cutoff = true;
load_results = true;

%% 2. Load results

load('../results/matfiles/productivity_gains_results');
R = (1/beta)-(1-delta);

% Compute delta HHI
delta_HHI = zeros(1,size(out_none.s_ij,2));
HHI = zeros(1,size(out_none.s_ij,2));
for jj=1:size(out_none.s_ij,2)
    idx1 = logical(merge_ij(:,jj));
    delta_HHI(1,jj) = sum(out_none.s_ij(idx1, jj))^2 - sum(out_none.s_ij(idx1, jj).^2);
    HHI(1,jj) = sum(out_none.s_ij(idx1, jj))^2 + sum(out_none.s_ij(~idx1, jj).^2);
end


%% 3. SOLVE MODEL AFTER INCREASING MERGING FIRMS Z
if load_results == false

    % Define efficiency gains
    efficiency_gains = linspace(0., 0.1, 11);

    % Initialize an array or cell array to store outputs
    out_merges = cell(1, length(efficiency_gains));

    % Loop over efficiency gains
    for i = 1:length(efficiency_gains)
        % Copy original z_ij for each iteration to modify
        z_ij_modified = z_ij;

        % Apply the efficiency gain to the selected elements
        for jj = 1:J
            if Mj(jj) >= Mjcut
                z_ij_modified(logical(merge_ij(:,jj)), jj) = exp(efficiency_gains(i)) * z_ij(logical(merge_ij(:,jj)), jj);
            end
        end

        % Print update before solving the model
        fprintf('Solving model for efficiency gain of %d%%...\n', efficiency_gains(i) * 100);

        % Solve the model with the modified z_ij
        out_merges{i} = solve_model('select', z_ij_modified, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none);

        % Print update after solving the model
        fprintf('Model solved for efficiency gain of %d%%. Result stored in out_merges{%d}.\n', efficiency_gains(i) * 100, i);
    end
    
    %% Store results

    save('../results/matfiles/aggregates.mat', '-v7.3');

else
    
    load('../results/matfiles/aggregates.mat');
    
end

%% Arnold cutoffs

if arnold_cutoff == true
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
else
    jstar = zeros(1, J);
    jstar(Mj>=Mjcut) = 1;
end

%% Blocking Policies
HHI_grid = [[0.1, 0.01]; [0.18, 0.01]; [0.15, 0.01]; [0.25, 0.02]];

% Prepare a table or structure to hold the results
num_HHI_points = size(HHI_grid, 1);
tab = zeros(length(out_merges), num_HHI_points, 4);

% Create a structured array to store additional details
detailed_results = struct();

for i = 1:num_HHI_points
    HHI_cutoff = HHI_grid(i, 1);
    deltaHHI_cutoff = HHI_grid(i, 2);

    % Mergers that get blocked
    idx_blocked = logical((HHI >= HHI_cutoff) & (delta_HHI >= deltaHHI_cutoff));
    
    % Permitted mergers
    idx_permitted = ~idx_blocked;
    
    % Arnold cutoff
    idx_blocked(~jstar) = 0;
    idx_permitted(~jstar) = 0;
    
    assert(sum(idx_blocked + idx_permitted) == sum(jstar), 'Mismatch in permitted and blocked indexing');
    
    % Store results for out_none (No Merger Case)
    detailed_results(i).HHI = HHI(idx_permitted);
    detailed_results(i).deltaHHI = delta_HHI(idx_permitted);

    detailed_results(i).Y_none = sum(out_none.y_j(idx_permitted));
    detailed_results(i).Payroll_none = sum(out_none.w_j(idx_permitted) .* out_none.n_j(idx_permitted));
    detailed_results(i).K_none = sum(out_none.k_j(idx_permitted));
    detailed_results(i).Pi_none = sum(out_none.pi_j(idx_permitted));

    % Initialize storage for efficiency gains and corresponding Y, Payroll, K, Pi
    detailed_results(i).efficiency_gain = efficiency_gains;
    detailed_results(i).Y = zeros(length(out_merges), 1);
    detailed_results(i).Payroll = zeros(length(out_merges), 1);
    detailed_results(i).K = zeros(length(out_merges), 1);
    detailed_results(i).Pi = zeros(length(out_merges), 1);

    % Loop over all merge cases, including out_merge and out_merges
    for j = 1:length(out_merges)
        out_current = out_merges{j};
        
        % Compute economic variables for permitted mergers
        Y = sum(out_current.y_j(idx_permitted));
        Payroll = sum(out_current.w_j(idx_permitted) .* out_current.n_j(idx_permitted));
        K = sum(out_current.k_j(idx_permitted));
        Pi = sum(out_current.pi_j(idx_permitted));

        % Store market level results
        detailed_results(i).y_j(j, :) = out_current.y_j(idx_permitted);
        detailed_results(i).payroll_j(j, :) = out_current.w_j(idx_permitted) .* out_current.n_j(idx_permitted);
        detailed_results(i).k_j(j, :) = out_current.k_j(idx_permitted);
        detailed_results(i).pi_j(j, :) = out_current.pi_j(idx_permitted);
          
        % Store raw values for this efficiency gain level
        detailed_results(i).Y(j) = Y;
        detailed_results(i).Payroll(j) = Payroll;
        detailed_results(i).K(j) = K;
        detailed_results(i).Pi(j) = Pi;

        % Compute percentage changes and store them in `tab`
        tab(j, i, 1) = 100 * (log(Y) - log(detailed_results(i).Y_none));
        tab(j, i, 2) = 100 * ((Payroll / Y) - (detailed_results(i).Payroll_none / detailed_results(i).Y_none));
        tab(j, i, 3) = 100 * ((R * K / Y) - (R * detailed_results(i).K_none / detailed_results(i).Y_none));
        tab(j, i, 4) = 100 * ((Pi / Y) - (detailed_results(i).Pi_none / detailed_results(i).Y_none));
    end
end


%% Latex table
target_values = [0, 0.01, 0.02, 0.03, 0.04, 0.05];
indices = find(arrayfun(@(x) any(abs(efficiency_gains - x) < 1e-5), target_values));

names = {'0 Percent Efficiency Gain', '1 Percent Efficiency Gain', '2 Percent Efficiency Gain', ...
         '3 Percent Efficiency Gain', '4 Percent Efficiency Gain', '5 Percent Efficiency Gain'};
     
fid = fopen('../results/tables/table5_aggregates.tex','w');

fprintf(fid,'\\begin{tabular}{l @{\\hspace{2em}} c @{\\hspace{1.5em}}c}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,' & \\multicolumn{1}{l}{\\textbf{{------1982 guidelines------ $\\quad$ }}} & \\multicolumn{1}{l}{\\textbf{{------2010 guidelines------}}} \\\\ \n');
fprintf(fid,'DOJ/FTC Market Classification & \\textit{Highly Concentrated} & \\textit{Highly Concentrated} \\\\ \n');
fprintf(fid,'Threshold (HHI, $\\Delta$HHI) & (1800, 100) & (2500, 200) \\\\ \n');
fprintf(fid,' & (1) & (2) \\\\ \n');
fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{3}{l}{\\textbf{I. Percentage Change in Output Relative to No-Merger Baseline}} \\\\ \n');
for i = 1:6
    j = indices(i);
    fprintf(fid,'%s & %3.2f & %3.2f \\\\ \n', names{i}, tab(j, [2, 4], 1));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{3}{l}{\\textbf{II. Change in Labor Share Relative to No-Merger Baseline}} \\\\ \n');
for i = 1:6
    j = indices(i);
    fprintf(fid,'%s & %3.2f & %3.2f \\\\ \n', names{i}, tab(j, [2, 4], 2));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{3}{l}{\\textbf{III. Change in Profit Share Relative to No-Merger Baseline}} \\\\ \n');
for i = 1:6
    j = indices(i);
    fprintf(fid,'%s & %3.2f & %3.2f \\\\ \n', names{i}, tab(j, [2, 4], 4));    
end

% CLOSE TABLE
fprintf(fid,'\\bottomrule \n');
fprintf(fid,'\\end{tabular}');
fclose(fid);


%% Output and Payroll Growth Plot

% Create and configure figure
F1 = figure(1);
set(F1, 'Position', [0 0 1600 1000]); % Set figure size

% Set font and line width
font = 24;
title_font = 28; % Larger font size for title
lw = 5;

hold on;

% Plot with different line styles
plot(efficiency_gains * 100, 100 * (log(detailed_results(2).Y) - log(detailed_results(2).Y_none(1))), ...
    '-', 'LineWidth', lw, 'DisplayName', 'Output Growth'); % Solid Line

plot(efficiency_gains * 100, 100 * (log(detailed_results(2).Payroll) - log(detailed_results(2).Payroll_none(1))), ...
    '--', 'LineWidth', lw, 'DisplayName', 'Payroll Growth'); % Dashed Line

plot(efficiency_gains * 100, 100 * (log(detailed_results(2).Pi) - log(detailed_results(2).Pi_none(1))), ...
    ':', 'LineWidth', lw, 'DisplayName', 'Profit Growth'); % Dotted Line

% Labels with LaTeX interpreter
xlabel('Efficiency Gain (\%)', 'FontSize', font, 'Interpreter', 'latex');
ylabel('Percentage Change Relative to Pre-Merger Case (\%)', 'FontSize', font+4, 'Interpreter', 'latex');

% Title with (HHI, ?HHI) format and "Threshold" wording
title(['Impact of Mergers on Output, Payroll, and Profits' newline ...
       'Relative to Pre-Merger Case, with $(HHI, \Delta HHI)$ Thresholds of (1800, 100)'], ...
       'FontSize', title_font, 'Interpreter', 'latex');

% Configure legend
legend('Location', 'Best', 'FontSize', font, 'Interpreter', 'latex');

% Configure axes
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', font);
grid on;
hold off;

% Save figure as PNG and FIG file
saveas(gcf, '../results/figures/fig3_aggregates.png');



