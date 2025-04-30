clear all;
close all;

% Add functions
addpath(genpath('../functions'));

arnold_cutoff = true;

%% 2. Load results

load('../results/matfiles/productivity_gains_results_alpha103');

% Compute delta HHI
delta_HHI = zeros(1,size(out_none.s_ij,2));
HHI = zeros(1,size(out_none.s_ij,2));
for jj=1:size(out_none.s_ij,2)
    idx1 = logical(merge_ij(:,jj));
    delta_HHI(1,jj) = sum(out_none.s_ij(idx1, jj))^2 - sum(out_none.s_ij(idx1, jj).^2);
    HHI(1,jj) = sum(out_none.s_ij(idx1, jj))^2 + sum(out_none.s_ij(~idx1, jj).^2);
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

%% Make table

HHI_grid = [[0.1, 0.01]; [0.18, 0.01]; [0.15, 0.01]; [0.25, 0.02]];
tab = zeros(2, size(HHI_grid,1));
Delta = [1, 2, 3, 4, 5];

for i=1:size(HHI_grid,1)
    
    HHI_cutoff = HHI_grid(i,1);
    deltaHHI_cutoff = HHI_grid(i,2);
    
    % Mergers that get blocked
    idx_blocked = logical(((HHI >= HHI_cutoff) & (delta_HHI >= deltaHHI_cutoff)));
    
    % Permitted mergers
    idx_permitted = ~idx_blocked;
    
    % Remove monopolies
    idx_blocked(~jstar) = 0;
    idx_permitted(~jstar) = 0;
    
    assert(sum(idx_blocked + idx_permitted) == sum(jstar));
    
    % AVERAGE WAGE OF ALL FIRMS
    % Probability of blocking a merger that would have generated WS gain
    tab(1,i) = mean(out.Delta_j(idx_permitted,1));
    
    % Probability of letting a merger through that generates WS loss
    tab(2,i) = mean(out.Delta_j(idx_blocked,1));
       
end


%% Latex table 1

names = {'Permitted mergers', 'Blocked mergers', 'Probability that a blocked merger yields WS gain', 'Probability that a permitted merger yields WS loss', 'Probability that a blocked merger yields WS gain', 'Probability that a permitted merger yields WS loss', 'Probability that a blocked merger yields WS gain', 'Probability that a permitted merger yields WS loss', 'Probability that a blocked merger yields WS gain', 'Probability that a permitted merger yields WS loss', 'Probability that a blocked merger yields WS gain', 'Probability that a permitted merger yields WS loss', 'Probability that a blocked merger yields WS gain', 'Probability that a permitted merger yields WS loss'};

fid = fopen('../results/tables/tableD4_alpha103.tex','w');

fprintf(fid,'\\begin{tabular}{l @{\\hspace{1em}} cc @{\\hspace{3em}} cc}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,' & \\multicolumn{2}{c}{\\textbf{A. 1982/2023 guidelines}} & \\multicolumn{2}{c}{\\textbf{B. 2010 guidelines}} \\\\ \n');
fprintf(fid,'\\cmidrule{2-3}\n');
fprintf(fid,'\\cmidrule{4-5}\n');

fprintf(fid,'DOJ/FTC market classification & Moderate & High & Moderate & High \\\\ \n');

fprintf(fid,'Threshold (HHI, $\\Delta$HHI) & (1000,  100) & (1800, 100) & (1500, 100) & (2500, 200) \\\\ \n');
fprintf(fid,' & (1) & (2) & (3) & (4) \\\\ \n');


fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{5}{l}{\\textbf{I. Average REG}} \\\\ \n');

for i = 1:2
    fprintf(fid,'%s & %3.2f & %3.2f & %3.2f & %3.2f \\\\ \n', names{i}, tab(i,:));    
end

% CLOSE TABLE
fprintf(fid,'\\bottomrule \n');
fprintf(fid,'\\end{tabular}');
fclose(fid);



















