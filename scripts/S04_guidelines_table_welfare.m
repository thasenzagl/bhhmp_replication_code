clear all;
close all;

% Add functions
addpath(genpath('../functions'));

arnold_cutoff = true;
load_results = true;

%% 2. Load results

load('../results/matfiles/productivity_gains_results.mat')

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

    z_ij1 = z_ij;
    z_ij2 = z_ij;
    z_ij3 = z_ij;
    z_ij4 = z_ij;    
    z_ij5 = z_ij;

    for jj=1:J
        if Mj(jj)>=Mjcut
            z_ij1(logical(merge_ij(:,jj)), jj) = exp(0.01)*z_ij(logical(merge_ij(:,jj)), jj);
            z_ij2(logical(merge_ij(:,jj)), jj) = exp(0.02)*z_ij(logical(merge_ij(:,jj)), jj);        
            z_ij3(logical(merge_ij(:,jj)), jj) = exp(0.03)*z_ij(logical(merge_ij(:,jj)), jj);
            z_ij4(logical(merge_ij(:,jj)), jj) = exp(0.04)*z_ij(logical(merge_ij(:,jj)), jj);            
            z_ij5(logical(merge_ij(:,jj)), jj) = exp(0.05)*z_ij(logical(merge_ij(:,jj)), jj);       
        end
    end

    out_merge1 = solve_model('select', z_ij1, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none);
    out_merge2 = solve_model('select', z_ij2, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none);
    out_merge3 = solve_model('select', z_ij3, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none);
    out_merge4 = solve_model('select', z_ij4, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none);    
    out_merge5 = solve_model('select', z_ij5, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none);
    
    %% Store results

    save('../results/matfiles/guidelines_results.mat', '-v7.3');

else
    
    load('../results/matfiles/guidelines_results.mat')
    
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

%% Compute welfare
[baseline_welfare1, new_welfare1] = compute_welfare(out_none, out_merge1, J, glob, param);
[baseline_welfare2, new_welfare2] = compute_welfare(out_none, out_merge2, J, glob, param);
[baseline_welfare3, new_welfare3] = compute_welfare(out_none, out_merge3, J, glob, param);
[baseline_welfare4, new_welfare4] = compute_welfare(out_none, out_merge4, J, glob, param);
[baseline_welfare5, new_welfare5] = compute_welfare(out_none, out_merge5, J, glob, param);

%% Make table

HHI_grid = [[0.1, 0.01]; [0.18, 0.01]; [0.15, 0.01]; [0.25, 0.02]; [Inf, Inf]];
tab = zeros(10, size(HHI_grid,1));

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
    
    % Average Efficiency gain of consumated mergers
    tab(1,i) = mean(out.Delta_j(idx_permitted,1));
    tab(3,i) = 100*(log(mean(new_welfare1(1,idx_permitted))) - log(mean(baseline_welfare1(1,idx_permitted))));
    tab(5,i) = 100*(log(mean(new_welfare2(1,idx_permitted))) - log(mean(baseline_welfare2(1,idx_permitted)))); 
    tab(7,i) = 100*(log(mean(new_welfare3(1,idx_permitted))) - log(mean(baseline_welfare3(1,idx_permitted))));
    tab(9,i) = 100*(log(mean(new_welfare4(1,idx_permitted))) - log(mean(baseline_welfare4(1,idx_permitted))));
    tab(11,i) = 100*(log(mean(new_welfare5(1,idx_permitted))) - log(mean(baseline_welfare5(1,idx_permitted))));
    
    % Average Efficiency gain of blocked mergers
    tab(2,i) = mean(out.Delta_j(idx_blocked,1));
    tab(4,i) = 100*(log(mean(new_welfare1(1,idx_blocked))) - log(mean(baseline_welfare1(1,idx_blocked))));
    tab(6,i) = 100*(log(mean(new_welfare2(1,idx_blocked))) - log(mean(baseline_welfare2(1,idx_blocked))));
    tab(8,i) = 100*(log(mean(new_welfare3(1,idx_blocked))) - log(mean(baseline_welfare3(1,idx_blocked))));
    tab(10,i) = 100*(log(mean(new_welfare4(1,idx_blocked))) - log(mean(baseline_welfare4(1,idx_blocked))));
    tab(12,i) = 100*(log(mean(new_welfare5(1,idx_blocked))) - log(mean(baseline_welfare5(1,idx_blocked))));
    
end

%% Latex table

names = {'Permitted mergers', 'Blocked mergers', 'Permitted mergers','Blocked mergers', 'Permitted mergers', 'Blocked mergers', 'Permitted mergers', 'Blocked mergers', 'Permitted mergers', 'Blocked mergers', 'Permitted mergers', 'Blocked mergers'};

fid = fopen('../results/tables/table4_guidelines_welfare.tex','w');

fprintf(fid,'\\begin{tabular}{l @{\\hspace{2em}} c @{\\hspace{1.5em}}c}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,' & \\multicolumn{1}{l}{\\textbf{{------1982 guidelines------ $\\quad$ }}} & \\multicolumn{1}{l}{\\textbf{{------2010 guidelines------}}} \\\\ \n');
fprintf(fid,'DOJ/FTC Market Classification & \\textit{Highly Concentrated} & \\textit{Highly Concentrated} \\\\ \n');
fprintf(fid,'Threshold (HHI, $\\Delta$HHI) & (1800, 100) & (2500, 200) \\\\ \n');
fprintf(fid,' & (1) & (2) \\\\ \n');
fprintf(fid,'\\midrule\n');

fprintf(fid, '\\multicolumn{3}{l}{\\textbf{I. Change in average welfare assuming 1 percent efficiency gain ($\\%%$)}} \\\\ \n');
for i = 1:2
    fprintf(fid, '%s & %3.1f & %3.1f \\\\ \n', names{i}, tab(i+2,2), tab(i+2,4));
end

fprintf(fid, '\\midrule\n');

fprintf(fid, '\\multicolumn{3}{l}{\\textbf{II. Change in average welfare assuming 2 percent efficiency gain ($\\%%$)}} \\\\ \n');
for i = 3:4
    fprintf(fid, '%s & %3.1f & %3.1f \\\\ \n', names{i}, tab(i+2,2), tab(i+2,4));
end

fprintf(fid, '\\midrule\n');

fprintf(fid, '\\multicolumn{3}{l}{\\textbf{III. Change in average welfare assuming 3 percent efficiency gain ($\\%%$)}} \\\\ \n');
for i = 5:6
    fprintf(fid, '%s & %3.1f & %3.1f \\\\ \n', names{i}, tab(i+2,2), tab(i+2,4));
end

fprintf(fid, '\\midrule\n');

fprintf(fid, '\\multicolumn{3}{l}{\\textbf{IV. Change in average welfare assuming 4 percent efficiency gain ($\\%%$)}} \\\\ \n');
for i = 7:8
    fprintf(fid, '%s & %3.1f & %3.1f \\\\ \n', names{i}, tab(i+2,2), tab(i+2,4));
end

fprintf(fid, '\\midrule\n');

fprintf(fid, '\\multicolumn{3}{l}{\\textbf{V. Change in average welfare assuming 5 percent efficiency gain ($\\%%$)}} \\\\ \n');
for i = 9:10
    fprintf(fid, '%s & %3.1f & %3.1f \\\\ \n', names{i}, tab(i+2,2), tab(i+2,4));
end

fprintf(fid, '\\bottomrule \n');
fprintf(fid, '\\end{tabular}');
fclose(fid);

