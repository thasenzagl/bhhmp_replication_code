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

%% Make table

HHI_grid = [[0.1, 0.01]; [0.18, 0.01]; [0.15, 0.01]; [0.25, 0.02]];
tab = zeros(10, size(HHI_grid,1));

for i=1:size(HHI_grid,1)
    
    HHI_cutoff = HHI_grid(i,1);
    deltaHHI_cutoff = HHI_grid(i,2);
    
    % Mergers that get blocked
    idx_blocked = logical(((HHI >= HHI_cutoff) & (delta_HHI >= deltaHHI_cutoff)));
    
    % Permitted mergers
    idx_permitted = ~idx_blocked;
    
    % Arnold cutoff
    idx_blocked(~jstar) = 0;
    idx_permitted(~jstar) = 0;
    
    assert(sum(idx_blocked + idx_permitted) == sum(jstar));
    
    % AVERAGE WAGE OF ALL FIRMS
    % Average Efficiency gain of consumated mergers
    tab(1,i) = mean(out.Delta_j(idx_permitted,1));
    tab(3,i) = 100*(log(mean(out_merge1.w_j(1,idx_permitted))) - log(mean(out_none.w_j(1,idx_permitted))));
    tab(5,i) = 100*(log(mean(out_merge2.w_j(1,idx_permitted))) - log(mean(out_none.w_j(1,idx_permitted)))); 
    tab(7,i) = 100*(log(mean(out_merge3.w_j(1,idx_permitted))) - log(mean(out_none.w_j(1,idx_permitted)))); 
    tab(9,i) = 100*(log(mean(out_merge4.w_j(1,idx_permitted))) - log(mean(out_none.w_j(1,idx_permitted)))); 
    tab(11,i) = 100*(log(mean(out_merge5.w_j(1,idx_permitted))) - log(mean(out_none.w_j(1,idx_permitted)))); 
    
    % Average Efficiency gain of blocked mergers
    tab(2,i) = mean(out.Delta_j(idx_blocked,1));
    tab(4,i) = 100*(log(mean(out_merge1.w_j(1,idx_blocked))) - log(mean(out_none.w_j(1,idx_blocked))));
    tab(6,i) = 100*(log(mean(out_merge2.w_j(1,idx_blocked))) - log(mean(out_none.w_j(1,idx_blocked))));
    tab(8,i) = 100*(log(mean(out_merge3.w_j(1,idx_blocked))) - log(mean(out_none.w_j(1,idx_blocked))));
    tab(10,i) = 100*(log(mean(out_merge4.w_j(1,idx_blocked))) - log(mean(out_none.w_j(1,idx_blocked))));
    tab(12,i) = 100*(log(mean(out_merge5.w_j(1,idx_blocked))) - log(mean(out_none.w_j(1,idx_blocked))));
    
end


%% Make Table 3

names = {'Permitted mergers', 'Blocked mergers', 'Permitted mergers','Blocked mergers', 'Permitted mergers', 'Blocked mergers', 'Permitted mergers', 'Blocked mergers', 'Permitted mergers', 'Blocked mergers', 'Permitted mergers', 'Blocked mergers'};

fid = fopen('../results/tables/table3_guidelines_wages.tex','w');

fprintf(fid,'\\begin{tabular}{l @{\\hspace{2em}} c @{\\hspace{1.5em}}c}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,' & \\multicolumn{1}{l}{\\textbf{{------1982 guidelines------ $\\quad$ }}} & \\multicolumn{1}{l}{\\textbf{{------2010 guidelines------}}} \\\\ \n');
fprintf(fid,'DOJ/FTC Market Classification & \\textit{Highly Concentrated} & \\textit{Highly Concentrated} \\\\ \n');
fprintf(fid,'Threshold (HHI, $\\Delta$HHI) & (1800, 100) & (2500, 200) \\\\ \n');
fprintf(fid,' & (1) & (2) \\\\ \n');
fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{3}{l}{\\textbf{I. Average REG}} \\\\ \n');

for i = 1:2
    fprintf(fid,'%s & %3.2f & %3.2f  \\\\ \n', names{i}, tab(i,[2,4]));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{3}{l}{\\textbf{II. Change in average $\\mathbf{W}_j$ assuming 1 percent efficiency gain ($\\%%$)}} \\\\ \n');

for i = 3:4
    fprintf(fid,'%s & %3.2f & %3.2f \\\\ \n', names{i}, tab(i,[2,4]));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{3}{l}{\\textbf{III. Change in average $\\mathbf{W}_j$ assuming 2 percent efficiency gain ($\\%%$)}} \\\\ \n');

for i = 5:6
    fprintf(fid,'%s & %3.2f & %3.2f \\\\ \n', names{i}, tab(i,[2,4]));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{3}{l}{\\textbf{IV. Change in average $\\mathbf{W}_j$ assuming 3 percent efficiency gain ($\\%%$)}} \\\\ \n');

for i = 7:8
    fprintf(fid,'%s & %3.2f & %3.2f \\\\ \n', names{i}, tab(i,[2,4]));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{3}{l}{\\textbf{V. Change in average $\\mathbf{W}_j$ assuming 4 percent efficiency gain ($\\%%$)}} \\\\ \n');

for i = 9:10
    fprintf(fid,'%s & %3.2f & %3.2f \\\\ \n', names{i}, tab(i,[2,4]));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{3}{l}{\\textbf{VI. Change in average $\\mathbf{W}_j$ assuming 5 percent efficiency gain ($\\%%$)}} \\\\ \n');

for i = 11:12
    fprintf(fid,'%s & %3.2f & %3.2f \\\\ \n', names{i}, tab(i,[2,4]));    
end

% CLOSE TABLE
fprintf(fid,'\\bottomrule \n');
fprintf(fid,'\\end{tabular}');
fclose(fid);

%% Make Table D3

names = {'Permitted mergers', 'Blocked mergers', 'Permitted mergers','Blocked mergers', 'Permitted mergers', 'Blocked mergers', 'Permitted mergers', 'Blocked mergers', 'Permitted mergers', 'Blocked mergers', 'Permitted mergers', 'Blocked mergers'};

fid = fopen('../results/tables/table3D_guidelines_wages.tex','w');

fprintf(fid,'\\begin{tabular}{l @{\\hspace{2em}} c @{\\hspace{1.5em}}c}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,' & \\multicolumn{2}{l}{\\textbf{{------1982 guidelines------ $\\quad$ }}} & \\multicolumn{2}{l}{\\textbf{{------2010 guidelines------}}} \\\\ \n');
fprintf(fid,'\\cmidrule{2-3}\n');
fprintf(fid,'\\cmidrule{4-5}\n');

fprintf(fid,'DOJ/FTC market classification & Moderate & High & Moderate & High \\\\ \n');

fprintf(fid,'Threshold (HHI, $\\Delta$HHI) & (1000,  100) & (1800, 100) & (1500, 100) & (2500, 200) \\\\ \n');
fprintf(fid,' & (1) & (2) & (3) & (4) \\\\ \n');

fprintf(fid,'\\multicolumn{5}{l}{\\textbf{I. Average REG}} \\\\ \n');

for i = 1:2
    fprintf(fid,'%s & %3.2f & %3.2f %3.2f & %3.2f  \\\\ \n', names{i}, tab(i,:));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{5}{l}{\\textbf{II. Change in average $\\mathbf{W}_j$ assuming 1 percent efficiency gain ($\\%%$)}} \\\\ \n');

for i = 3:4
    fprintf(fid,'%s & %3.2f & %3.2f %3.2f & %3.2f\\\\ \n', names{i}, tab(i,:));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{5}{l}{\\textbf{III. Change in average $\\mathbf{W}_j$ assuming 2 percent efficiency gain ($\\%%$)}} \\\\ \n');

for i = 5:6
    fprintf(fid,'%s & %3.2f & %3.2f %3.2f & %3.2f\\\\ \n', names{i}, tab(i,:));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{5}{l}{\\textbf{IV. Change in average $\\mathbf{W}_j$ assuming 3 percent efficiency gain ($\\%%$)}} \\\\ \n');

for i = 7:8
    fprintf(fid,'%s & %3.2f & %3.2f %3.2f & %3.2f\\\\ \n', names{i}, tab(i,:));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{5}{l}{\\textbf{V. Change in average $\\mathbf{W}_j$ assuming 4 percent efficiency gain ($\\%%$)}} \\\\ \n');

for i = 9:10
    fprintf(fid,'%s & %3.2f & %3.2f %3.2f & %3.2f\\\\ \n', names{i}, tab(i,:));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{5}{l}{\\textbf{VI. Change in average $\\mathbf{W}_j$ assuming 5 percent efficiency gain ($\\%%$)}} \\\\ \n');

for i = 11:12
    fprintf(fid,'%s & %3.2f & %3.2f %3.2f & %3.2f\\\\ \n', names{i}, tab(i,:));    
end

% CLOSE TABLE
fprintf(fid,'\\bottomrule \n');
fprintf(fid,'\\end{tabular}');
fclose(fid);


