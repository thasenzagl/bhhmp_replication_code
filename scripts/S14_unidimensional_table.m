clear all;
close all;

% Add functions
addpath(genpath('../functions'))

arnold_cutoff = true;

%% 2. LOAD RESULTS

load('../results/matfiles/productivity_gains_results.mat')

delta_HHI = zeros(1,size(out_none.s_ij,2));
HHI = zeros(1,size(out_none.s_ij,2));
for jj=1:size(out_none.s_ij,2)
    idx1 = logical(merge_ij(:,jj));
    delta_HHI(1,jj) = sum(out_none.s_ij(idx1, jj))^2 - sum(out_none.s_ij(idx1, jj).^2);
    HHI(1,jj) = sum(out_none.s_ij(idx1, jj))^2 + sum(out_none.s_ij(~idx1, jj).^2);
end

%% 3. FIND ARNOLD CUTOFF
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

delta_HHI = delta_HHI(:,logical(jstar));
HHI = HHI(:,logical(jstar));
Delta_j = out.Delta_j(logical(jstar));

%% 4. INITILIZE TABLE

Delta_grid = [1, 2, 3, 4, 5];
Confidence_grid = [0.3, 0.35, 0.4, 0.45, 0.5];
tab = zeros(size(Delta_grid,2), size(Confidence_grid,2)*2);

%% 5. SOLVE FOR HHI/delta HHI

HHI_grid = linspace(min(HHI), 1, 1000);
deltaHHI_grid = linspace(min(delta_HHI), 0.5, 1000);
Confidence = zeros(size(HHI_grid,2), 1);

for i=1:size(Delta_grid,2)
    for j=1:size(Confidence_grid,2)
        for k=1:size(HHI_grid,2)
            Confidence(k,1) = find_confidence(HHI_grid(k), 0, Delta_grid(i), delta_HHI, HHI, Delta_j);     
        end
        
        start = find(diff(flipud(Confidence)) > 0, 1, 'first');
        Confidence(1:start) = Confidence(start);
        
        idx = find(Confidence >= Confidence_grid(j), 1, 'last');
        if ~isempty(idx)
            tab(i, j) = 10000*HHI_grid(idx);
        else
            tab(i, j) = NaN;
        end
        
    end
end


for i=1:size(Delta_grid,2)
    for j=1:size(Confidence_grid,2)
        for k=1:size(deltaHHI_grid,2)
            Confidence(k,1) = find_confidence(0, deltaHHI_grid(k), Delta_grid(i), delta_HHI, HHI, Delta_j);     
        end
        
        idx = find(Confidence >= Confidence_grid(j), 1, 'last');
        if ~isempty(idx)
            tab(i, j+size(Confidence_grid,2)) = 10000*deltaHHI_grid(idx);
        else
            tab(i, j+size(Confidence_grid,2)) = NaN;
        end
        
    end
end

%% 6. SAVE LATEX TABLE

names = {'1$\%$ Efficiency gain', '2$\%$ Efficiency gain', '3$\%$ Efficiency gain', '4$\%$ Efficiency gain', '5$\%$ Efficiency gain'};

fid = fopen('../results/tables/tabled1_unidimensional_table.tex','w');

fprintf(fid,'\\begin{tabular}{lcccccc}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,' & \\multicolumn{5}{c}{\\textbf{Probability of WS gain}} \\\\ \n');
fprintf(fid,' & 30$\\%%$ & 35$\\%%$ & 40$\\%%$ & 45$\\%%$ & 50$\\%%$ \\\\ \n');
fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{6}{l}{\\textbf{HHI}} \\\\ \n');

for i=1:size(Delta_grid,2)
    fprintf(fid,'%s & %4.0f & %4.0f & %4.0f & %4.0f & %4.0f\\\\ \n', names{i}, tab(i,1:5));    
end

fprintf(fid,'\\midrule\n');

fprintf(fid,'\\multicolumn{6}{l}{\\textbf{$\\Delta$ HHI}} \\\\ \n');

for i=1:size(Delta_grid,2)
    fprintf(fid,'%s & %4.0f & %4.0f & %4.0f & %4.0f & %4.0f \\\\ \n', names{i}, tab(i,6:10));    
end

% CLOSE TABLE
fprintf(fid,'\\bottomrule \n');
fprintf(fid,'\\end{tabular}');
fclose(fid);


