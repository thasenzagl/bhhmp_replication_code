function [table, table_freq, DWP] = make_chart_dwp(bins, merge_ij, Delta_j, Mj, Mjcut, J, out_none, eta, theta, percentile)

    % Compute DWP and generate table
    DWP = zeros(size(out_none.s_ij,2), 2);
    table_temp1 = zeros(size(bins,2), size(bins,2), size(out_none.s_ij,2));
    table_temp2 = zeros(size(bins,2), size(bins,2), size(out_none.s_ij,2));  
    table_freq = zeros(size(bins,2), size(bins,2));
        
    for jj=(1:size(out_none.s_ij,2))
        if (Mj(jj)>=Mjcut)
            merge = merge_ij(1:Mj(jj),jj);
            
            s_j = out_none.s_ij(1:Mj(jj),jj);
            s_j_merge = s_j(logical(merge));

            DWP(jj,1) = compute_dwp(s_j_merge(2), eta, theta);
            DWP(jj,2) = compute_dwp(s_j_merge(1), eta, theta);

            idx1 = find(bins < DWP(jj,1), 1, 'last');
            idx2 = find(bins < DWP(jj,2), 1, 'last');
            
            table_temp1(idx1, idx2, jj) = Delta_j(jj);  
            table_temp2(idx2, idx1, jj) = Delta_j(jj);  
            
            table_freq(idx1, idx2) = table_freq(idx1, idx2) + 1;
            table_freq(idx2, idx1) = table_freq(idx2, idx1) + 1;
            
        end
    end
    
    table_temp = cat(3, table_temp1, table_temp2);
    table_temp(logical(table_temp==0)) = NaN;
    
    table = prctile(table_temp, percentile, 3);

end



