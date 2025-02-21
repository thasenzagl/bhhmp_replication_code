function table = make_chart_share_delta(bins, merge_ij, Delta_j, Mj, Mjcut, J, out_none, percentile)

    table_temp = zeros(size(bins,2), size(bins,2), size(out_none.s_ij,2));
    for jj=(1:size(out_none.s_ij,2))
        if (Mj(jj)>=Mjcut)
            merge = merge_ij(1:Mj(jj),jj);
            s_j = out_none.s_ij(1:Mj(jj),jj);
            s_j_merge = sort(s_j(logical(merge)));

            idx_small = find(bins > s_j_merge(1), 1, 'first');
            idx_large = find(bins > s_j_merge(2), 1, 'first');

            table_temp(idx_large, idx_small, jj) = Delta_j(jj);
        end
    end

    table_temp(logical(table_temp==0)) = NaN;
    
    table = prctile(table_temp, percentile, 3);
    
end
