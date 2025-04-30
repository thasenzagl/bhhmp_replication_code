function table = make_chart_share_fraction(bins_small, bins_large, merge_ij, Delta_j, Mj, out_none, cutoff)

    table = zeros(size(bins_large,2), size(bins_small,2));
    table_freq = zeros(size(bins_large,2), size(bins_small,2));
    
    for jj=(1:size(out_none.s_ij,2))
        merge = merge_ij(1:Mj(jj),jj);
        s_j = out_none.s_ij(1:Mj(jj),jj);
        s_j_merge = sort(s_j(logical(merge)));

        idx_small = find(bins_small > s_j_merge(1), 1, 'first');
        idx_large = find(bins_large > s_j_merge(2), 1, 'first');

        table_freq(idx_large, idx_small) = table_freq(idx_large, idx_small) + 1;

        if Delta_j(jj) < cutoff
            table(idx_large, idx_small) = table(idx_large, idx_small) + 1;
        end
    end

    table = table ./ table_freq;
    
end
