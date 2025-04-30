function table = make_chart_hhi_fraction(HHI_bins, delta_HHI_bins, merge_ij, Delta_j, out_none, cutoff)

    % Compute HHI
    delta_HHI = zeros(1,size(out_none.s_ij,2));
    HHI = zeros(1,size(out_none.s_ij,2));
    for jj=1:size(out_none.s_ij,2)
        idx1 = logical(merge_ij(:,jj));
        delta_HHI(1,jj) = sum(out_none.s_ij(idx1, jj))^2 - sum(out_none.s_ij(idx1, jj).^2);
        HHI(1,jj) = sum(out_none.s_ij(idx1, jj))^2 + sum(out_none.s_ij(~idx1, jj).^2);
    end
    
    table = zeros(size(HHI_bins,2)-1, size(delta_HHI_bins,2)-1);
    for row = 2:size(HHI_bins,2)
        for col = 2:size(delta_HHI_bins,2)
            idx = (HHI >= (HHI_bins(row-1)) & ...
                    (HHI < HHI_bins(row)) & ...
                    (delta_HHI >= delta_HHI_bins(col-1)) & ...
                    (delta_HHI < delta_HHI_bins(col)));
            
            if sum(idx) == 0
                table(row-1, col-1) = NaN;
            else
                table(row-1, col-1) = 1-(sum(Delta_j(logical(idx))>=cutoff)/sum(idx)); 
            end
            
        end
    end
    
end
