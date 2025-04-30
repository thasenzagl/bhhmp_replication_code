function [HHI_fig, delta_HHI] = make_chart_hhi_delta(out, out_none, merge_ij, HHI_bins, delta_HHI_bins, percentile)

    % Compute HHI
    delta_HHI = zeros(1,size(out_none.s_ij,2));
    HHI = zeros(1,size(out_none.s_ij,2));
    for jj=1:size(out_none.s_ij,2)
        idx1 = logical(merge_ij(:,jj));
        delta_HHI(1,jj) = sum(out_none.s_ij(idx1, jj))^2 - sum(out_none.s_ij(idx1, jj).^2);
        HHI(1,jj) = sum(out_none.s_ij(idx1, jj))^2 + sum(out_none.s_ij(~idx1, jj).^2);
    end
    
    % Table for heatmap figure
    HHI_fig = zeros(size(HHI_bins,2)-1, size(delta_HHI_bins,2)-1);
    for row = 2:size(HHI_bins,2)
        for col = 2:size(delta_HHI_bins,2)
            idx = (HHI >= (HHI_bins(row-1)) & ...
                    (HHI < HHI_bins(row)) & ...
                    (delta_HHI >= delta_HHI_bins(col-1)) & ...
                    (delta_HHI < delta_HHI_bins(col)));
            
            if sum(idx) == 0
                HHI_fig(row-1, col-1) = NaN;                
            else
                HHI_fig(row-1, col-1) = prctile(out.Delta_j(logical(idx)), percentile); 
            end
            
        end
    end
    
end