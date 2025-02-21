function confidence = find_confidence(HHI_cutoff, deltaHHI_cutoff, Delta, delta_HHI, HHI, Delta_j)

    % Blocked mergers
    idx_blocked = logical(((HHI > HHI_cutoff) & (delta_HHI > deltaHHI_cutoff)));
    
    % Permitted mergers
    idx = ~idx_blocked;
    
    if sum(idx) ~= 0
        confidence = sum(Delta_j(idx)<=Delta)/sum(idx); 
    else
        confidence = NaN;
    end
end
