function [z_ij_LB,  z_ij_UB] = find_bounds(bound_lower, bound_upper, z_ij, merge_ij, Mj, J)

    z_ij_LB = z_ij; % lower bound for z
    z_ij_UB = z_ij; % upper bound for z

    for jj=(1:J)   

        merge = merge_ij(1:Mj(jj),jj);
        merge_idx = find(merge==1);

        % Upper and lower bounds for the bisection
        z_ij_LB(merge_idx, jj) = z_ij(merge_idx,jj)*exp(-bound_lower);
        z_ij_UB(merge_idx, jj) = z_ij(merge_idx,jj)*exp(+bound_upper);    
    end
end