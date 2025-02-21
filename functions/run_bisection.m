function out = run_bisection(variable, tol_bisect, max_iter, out_none, z_ij, merge_ij, J, Mj, Mj_max, Mjcut, bound_upper, bound_lower, param, val, glob, options)

Delta_j  = zeros(1, J);

% Find bounds
[z_ij_LB,  z_ij_UB] = find_bounds(bound_lower, bound_upper, z_ij, merge_ij, Mj, J);   
 
    
% Check that the bounds are large enough
out_rand_LB = solve_model('select', z_ij_LB, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none);
out_rand_UB = solve_model('select', z_ij_UB, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none);

if strcmp(variable, 'n_j') == 1    
    assert(all(out_rand_UB.n_j(logical(Mj>=Mjcut)) > out_none.n_j(logical(Mj>=Mjcut))), 'Need to increase upper bound for bisection.')
    assert(all(out_rand_LB.n_j(logical(Mj>=Mjcut)) < out_none.n_j(logical(Mj>=Mjcut))), 'Need to decrease lower bound for bisection.')
elseif strcmp(variable, 'w_j') == 1    
    assert(all(out_rand_UB.w_j(logical(Mj>=Mjcut)) > out_none.w_j(logical(Mj>=Mjcut))), 'Need to increase upper bound for bisection.')
    assert(all(out_rand_LB.w_j(logical(Mj>=Mjcut)) < out_none.w_j(logical(Mj>=Mjcut))), 'Need to decrease lower bound for bisection.')
elseif strcmp(variable, 'wmerge') == 1    
    assert(all(out_rand_UB.wmerge(logical(Mj>=Mjcut)) > out_none.wmerge(logical(Mj>=Mjcut))), 'Need to increase upper bound for bisection.')
    assert(all(out_rand_LB.wmerge(logical(Mj>=Mjcut)) < out_none.wmerge(logical(Mj>=Mjcut))), 'Need to decrease lower bound for bisection.')   
elseif strcmp(variable, 'nmerge') == 1    
    assert(all(out_rand_UB.nmerge(logical(Mj>=Mjcut)) > out_none.nmerge(logical(Mj>=Mjcut))), 'Need to increase upper bound for bisection.')
    assert(all(out_rand_LB.nmerge(logical(Mj>=Mjcut)) < out_none.nmerge(logical(Mj>=Mjcut))), 'Need to decrease lower bound for bisection.')              
end

% Run bisection
z_ij1 = z_ij;     % Initial guess
metric = 1;
i=1;
while metric > tol_bisect

    out_rand1 = solve_model('select', z_ij1, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none);

    if strcmp(variable, 'n_j') == 1
        resid = out_rand1.n_j - out_none.n_j;
    elseif strcmp(variable, 'w_j') == 1
        resid = out_rand1.w_j - out_none.w_j;
    elseif strcmp(variable, 'wmerge') == 1
        resid = out_rand1.wmerge - out_none.wmerge;        
    elseif strcmp(variable, 'nmerge') == 1
        resid = out_rand1.nmerge - out_none.nmerge;                
    end

    metric = max(abs(resid));

    fprintf('Iteration %i of Bisection: Error is %f and %i/%i markets have converged. \n' ,i, metric, sum(abs(resid)<=tol_bisect)-sum(Mj<Mjcut), sum(Mj>=Mjcut));

    % Loop over the markets
    for jj=(1:J)

        idx = logical(merge_ij(:,jj));

        if ((abs(resid(jj)) > tol_bisect) && (resid(jj) < 0))
            z_ij_LB(idx,jj) = z_ij1(idx,jj);
            z_ij1(idx,jj) = z_ij1(idx,jj) + 0.5*(z_ij_UB(idx,jj)-z_ij1(idx,jj));
        elseif ((abs(resid(jj)) > tol_bisect) && (resid(jj) > 0))   
            z_ij_UB(idx,jj) = z_ij1(idx,jj);                 
            z_ij1(idx,jj) = z_ij1(idx,jj) - 0.5*(z_ij1(idx,jj)-z_ij_LB(idx,jj));
        end
    end

    i = i + 1;

    if i > max_iter
        fprintf("Bisecton did not converge in every market.\n")
        break;
    end
end

Delta_ij = zeros(options.Mj_max, J);
Delta_ij(logical(merge_ij)) = 100*(log(z_ij1(logical(merge_ij)))-log(z_ij(logical(merge_ij))));
for jj=1:J
    if (Mj(jj)>=Mjcut)
        temp = Delta_ij(logical(merge_ij(:,jj)),jj);
        assert(abs(temp(1)-temp(2))<1e-8);
        Delta_j(jj) = mean(temp);
    end
end

%% 4. Store results

out.Delta_j = Delta_j';
out.z_ij = z_ij1;
out.n_ij = out_rand1.n_ij;
out.w_ij = out_rand1.w_ij;
out.s_ij = out_rand1.s_ij;
out.merge_ij = merge_ij;
out.mu_ij = out_rand1.mu_ij;
out.n_j = out_rand1.n_j;
out.w_j = out_rand1.w_j;
%out.nmerge = out_rand1.nmerge;
%out.wmerge = out_rand1.wmerge;
out.Nbodies_j = out_rand1.Nbodies_j;
out.HHIwn_j = out_rand1.HHIwn_j;

end



