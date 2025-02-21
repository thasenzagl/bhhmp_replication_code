function [baseline_welfare, new_welfare] = compute_welfare(out_none, out_merge, J, glob, param)

    varphibar = out_none.varphibar;
    varphi = glob.varphi;
    theta = param.theta;

    C0 = out_none.C;
    N0 = out_none.N;
    U0 = C0 - (varphibar^(-(1/varphi))) * ((N0^(1+(1/varphi)))/(1+(1/varphi)));
    const = (varphibar^(-(1/varphi))) * N0^(1/varphi) * N0^(-(1/theta)) * (theta/(theta+1));

    w_j = out_merge.w_j;
    n_j = out_merge.n_j;
    w0_j = out_none.w_j;
    n0_j = out_none.n_j;

    % Welfare approximation in percentage changes
    baseline_welfare = zeros(1,J);
    new_welfare = zeros(1,J); 
    for jj = 1:J
        baseline_welfare(1,jj) = (w0_j(1,jj) * n0_j(1,jj)) - const * (n0_j(1,jj)^((theta+1)/theta));
        new_welfare(1,jj) = ((w_j(1,jj) * n_j(1,jj)) - const * (n_j(1,jj)^((theta+1)/theta)));
    end
    
end
