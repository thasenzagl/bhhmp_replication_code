function D = compute_diversion(s1, n1, n2, params)

    eta = params.eta;
    theta = params.theta;
    
    D = (((eta-theta)*s1)/(eta - (eta - theta)*s1)) * (n2/n1);

end