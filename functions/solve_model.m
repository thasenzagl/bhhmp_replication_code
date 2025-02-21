function out = solve_model(merger_type, z_ij, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none)

% Storage for outputs
mu_ij           = zeros(Mj_max,J);      % Markdown
s_ij           = zeros(Mj_max,J);       % Shares
for jj=(1:J)     % J regions
    %______________________________________________________________________
    % MERGER
    merge = merge_ij(1:Mj(jj),jj);
    
    %______________________________________________________________________
    z_i                     = z_ij(1:Mj(jj),jj); % Isolate productivity vector in region
    eq0                     = solve_eq_DRS_merge(z_i,param,val,glob,options,merge);
    mu_ij(1:Mj(jj),jj)      = eq0.mu_i;
    s_ij(1:Mj(jj),jj)      = eq0.s_i;   
    
    % Throw error if shares don't converge
    assert(eq0.non_converge == 0, 'Shares for market %i did not converge. Increase iter_S_max. \n' ,jj)
     
end

%% Solve for equilibrium
param.sigma                     = 0;
param.varphi                    = glob.varphi;
options.mu_omega_only           = 'N';
switch merger_type 
    case {'none'}
        options.solve_scale     = 'Y';
        options.PE              = 'N';
        out                     = Welfare_closed_form_SM(mu_ij,z_ij,merge_ij,Mj,Mjcut,[],[],param,glob,options);
        out.s_ij                = s_ij;
        out.z_ij                = z_ij;      
    case {'select'}
        options.solve_scale     = 'N';
        options.PE              = 'Y';
        glob.zbar               = out_none.Ztilde;
        glob.varphibar          = out_none.varphibar;
        mu_base                 = out_none.mu;
        omega_base              = out_none.omega;
        glob.W                  = out_none.W;
        glob.N                  = out_none.N;
        glob.ztilde             = out_none.ztilde;
        out                     = Welfare_closed_form_SM(mu_ij,z_ij,merge_ij,Mj,Mjcut,mu_base,omega_base,param,glob,options);    
        out.s_ij                = s_ij; 
    case {'select_GE'}
        options.solve_scale     = 'N';
        options.PE              = 'N';
        glob.zbar               = out_none.Ztilde;
        glob.varphibar          = out_none.varphibar;        
        out                     = Welfare_closed_form_SM(mu_ij,z_ij,merge_ij,Mj,Mjcut,[],[],param,glob,options);
        out.s_ij                = s_ij;
        out.z_ij                = z_ij;              
end

end
