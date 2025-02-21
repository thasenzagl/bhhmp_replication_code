function out = solve_eq_DRS_merge(z_i,param,val,glob,options,merge)       
%__________________________________________________________________________
% UNPACK
theta           = param.theta;
eta             = param.eta;
alpha           = param.alpha; 
Cournot         = options.Cournot;
upsilon         = options.upsilon;
%__________________________________________________________________________
% Starting guess (equal shares)
s_i             = z_i/sum(z_i);
for iter_S = (1:options.iter_S_max)
  % Firm level elasticity of labor supply
  s_iold      = s_i;
  stilde_i 	  = s_iold; 	
  if ~(all(merge==0))
		% Replace with the sum for the merging firms
        i1          		= find(merge==1,1,'first');
        i2          		= find(merge==1,1,'last'); 
		stilde_i([i1;i2]) 	= stilde_i(i1) + stilde_i(i2); 	
  end
  % NEW: Use stilde_i to compute shares for epsilon
  eps_i         =   (Cournot==1).*((stilde_i.*(1/theta) + (1-stilde_i).*(1/eta)).^(-1)) + ...
						   (Cournot==0).*(stilde_i.*theta+(1-stilde_i).*eta);
  % Markdown
  mu_i          = eps_i./(eps_i+1);
  %________________________________________________________________________
  % Wage
  a1            = 1/(1+(1-alpha).*theta);
  a2            = -(1-alpha)*(eta-theta)/(eta+1);
  % NEW: Use s_i to compute wages
  w_i           = (mu_i.*alpha.*z_i.*s_i.^a2).^a1; 
  %________________________________________________________________________
  % Sectoral Wage (CES index)
  W_j           = sum(w_i.^(1+eta)).^(1/(1+eta));
  % Implied shares
  s_i_new       = (w_i./W_j).^(1+eta);
  % Distance of shares and new shares
  dist_S        = max(abs(s_i_new-s_iold));
  if (dist_S<options.tol_S)
		break;
  end
  % Update S slowly
  s_i           = upsilon*s_i_new + (1-upsilon)*s_iold;
  s_i           = s_i/sum(s_i); 
end
%__________________________________________________________________________
% Packup output
out.s_i         = s_i;
out.eps_i       = eps_i;
out.mu_i        = mu_i;
out.w_i         = w_i; 
%__________________________________________________________________________
% Diagnostics
out.non_converge    = (iter_S==options.iter_S_max);
out.iterations      = (iter_S);

end