function out = Welfare_closed_form_SM_arnold(mu_ij,z_ij,Mj,mu1,omega1,param,glob,options)

out = [];
try 
    PE = strcmp(options.PE,'Y');
catch
    PE = false;
end
try 
    mu_omega_only = strcmp(options.mu_omega_only,'Y');
catch
    mu_omega_only = false;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Note: {mu1,omega1} are the aggregate markdown and wedge from the economy
% that we are comparing the baseline economy to.
% FOR EXAMPLE
% If the baseline economy is oligopsony, and we want to compare to the
% merger economy, then we need to (i) compute {mu1,omega1} from the
% *merger* economy, (ii) feed these in here to compute 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Unpack:
eta                     = param.eta; 
theta                   = param.theta;
alpha                   = param.alpha;
Delta                   = param.Delta;
sigma                   = param.sigma;
varphi                  = param.varphi;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta                    = glob.beta;
delta                   = glob.delta;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Rename 
ztilde_ij               = z_ij;
alphatilde              = alpha;
alpha                   = (1-Delta)*alphatilde + Delta;
gamma                   = 1-Delta/alpha;

%% CONSTRUCT WEDGES: {ztilde, mu, omega}
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Productivity: ztilde_j, ztilde
coeff_eta               = (1+eta)/(1+eta*(1-alphatilde));
ztilde_j                = sum(ztilde_ij.^coeff_eta).^(1/coeff_eta);
coeff_theta             = (1+theta)/(1+theta*(1-alphatilde));
ztilde                  = sum(ztilde_j.^coeff_theta).^(1/coeff_theta);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Markdowns: mu_j, mu
zratio_ij               = bsxfun(@rdivide,ztilde_ij,ztilde_j);
zratio_ij(z_ij==0)      = 0;
zratio_j                = bsxfun(@rdivide,ztilde_j,ztilde);
zratio_j(ztilde_j==0)   = 0;
mu_j                    = sum(zratio_ij.^coeff_eta .* mu_ij.^coeff_eta).^(1/coeff_eta);  
mu                      = sum(zratio_j.^coeff_theta .* mu_j.^coeff_theta).^(1/coeff_theta);
if PE
    % If considering a partial equilibrium experiment then take the
    % previous equilibrium mu
    mu      = mu1; 
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Omegas: omega_j, omega
muratio_ij              = bsxfun(@rdivide,mu_ij,mu_j);
muratio_ij(z_ij==0)     = 0;
muratio_j               = bsxfun(@rdivide,mu_j,mu);
muratio_j(ztilde_j==0)  = 0;
coeff_eta2              = eta*alphatilde/(1+eta*(1-alphatilde));
coeff_theta2            = theta*alphatilde/(1+theta*(1-alphatilde));
omega_j                 = sum(zratio_ij.^coeff_eta .* muratio_ij.^coeff_eta2);  
omega                   = sum(zratio_j.^coeff_theta .* muratio_j.^coeff_theta2 .* omega_j);
if PE
    % If considering a partial equilibrium experiment then take the
    % previous equilibrium omega
    omega      = omega1; 
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out.mu              = mu;
out.omega           = omega;
out.ztilde          = ztilde;
if (mu_omega_only)
    return
end

%% Solve baseline economy
%__________________________________________________________________________
% A. SOLVE GE: {W,N,Y,K,C,R,Ztilde,varphibar}
% Common exogenous terms
R                       = (1/beta)-(1-delta);
s_c                     = (1- (delta/R)*(1-gamma)*alpha);
stilde_c                = s_c/(1-(1-gamma)*alpha);

if strcmp(options.solve_scale,'Y')
    %______________________________________________________________________
    % Compute nu_ij
    coeff_theta_nu      = theta/(1 + theta*(1-alphatilde));
    coeff_eta_nu        = eta/(1 + eta*(1-alphatilde));
    nu_j                = ((mu_j.*ztilde_j)./(mu.*ztilde)).^coeff_theta_nu;
    nu_ij               = bsxfun(@times,bsxfun(@rdivide,mu_ij.*ztilde_ij,mu_j.*ztilde_j).^coeff_eta_nu,nu_j);
    nu_ij(mu_ij==0)     = 0;
    %______________________________________________________________________
    % [SYSTEM 1] 
    W                   = glob.AveEarnings_data*sum(sum(nu_ij));
    N                   = glob.AveFirmSize_data*sum(Mj)/sum(sum(nu_ij));
    Ytilde              = W*N/(alphatilde*mu/omega);
    C                   = stilde_c*Ytilde;
    Y                   = Ytilde/(1-(1-gamma)*alpha);
    K                   = (1-gamma)*alpha*Y/R;
    Ztilde              = Ytilde/(omega*ztilde*N^alphatilde);
    varphibar           = N/(W^varphi * C^(-sigma*varphi));
else
    %______________________________________________________________________
    % Load in Ztilde and varphibar
    Ztilde              = glob.zbar;
    varphibar           = glob.varphibar;
    %______________________________________________________________________
    % [SYSTEM 2]
    N       = Nfun(omega,mu,varphibar,stilde_c,sigma,varphi,alphatilde,Ztilde,ztilde);
    W       = mu*Ztilde*ztilde*alphatilde*N^(alphatilde-1);
    Ytilde  = omega*Ztilde*ztilde*N^alphatilde;
    C       = stilde_c*Ytilde;
    Y       = Ytilde/(1-(1-gamma)*alpha);
    K       = ((1-gamma)*alpha)*Y/R;
end
%__________________________________________________________________________
% PARTIAL EQUILIBRIUM
if (PE)
    N       = glob.N; 
    W       = glob.W;
    Ytilde  = [];
    C       = [];
    Y       = [];
    K       = [];
end
%__________________________________________________________________________
% B. Compute firm level variables
% 1. Market
coeff_market            = 1/(1+theta*(1-alphatilde));
n_j                     = (mu_j.*alphatilde.*Ztilde.*ztilde_j./W).^(theta*coeff_market)*N^coeff_market;
w_j                     = mu_j.*alphatilde.*Ztilde.*ztilde_j.*n_j.^(alphatilde-1);
% 2. Firm
coeff_firm              = 1/(1+eta*(1-alphatilde));
n_ij                    = bsxfun(@times,bsxfun(@rdivide,mu_ij.*alphatilde.*Ztilde.*ztilde_ij,w_j).^(eta*coeff_firm),n_j.^coeff_firm);
w_ij                    = mu_ij.*alphatilde.*Ztilde.*ztilde_ij.*n_ij.^(alphatilde-1);
n_ij(mu_ij==0)          = 0;
w_ij(mu_ij==0)          = 0;
%__________________________________________________________________________
% C. CHECK AGGREGATION
w_j1        = sum(w_ij.^(eta+1)).^(1/(eta+1));
n_j1        = sum(n_ij.^((eta+1)/eta)).^(eta/(eta+1));
W1          = sum(w_j1.^(theta+1)).^(1/(theta+1));
N1          = sum(n_j1.^((theta+1)/theta)).^(theta/(theta+1));
if ~PE
    assert(abs((W1-W)/W)<1e-9);
    assert(abs((N1-N)/N)<1e-9);
end
%__________________________________________________________________________
% D. OUTPUT
out.varphibar       = varphibar;
out.Ztilde          = Ztilde;
out.Y               = Y;
out.K               = K;
out.N               = N;
out.C               = C;
out.W               = W;
out.Nbodies         = sum(sum(n_ij));
out.w_ij            = w_ij;
out.n_ij            = n_ij;

out.AveWage         = sum(sum(w_ij.*n_ij))/sum(sum(n_ij));
out.AggEmp          = sum(sum(n_ij));

swn_ij              = bsxfun(@rdivide,w_ij.*n_ij,sum(w_ij.*n_ij));
swn_j               = sum(w_ij.*n_ij)/sum(sum(w_ij.*n_ij));
HHIwn_j             = sum(swn_ij.^2);
HHIwn_uw            = mean(HHIwn_j);
HHIwn_w             = sum(swn_j.*HHIwn_j);

out.HHIwn_uw        = HHIwn_uw;
out.HHIwn_w         = HHIwn_w;


%% COMPETITIVE / ALTERNATIVE MODEL
if ~isempty(mu1) && ~PE

% A. SOLVE GE: {Nc,Wc,Cc,Yc,Kc}
Nc              = Nfun(omega1,mu1,varphibar,stilde_c,sigma,varphi,alphatilde,Ztilde,ztilde);
Wc              = mu1*Ztilde*ztilde*alphatilde*Nc^(alphatilde-1);
Ytildec         = omega1*Ztilde*ztilde*Nc^alphatilde;
Cc              = stilde_c*Ytildec;
Yc              = Ytildec/(1-(1-gamma)*alpha);
Kc              = (1-gamma)*alpha*Yc/R;

% B. Welfare
% CH: I think sigma=0 nests the GHH proference, so can get rid of the
% original lines that consider two cases separately
xxx        = 1+1/varphi;
lambdafun  = @(C0,Cc,N0,Nc) (1/C0)*(Cc.^(1-sigma) - (1-sigma)*varphibar^(-1/varphi)*(Nc^xxx - N0^xxx)/xxx).^(1/(1-sigma)) - 1;

lambdaE             = lambdafun(C,Cc,N,Nc);

% C. Compute competitive firm level variables
muc_ij              = mu_ij;    
muc_ij(mu_ij>0)     = mu1;      % Efficient 
muc_j               = mu_j;
muc_j(mu_j>0)       = mu1;      % Efficient
coeff               = 1/(1+theta*(1-alphatilde));
nc_j                = (muc_j.*alphatilde.*Ztilde.*ztilde_j./Wc).^(theta*coeff)*Nc^coeff;
wc_j                = muc_j.*alphatilde.*Ztilde.*ztilde_j.*nc_j.^(alphatilde-1);
% 2. Firm
coeff               = 1/(1+eta*(1-alphatilde));
nc_ij               = bsxfun(@times,bsxfun(@rdivide,muc_ij.*alphatilde.*Ztilde.*ztilde_ij,wc_j).^(eta*coeff),nc_j.^coeff);
wc_ij               = muc_ij.*alphatilde.*Ztilde.*ztilde_ij.*nc_ij.^(alphatilde-1);
nc_ij(muc_ij==0)    = 0;
wc_ij(muc_ij==0)    = 0;

% D. Only change mu -> mu1
Nc_temp         = Nfun(omega,mu1,varphibar,stilde_c,sigma,varphi,alphatilde,Ztilde,ztilde);
Cc_temp         = stilde_c*(omega*Ztilde*ztilde*Nc_temp^alphatilde);
lambdaE_mu      = lambdafun(C,Cc_temp,N,Nc_temp);

% E. Only change omega -> omega1
Nc_temp         = Nfun(omega1,mu,varphibar,stilde_c,sigma,varphi,alphatilde,Ztilde,ztilde);
Cc_temp         = stilde_c*(omega1*Ztilde*ztilde*Nc_temp^alphatilde);
lambdaE_omega   = lambdafun(C,Cc_temp,N,Nc_temp);

% F. Check aggregation
wc_j1        = sum(wc_ij.^(eta+1)).^(1/(eta+1));
nc_j1        = sum(nc_ij.^((eta+1)/eta)).^(eta/(eta+1));
Wc1          = sum(wc_j1.^(theta+1)).^(1/(theta+1));
Nc1          = sum(nc_j1.^((theta+1)/theta)).^(theta/(theta+1));
assert(abs((Wc1-Wc)/Wc)<1e-9);
assert(abs((Nc1-Nc)/Nc)<1e-9);

% G. OUTPUT - Competitive economy
out.lambdaE         = lambdaE;
out.lambdaE_mu      = lambdaE_mu;
out.lambdaE_omega   = lambdaE_omega;

out.Yc              = Yc;
out.Kc              = Kc;
out.Nc              = Nc;
out.Cc              = Cc;
out.Wc              = Wc;
out.Nbodiesc        = sum(sum(nc_ij));
out.wc_ij           = wc_ij;
out.nc_ij           = nc_ij;

out.AveWagec        = sum(sum(wc_ij.*nc_ij))/sum(sum(nc_ij));
out.AggEmpc         = sum(sum(nc_ij));

swnc_ij             = bsxfun(@rdivide,wc_ij.*nc_ij,sum(wc_ij.*nc_ij));
swnc_j              = sum(wc_ij.*nc_ij)/sum(sum(wc_ij.*nc_ij));
HHIwnc_j            = sum(swnc_ij.^2);
HHIwnc_uw           = mean(HHIwnc_j);
HHIwnc_w            = sum(swnc_j.*HHIwnc_j);

out.HHIwnc_uw       = HHIwnc_uw;
out.HHIwnc_w        = HHIwnc_w;

%% WELFARE APPROXIMATION
dlogmu              = (mu1-mu)/mu;
dlogomega           = (omega1-omega)/omega;
psi                 = alpha*gamma*(mu/omega)*(1/s_c);
denom               = 1+varphi*(1-alphatilde)+varphi*sigma*alphatilde;
C_mu                = ((alphatilde - psi)*varphi)/denom;
C_omega             = (1+varphi*(1-alphatilde)+psi*sigma*varphi)/denom;
lambdaA             = C_mu*dlogmu + C_omega*dlogomega;
lambdaA_mu          = C_mu*dlogmu;
lambdaA_omega       = C_omega*dlogomega;

out.lambdaA         = lambdaA;
out.lambdaA_mu      = lambdaA_mu;
out.lambdaA_omega   = lambdaA_omega;

end

%% NESTED FUNCTIONS
% Closed form solution for N
function N = Nfun(omega_in,mu_in,varphibar,stilde_c,sigma,varphi,alphatilde,Ztilde,ztilde)
    coeff   = 1/(1+varphi*(1-alphatilde)+sigma*varphi*alphatilde);
    N       = varphibar*(stilde_c*omega_in)^(-sigma*varphi)*(alphatilde*mu_in)^(varphi)*(Ztilde*ztilde)^((1-sigma)*varphi);
    N       = N^coeff;
end

end












    