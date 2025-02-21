clear all
close all

target_var = 'w_j';

% Add functions
addpath(genpath('../functions'));

X=load('../data/xmin_AER_revision_v20r2');
x      = X.xmin;
m_data = X.m_data;
%__________________________________________________________________________
% OPTIONS
options.iter_S_max  = 1e8      ;  % Max iterations over wage shares
options.Cournot     = 1        ;  % 1=cournot, 0=bertrand
options.tol_S       = 1e-12    ;  % Tolerance for wage shares
options.upsilon     = 0.20     ;  % Adjustment rate of shares in equilibrium solver
%__________________________________________________________________________
% OUTPUT OPTIONS
options.print       = 'Y';
options.minwage     = 'N';
options.Sobol       = 'N';
options.plot        = 'N';
options.save_fig    = 'N';
options.derenencourt= 'N';
options.tradeable   = 'N'; %%KFH 3.2.21 Calibrate to tradeables


glob.varphi         = 0.5;
glob.J              = 200000                          ; % Number of markets (j=1,...,J)
if strcmp(options.tradeable,'Y')
    glob.Ndist          = load('../data/Nfit')                  ; % Load frac_1 and Ndist parameters from fit_firm_dist_gp_II_2014.m
else
    glob.Ndist          = load('../data/Nfit_all')                  ; % Load frac_1 and Ndist parameters from fit_firm_dist_gp_II_2014.m
end

options.Mj_max      = max(max(glob.Ndist.M_j))      ; % Cap firms if region has more than X firms -- hard coded in fit_firm_dist_gp_II_2014
    
glob.beta           =1/(1.04)      ; % Discount factor
glob.delta          =.1            ; % Depreciation rate
glob.share_ccorp    =  0.418         ; % table_4_market_stats_v1_2014
glob.tauC_shock_size= 0        ; %KFH 3.2.21 %(b) (C:\Users\kyleh\Dropbox\Market Power\20_AER_Moments\tax\corp_tax_changes_giroud_rauh.txt)
glob.tauC           = 0        ; %KFH 3.2.21 % State corp tax avg rate 
glob.lambdaK 		= 0        ; %KFH 3.2.21 %(d) (C:\Users\kyleh\Dropbox\Market Power\20_AER_Moments\Compustat_debt_assets\1_DEM_Compustat_Moments.do –> target_moment_debt_assets.xls)

 
%Tradeable
if strcmp(options.tradeable,'Y')
    glob.AveFirmSize_data =34.63;
    glob.AveEarnings_data =(2018*1000)/34.63;
else
    glob.AveFirmSize_data =22.83;
    glob.AveEarnings_data =(1000*1000)/22.83;
end


%% A. CALIBRATION

options.Sobol           = 'Y';
options.minwage         = 'N';
options.solve_scale     = 'N'; 

% 1. Initialize storage
% Retrieve number of core that this code is running on
task = getCurrentTask;
if isempty(task)
    xxx=1;
else
    xxx=task.ID;
end
%__________________________________________________________________________
% FIX SEED
 rng('default');
%__________________________________________________________________________
% SETUP / UNPACK
param.eta           = x(1); 
param.theta         = x(2);
param.xi            = x(3);
param.alpha         = x(4);
param.lambdaC       = x(5);
param.Delta         = x(6);

share_ccorp     =glob.share_ccorp          ; %set this parameter in the sampling to yield actual_share_ccorp fraction of firms as ccorps
beta            =glob.beta                 ; %Discount factor
delta           =glob.delta                ; %Depreciation rate
tauC_shock_size =glob.tauC_shock_size      ; %Shock size
tauC            =glob.tauC                 ; %Corporate tax rate
Delta           =param.Delta ;

%Fixed
%KS=glob.KS               ; %Capital share
R=(1/beta) - 1 + delta   ; %Rental rate=r+delta 
%STOP EDITS!%

% Number of firms distribution parameters
param.m_Ndist       = glob.Ndist.m_Ndist;
param.sigma_Ndist   = glob.Ndist.sigma_Ndist;
param.theta_Ndist   = glob.Ndist.theta_Ndist;  
param.frac_1        = glob.Ndist.frac_1;        % Mass at 1

% Global parameters
varphi              = glob.varphi;
Mj_max              = options.Mj_max;
J                   = glob.J;           % Number of markets
val.P               = 1;                % Normalize to 1, final good is numeraire
%__________________________________________________________________________
% Take parameters out of param structure
eta                 = param.eta;
theta               = param.theta;
xi                  = param.xi; 
alpha               = param.alpha;

%Edits!%
lambdaK             = glob.lambdaK            ; %Fraction of capital that is debt financed
lambdaC             = param.lambdaC            ; %CCorp productivity premium
%Edits!%

%__________________________________________________________________________
% BREAK IF THETA>ETA;
if (theta>eta)
    fprintf('Help\n');
    return;
end

%% 1. DRAW RANDOM VARIABLES
% 1. Draw number of firms Mj in each market j=1,...,J 
u       = rand(J,1);
Mj      = zeros(1,J);
for jj = (1:J)
   ilow = find(glob.Ndist.F_Mj<u(jj),1,'last');
   if isempty(ilow)
       ilow = 1;
   end
   Mj(jj) = glob.Ndist.M_j(ilow);
end
Mj              = Mj'; 
Mj(Mj>=Mj_max)  = Mj_max;           % Truncate number of firms in a market -- hard coded in Nfit!!

% Only consider markets with two or more firms
Mjcut = 2;
 
%***Use this random number generator because it fixes the seeds on parfor
%https://www.mathworks.com/help/stats/prob.normaldistribution.random.html
%input parameter a=lower bound, b=upper bound
u_ij            = random('uniform',0,1,[options.Mj_max,J])  ; 
c_ij            = zeros(options.Mj_max,J);
for jj=1:J
    c_ij(1:Mj(jj),jj) = (u_ij(1:Mj(jj),jj)<share_ccorp)      ; 
end

% 2.  Draw firm level productivities z_ij
z_ij_notilde  = zeros(options.Mj_max,J); %Productivity draws
for jj=1:J
    % xi is the standard deviation of the lognormal
    z_ij_notilde(1:Mj(jj),jj) = random('Lognormal',1,xi,[1,Mj(jj)]);
end  
 
%Scale CCorp productivity
z_ij_notilde(c_ij==1)   = (1+lambdaC)*z_ij_notilde(c_ij==1);
 
%Recover the "tilde z_ij" which is ALWAYS z_ij!
%Firm's "tilde z_ij" is productivity after substituting out capital decision
z_ij  = zeros(options.Mj_max,J); %Productivity draws
for jj = (1:J)
    z_ij(1:Mj(jj),jj)=(c_ij(1:Mj(jj),jj)==1).*(1-Delta)*(Delta*(1-tauC)/((1-tauC*lambdaK)*R)).^(Delta/(1-Delta)).*z_ij_notilde(1:Mj(jj),jj).^(1/(1-Delta))...
                     +(c_ij(1:Mj(jj),jj)==0).*(1-Delta)*(Delta/R).^(Delta/(1-Delta)).*z_ij_notilde(1:Mj(jj),jj).^(1/(1-Delta)) ;
end

% Draw merging firms
merge_ij = zeros(Mj_max,J);
for jj=(1:J)   
    %______________________________________________________________________
    % MERGER
    merge = zeros(Mj(jj),1); 
    if Mj(jj) >= Mjcut
        merge(randsample(Mj(jj),2)) = 1;
    end
    merge_ij(1:Mj(jj),jj) = merge;

end

%% 2. SOLVE BASELINE MODEL

% No merger
out_none = solve_model('none', z_ij, zeros(Mj_max,J), J, Mj, Mj_max, Mjcut, param, val, glob, options, []);
out_none.wmerge(Mj>=Mjcut) = sum(reshape(out_none.w_ij(logical(merge_ij)), [2, sum(Mj>=Mjcut)]).^(eta+1)).^(1/(eta+1));
out_none.nmerge(Mj>=Mjcut) = sum(reshape(out_none.n_ij(logical(merge_ij)), [2, sum(Mj>=Mjcut)]).^((eta+1)/eta)).^(eta/(eta+1));

out_merge = solve_model('select', z_ij, merge_ij, J, Mj, Mj_max, Mjcut, param, val, glob, options, out_none);

%% 3. COMPUTE REQUIRED PRODUCTIVITY GAINS

tol_bisect = 1e-8;
max_iter = 100;

bound_lower = 0.1;
bound_upper = 0.5;

% Run bisections
out = run_bisection(target_var, tol_bisect, max_iter, out_none, z_ij, merge_ij, J, Mj, Mj_max, Mjcut, bound_upper, bound_lower, param, val, glob, options);

%% SAVE RESULTS

save('../results/matfiles/productivity_gains_results', '-v7.3');

