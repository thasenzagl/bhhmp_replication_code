clear all
close all
% clc
dbstop if error 

% Add functions
addpath(genpath('../functions'));


X=load('../data/xmin_AER_revision_v20r2');
x      = X.xmin;
m_data = X.m_data;
%__________________________________________________________________________
% OPTIONS
options.iter_S_max  = 100       ;  % Max iterations over wage shares
options.Cournot     = 1         ;  % 1=cournot, 0=bertrand
options.tol_S       = 1e-5      ;  % Tolerance for wage shares
options.upsilon     = 0.20      ;  % Adjustment rate of shares in equilibrium solver
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
Delta           = param.Delta ;

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
    return;
    fprintf('Help\n');
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
z_ij  = zeros(options.Mj_max,J); %Productivity draws
for jj = (1:J)
    z_ij(1:Mj(jj),jj)=(c_ij(1:Mj(jj),jj)==1).*(1-Delta)*(Delta*(1-tauC)/((1-tauC*lambdaK)*R)).^(Delta/(1-Delta)).*z_ij_notilde(1:Mj(jj),jj).^(1/(1-Delta))...
                     +(c_ij(1:Mj(jj),jj)==0).*(1-Delta)*(Delta/R).^(Delta/(1-Delta)).*z_ij_notilde(1:Mj(jj),jj).^(1/(1-Delta)) ;
end
%%%%%%STOP EDITS!%%%%%%%%%%%%%%%%

%% 2. SOLVE MODEL
    
%% A. Solve sectoral equilibria for different types of merger policies
merger_type_vec = {'none','random','top','bottom'};

for xxx = (1:4)
merger_type     = merger_type_vec{xxx};
fprintf('xxx = %2i,\t %s\n',xxx,merger_type);

switch merger_type
    case {'none','random'}
        Mjcut           = 2;    % Was 5 in what we had before (has to be at least 2)
    case {'top','bottom'}
        Mjcut           = 2;    % Was 5 in what we had before (has to be at least 2)
end

% Storage for outputs
mu_ij           = zeros(Mj_max,J);      % Markdown
merge_ij        = zeros(Mj_max,J);      % Merge id
for jj=(1:J)     % J regions
    %______________________________________________________________________
    % MERGER
    merge                       = zeros(Mj(jj),1); 
    if Mj(jj)>=Mjcut
        switch merger_type
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            case 'none'             % None
                % Do nothing
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            case 'random'           % Merge two firms at random
                merge(randsample(Mj(jj),2)) = 1;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            case 'top'              % Merge two most productive firms
                [~,zorder]                  = sort(z_ij(1:Mj(jj),jj),'descend');
                merge(zorder(1:2))          = 1;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            case 'bottom'           % Merge two least productive firms
                [~,zorder]                  = sort(z_ij(1:Mj(jj),jj),'ascend');
                merge(zorder(1:2))          = 1;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end
    end
    merge_ij(1:Mj(jj),jj)   = merge;
    %______________________________________________________________________
    z_i                     = z_ij(1:Mj(jj),jj); % Isolate productivity vector in region
    eq0                     = solve_eq_DRS_merge(z_i,param,val,glob,options,merge);
    mu_ij(1:Mj(jj),jj)      = eq0.mu_i;
end

%% Solve general equilibrium
param.sigma                     = 0;
param.varphi                    = glob.varphi;
options.mu_omega_only           = 'N';
switch merger_type 
    case {'none'}
        options.solve_scale     = 'Y';
        options.PE              = 'N';
        out                     = Welfare_closed_form_SM_arnold(mu_ij,z_ij,Mj,[],[],param,glob,options);
        glob.zbar               = out.Ztilde;
        glob.varphibar          = out.varphibar;
        varphibar               = out.varphibar;
        zbar                    = out.Ztilde;
        mu_base                 = out.mu;
        omega_base              = out.omega;
        glob.W                  = out.W;
        glob.N                  = out.N;
    case {'random'}
        options.solve_scale     = 'N';
        options.PE              = 'Y';
        out                     = Welfare_closed_form_SM_arnold(mu_ij,z_ij,Mj,mu_base,omega_base,param,glob,options);
    case {'top','bottom'}
        options.solve_scale     = 'N';
        options.PE              = 'N';
        out                     = Welfare_closed_form_SM_arnold(mu_ij,z_ij,Mj,[],[],param,glob,options);
end

n_ij        = out.n_ij;
w_ij        = out.w_ij;
N           = out.N;
C           = out.C;

wn_j        = sum(w_ij.*n_ij);
swn_ij      = bsxfun(@rdivide,w_ij.*n_ij,wn_j);
swn_j       = wn_j/sum(wn_j);

save(sprintf('../results/matfiles/baseline_merger_%s',merger_type));

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

