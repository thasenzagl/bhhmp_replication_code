%% TABLES

clear;
clc;
dbstop if error

run('../functions/baseline_merger_arnold')

% Only 29% of mergers are within same CZ
factor                          = 0.29;  

%% ARNOLD DATA
arnold.mean_dlogn_Y_w           = -0.144;
arnold.mean_dlogn_Y_uw          = -0.081;
arnold.mean_dlogwn_Y_w          = -0.121;

arnold.mean_dlogn_Y_uw_small_1  = -0.055;
arnold.mean_dlogn_Y_uw_big_1    = -0.162;

arnold.mean_dlogw_Y             = -0.008;

arnold.mean_dlogw_HI            = -0.031;
arnold.mean_dlogw_MI            = -0.008;
arnold.mean_dlogw_LI            = -0.005;

arnold.b_dlogw_dlogHHI          = -0.085;
arnold.b_dlogn_dlogHHI          =  0.309; 
% arnold.b_dlogn_dlogHHI          =  -0.20; 

arnold.b_hhi                    = 0.834;    % Table 8

arnold.dlogHHI                  = 0.17;     % Figure 8A

arnold.mean_n_pre_Y             = 250;
arnold.median_n_pre_Y           = 116;
arnold.mean_w_pre_Y             = 43940;

arnold.b_elast1_w               = -0.147;
arnold.b_elast2_w               = -0.176;
arnold.b_elast3_w               =  0.058;
arnold.b_elast1_n               = -0.450;   % By email
            
results_arnold     = [ arnold.mean_dlogn_Y_w,    arnold.mean_dlogn_Y_uw,     arnold.mean_dlogwn_Y_w ,    ... % Table 3
                arnold.mean_dlogn_Y_uw_small_1 , arnold.mean_dlogn_Y_uw_big_1  ,          ... % Table 4
                arnold.mean_dlogw_Y    ,                                           ... % Table 5
                arnold.mean_dlogw_HI     , arnold.mean_dlogw_MI,  arnold.mean_dlogw_LI,                   ... % Table 6 
                arnold.b_dlogw_dlogHHI     , arnold.b_dlogn_dlogHHI,                      ... % Table 7
                arnold.b_hhi ,                                                     ... % Table 8
                arnold.dlogHHI ,    ... % Figure 8A
                arnold.mean_n_pre_Y , arnold.median_n_pre_Y , arnold.mean_w_pre_Y/1000,               ... % Summary stats
                0,0,0, ... % welfare1
                0,0,0,0,0, ... % welfare2
                arnold.b_elast1_w, arnold.b_elast2_w, arnold.b_elast3_w, arnold.b_elast1_n,...
                ];
            
names_list = cell(14,1);
names_list{1} = 'Change in log employment (weighted)';
names_list{2} = 'Change in log employment (unweighted)';
names_list{3} = 'Change in log payroll (weighted)';
names_list{4} = 'Change in log employment (weighted, small)';
names_list{5} = 'Change in log employment (weighted, big)';
names_list{6} = 'Change in log worker earnings';
names_list{7} = '... high impact market';
names_list{8} = '... medium impact market';
names_list{9} = '... low impact market';
names_list{10}= '$\Delta\log{\overline{w}_j} = \alpha + \beta\Delta\log{HHI_j}$';
names_list{11}= '$\Delta\log{\overline{n}_j} = \alpha + \beta\Delta\log{HHI_j}$';
names_list{12}= '$\Delta HHI_j = \alpha + \beta \Delta \widehat{HHI}_j$';
names_list{13}= 'Average change in log $HHI_j$';

names_list{14}= 'Mean employment pre-merger';
names_list{15}= 'Median employment pre-merger';
names_list{16}= 'Mean worker earnings pre-merger (2011 dollars)';

names_list{17}= 'Welfare: $100\times\lambda_{SS}$';
names_list{18}= 'Change in concentration (weighted): $\Delta HHI^{wn}$';
names_list{19}= 'Change in concentration (unweighted) $\Delta HHI^{wn}$';

names_list{20}= 'Welfare: $100\times\lambda_{SS}$';
names_list{21}= 'Fraction due to change in only $\mu$';
names_list{22}= 'Fraction due to change in only $\omega$';
names_list{23}= 'Change in $\omega$ ($100\times \Delta\log\omega$)';
names_list{24}= 'Change in $\mu$ ($100\times \Delta\log\mu$)';

names_list{25}= 'Elasticity of market wage to $HHI$';
names_list{26}= '... above median $HHI$';
names_list{27}= '... below median $HHI$';
names_list{28}= 'Elasticity of market emp. to $HHI$';

% lambdaE, lambdaE_mu/lambdaE, lambdaE_omega/lambdaE, mu_merger, omega_merger, ...

Table_list = cell(14,1);
Table_list{1} = 'Table 3';
Table_list{2} = 'Table 3';
Table_list{3} = 'Table 3';
Table_list{4} = 'Table 4';
Table_list{5} = 'Table 4';
Table_list{6} = 'Table 5';
Table_list{7} = 'Table 6';
Table_list{8} = 'Table 6';
Table_list{9} = 'Table 6';
Table_list{10}= 'Table 7';
Table_list{11}= 'Table 7';
Table_list{12}= 'Table 8';
Table_list{13}= 'Figure 8A';

Table_list{14}= 'Table 1';
Table_list{15}= 'Table 1';
Table_list{16}= 'Table 1';

Table_list{17}= '';
Table_list{18}= '';
Table_list{19}= '';

Table_list{20}= '';
Table_list{21}= '';
Table_list{22}= '';
Table_list{23}= '';
Table_list{24}= '';

Table_list{25}= 'Table 10(3)';
Table_list{26}= 'Table 10(6)';
Table_list{27}= 'Table 10(6)';
Table_list{28}= 'Unpublished';

%% MERGER - NONE
KS              = 0.1800;
%__________________________________________________________________________
%load('baseline_merger_new_none','z_ij','n_ij','w_ij','tauC','KS','beta','eta','theta','zbar','varphibar','alpha',...
%    'lambdaK','delta','varphi','Delta','KS','Mj','Mjcut','glob','param');
load('../results/matfiles/baseline_merger_none','z_ij','n_ij','w_ij','tauC','beta','eta','theta','zbar','varphibar','alpha',...
    'lambdaK','delta','varphi','Delta','Mj','Mjcut','glob','param');
% This is the same as `baseline_no_tax' since merger experiments ran under
% no corporate taxes
%__________________________________________________________________________
clear hhi*
X0              = load('../results/matfiles/baseline_merger_none');
n_ij0           = X0.n_ij;
w_ij0           = X0.w_ij;
N0_index        = X0.N;
N0_bodies       = sum(sum(n_ij0));
W0_mean         = sum(sum(w_ij0.*n_ij0))/N0_bodies;
swn_ij          = X0.swn_ij;
hhi_j           = sum(swn_ij.^2);

% CPI_2014/CPI_2011
dcpi            = 1.062177883;
hhi_j0          = hhi_j;

% glob
% return

%% MERGER - TOP TWO
fprintf('Top 2\n');
n_ij0           = X0.n_ij;
w_ij0           = X0.w_ij/dcpi;
z_ij            = X0.z_ij;
mu_ij0          = X0.mu_ij;
X               = load('../results/matfiles/baseline_merger_top');
nc_ij           = X.n_ij;
wc_ij           = X.w_ij;
muc_ij          = X.mu_ij;

merge_ij        = X.merge_ij;
merge_j         = sum(merge_ij)>0;
n_ij0           = n_ij0(:,merge_j);
w_ij0           = w_ij0(:,merge_j);
mu_ij0          = mu_ij0(:,merge_j);
z_ij0           = z_ij(:,merge_j);
nc_ij           = nc_ij(:,merge_j);
wc_ij           = wc_ij(:,merge_j);
muc_ij          = muc_ij(:,merge_j);
zc_ij           = X0.z_ij(:,merge_j);
merge_ij        = merge_ij(:,merge_j);
swn_ij          = X.swn_ij(:,merge_j);

% return

Nc_bodies       = sum(sum(nc_ij));
Wc_mean         = sum(sum(wc_ij.*nc_ij))/Nc_bodies;

% Compute new HHIs
Y_ij            = logical(merge_ij) & (n_ij0>0);     % YES - merge
N_ij            = logical(~merge_ij) & (n_ij0>0);    % NO  - not merge
swn_ij_Y        = sum(Y_ij.*swn_ij);    % Combined share of merged firms
swn_ij_N        = N_ij.*swn_ij;         % Shares of non-merged firms
swn_ij_new      = [swn_ij_Y;swn_ij_N];  % Concatenate
hhi_j1          = sum(swn_ij_new.^2);

% Table 8 - Compute counterfactual initial HHIs with just firm removal
swn_ij0             = X0.swn_ij(:,merge_j);
swn_ij0_Y           = sum(Y_ij.*swn_ij0);                               % Combined share of merged firms
swn_ij0_N           = N_ij.*swn_ij0;                                    % Shares of non-merged firms
swn_ij0_new         = [swn_ij0_Y;swn_ij0_N];                            % Concatenate
swn_ij0_new         = bsxfun(@rdivide,swn_ij0_new,sum(swn_ij0_new));    % Renormalize
hhi_j_new           = sum(swn_ij0_new.^2);
dhhi_j_hat          = hhi_j_new - hhi_j0(merge_j);
dhhi_j_actual       = hhi_j1    - hhi_j0(merge_j);

Nmerge              = numel(dhhi_j_actual);
Ntotal              = floor(Nmerge/factor);
Zerovec             = zeros(Ntotal-Nmerge,1);
Onevec              = ones(Ntotal-Nmerge,1);
b_hhi               = regress([dhhi_j_actual';Zerovec],[[ones(size(dhhi_j_hat))';Onevec],[dhhi_j_hat';Zerovec]]);
b_hhi               = b_hhi(2);

n_ij0_Y_vec     = reshape(n_ij0(Y_ij),sum(sum(Y_ij)),1);
n_ij0_N_vec     = reshape(n_ij0(N_ij),sum(sum(N_ij)),1);
w_ij0_Y_vec     = reshape(w_ij0(Y_ij),sum(sum(Y_ij)),1);
w_ij0_N_vec     = reshape(w_ij0(N_ij),sum(sum(N_ij)),1);

n_ij1_Y_vec     = reshape(nc_ij(Y_ij),sum(sum(Y_ij)),1);
n_ij1_N_vec     = reshape(nc_ij(N_ij),sum(sum(N_ij)),1);
w_ij1_Y_vec     = reshape(wc_ij(Y_ij),sum(sum(Y_ij)),1);
w_ij1_N_vec     = reshape(wc_ij(N_ij),sum(sum(N_ij)),1);

mean_n_pre_Y    = mean(n_ij0_Y_vec);
mean_n_pre_N    = mean(n_ij0_N_vec);
median_n_pre_Y  = prctile(n_ij0_Y_vec,50);
median_n_pre_N  = prctile(n_ij0_N_vec,50);
mean_w_pre_Y    = mean(w_ij0_Y_vec);
mean_w_pre_N    = mean(w_ij0_N_vec);

dlogn_ij_Y      = log(n_ij1_Y_vec)-log(n_ij0_Y_vec);
dlogw_ij_Y      = log(w_ij1_Y_vec)-log(w_ij0_Y_vec);
dlogwn_ij_Y     = log(w_ij1_Y_vec.*n_ij1_Y_vec) - log(w_ij0_Y_vec.*n_ij0_Y_vec);

omega_ij_Y      = n_ij0_Y_vec/sum(n_ij0_Y_vec);
ones_ij_Y       = ones(size(n_ij0_Y_vec))/sum(ones(size(n_ij0_Y_vec)));

% Table 3
mean_dlogn_Y_w  = sum(omega_ij_Y.*dlogn_ij_Y);
mean_dlogn_Y_uw = sum(ones_ij_Y.*dlogn_ij_Y);
mean_dlogwn_Y_w = sum(omega_ij_Y.*dlogwn_ij_Y);

% Table 4 - Weighted
i_big                   = (n_ij0_Y_vec>prctile(n_ij0_Y_vec,50));
% i_big                   = (n_ij0_Y_vec>prctile(n_ij0(n_ij0>0),50));
omega_ij_Y_small        = n_ij0_Y_vec(~i_big)/sum(n_ij0_Y_vec(~i_big));
mean_dlogn_Y_w_small    = sum(omega_ij_Y_small.*dlogn_ij_Y(~i_big));
omega_ij_Y_big          = n_ij0_Y_vec(i_big)/sum(n_ij0_Y_vec(i_big));
mean_dlogn_Y_w_big      = sum(omega_ij_Y_big.*dlogn_ij_Y(i_big));

% Table 5 - Earnings 
% Here we weight by the number of workers, since its estimated by Arnold i
% worker level data. Therefore, e.g. if there are 100 workers in one firm
% and 200 workers in another firm, he would have 300 observations.
% Weighting by employment blows up the observations at the bigger firms,
% which is exactly what worker level data does.
mean_dlogw_Y            = sum(omega_ij_Y.*dlogw_ij_Y);  % (-0.008)

% Table 6 - By Initial concentration
dhhi_j              = log(hhi_j1) - log(hhi_j0(merge_j));
dhhi_cut            = prctile(dhhi_j_hat,75);
LI_j                = (dhhi_j_hat < dhhi_cut);              % Below top quantile for dhhi
hhi_j0_merge        = hhi_j0(merge_j);
hhi_cut             = prctile(hhi_j0_merge(~LI_j),50);
MI_j                = (~LI_j) & (hhi_j0_merge<hhi_cut);     % Above top quantile for dhhi & below median hhi in these
HI_j                = (~LI_j) & (hhi_j0_merge>hhi_cut);     % Above top quantile for dhhi & above median hhi in these

% (a) High impact
n_ij1_conc          = nc_ij(:,HI_j);
w_ij0_conc          = w_ij0(:,HI_j);
w_ij1_conc          = wc_ij(:,HI_j );
Y_ij_conc           = Y_ij(:,HI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_HI    = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_HI       = sum(omega_ij_conc.*dlogw_ij_conc_HI);    % (-0.031)

% (b) Medium impact
n_ij1_conc          = nc_ij(:,MI_j);
w_ij0_conc          = w_ij0(:,MI_j);
w_ij1_conc          = wc_ij(:,MI_j );
Y_ij_conc           = Y_ij(:,MI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_MI       = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_MI       = sum(omega_ij_conc.*dlogw_ij_conc_MI);    % (-0.031)

% (c) Low impact
n_ij1_conc          = nc_ij(:,LI_j);
w_ij0_conc          = w_ij0(:,LI_j);
w_ij1_conc          = wc_ij(:,LI_j );
Y_ij_conc           = Y_ij(:,LI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_LI    = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_LI       = sum(omega_ij_conc.*dlogw_ij_conc_LI);    % (-0.031)

% Table 7 - Correlations
dlogn_j             = log(sum(nc_ij))-log(sum(n_ij0));
wbar_j0             = sum(n_ij0.*w_ij0)./sum(n_ij0);
wbar_j1             = sum(nc_ij.*wc_ij)./sum(nc_ij);
dlogw_j             = log(wbar_j1) - log(wbar_j0); 
dlogHHI_j           = log(hhi_j1) - log(hhi_j0(merge_j));

ones_j              = ones(size(dlogw_j));

Nmerge              = numel(dlogw_j);
Ntotal              = floor(Nmerge/factor);
Zerovec             = zeros(Ntotal-Nmerge,1);
Onevec              = ones(Ntotal-Nmerge,1);

% Unweighted
b_dlogw_dlogHHI     = regress([dlogw_j';Zerovec],[[ones_j';Onevec],[dlogHHI_j';Zerovec]]);
b_dlogw_dlogHHI     = b_dlogw_dlogHHI(2);

b_dlogn_dlogHHI     = regress([dlogn_j';Zerovec],[[ones_j';Onevec],[dlogHHI_j';Zerovec]]);
b_dlogn_dlogHHI     = b_dlogn_dlogHHI(2);

% Figure 8A
dlogHHI             = mean(dlogHHI_j);
% dlogHHI             = mean(dlogHHI_j(dlogHHI_j>prctile(dlogHHI_j,95)));

% Table 10 - Again for top ventile
wn_j0          = sum(n_ij0.*w_ij0);
n_j0           = sum(n_ij0);
w_j0           = wn_j0./n_j0;
wn_j1          = sum(nc_ij.*wc_ij);
n_j1           = sum(nc_ij);
w_j1           = wn_j1./n_j1;
dlogw_j        = log(w_j1)-log(w_j0);
p95dloghhi     = prctile(dlogHHI_j,95);
i_j            = (dlogHHI_j>p95dloghhi);
hhi_j_merge    = hhi_j(merge_j);

Nmerge         = numel(dlogw_j(i_j));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j)';Zerovec],[[ones(size(dlogw_j(i_j)))';Onevec],[dlogHHI_j(i_j)';Zerovec]]);
b_elast1       = b(2);
i_j2           = i_j & (hhi_j_merge>prctile(hhi_j_merge(i_j),50));
Nmerge         = numel(dlogw_j(i_j2));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j2)';Zerovec],[[ones(size(dlogw_j(i_j2)))';Onevec],[dlogHHI_j(i_j2)';Zerovec]]);
b_elast2       = b(2);
i_j3           = i_j & (hhi_j_merge<prctile(hhi_j_merge(i_j),50));
Nmerge         = numel(dlogw_j(i_j3));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j3)';Zerovec],[[ones(size(dlogw_j(i_j3)))';Onevec],[dlogHHI_j(i_j3)';Zerovec]]);
b_elast3       = b(2);

fprintf('All = %4.3f \t,\t Above median = %4.3f \t,\t Below median = %4.3f\n',[b_elast1,b_elast2,b_elast3]);
%______________________________________________________
% Figure 7
% wn_j0N          = sum(n_ij0.*w_ij0.*N_ij);
% n_j0N           = sum(n_ij0.*N_ij);
% w_j0N           = wn_j0N./n_j0N;
% wn_j1N          = sum(nc_ij.*wc_ij.*N_ij);
% n_j1N           = sum(nc_ij.*N_ij);
% w_j1N           = wn_j1N./n_j1N;
% dlogw_j         = log(w_j1N)-log(w_j0N);
% 
% dhhibins        = [0:0.05:0.7];
% Nbins           = numel(dhhibins);
% mean_dlogw_bin  = zeros(Nbins-1,1);
% for bb = (1:Nbins-1)
%     id_j                = (dlogHHI_j>=dhhibins(bb)) & (dlogHHI_j<=dhhibins(bb+1));
%     mean_dlogw_bin(bb)  = mean(dlogw_j(id_j & ~isnan(dlogw_j)));
% end
% p95 = prctile(dlogHHI_j,95);
% plot(dhhibins(2:end),mean_dlogw_bin,'bo-');grid on;hold on;
% plot([-10,10],[0,0],'k:');hold on;
% plot([p95,p95],[-0.3,0.2],'k-');hold on;
% xlim([0,0.7]);ylim([-0.3,0.2]);
%______________________________________________________

% lambda calculation
r           = (1/beta) - 1;
nc_j        = sum(nc_ij.^((eta+1)/eta)).^(eta/(eta+1));
Nc          = sum(nc_j.^((theta+1)/theta)).^(theta/(theta+1));
agg_y       = sum(sum(zbar.*z_ij0.*(n_ij0.^(alpha))))/(1-Delta);   
agg_yc      = sum(sum(zbar.*z_ij0.*(nc_ij.^(alpha))))/(1-Delta);
agg_cc      = agg_yc*(1-(delta/(r+delta))*KS); 
agg_c       = agg_y*(1-(delta/(r+delta))*KS);
lambda      = (agg_cc - ((1/varphibar).^(1/varphi)*(Nc^(1+1/varphi))/(1+1/varphi)) + ((1/varphibar).^(1/varphi)*(N0_index.^(1+1/varphi))/(1+1/varphi)) )/agg_c ;
lambda      = 100*(lambda-1);  

% change in concentration
wn_j0       = sum(w_ij0.*n_ij0);
wn_j1       = sum(wc_ij.*nc_ij);
omega_j0    = wn_j0/sum(wn_j0);
omega_j1    = wn_j1/sum(wn_j1);

hhi_0_w     = sum(omega_j0.*hhi_j_merge);
hhi_1_w     = sum(omega_j1.*hhi_j1);
hhi_0_uw    = mean(hhi_j_merge);
hhi_1_uw    = mean(hhi_j1);

dhhi_w      = hhi_1_w   - hhi_0_w;
dhhi_uw     = hhi_1_uw  - hhi_0_uw;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Welfare
param.sigma             = 0;
param.varphi            = glob.varphi;
X0                      = load('../results/matfiles/baseline_merger_none');
z_ij0                   = X0.z_ij;
mu_ij0                  = X0.mu_ij;
Mj                      = sum(z_ij0>0);
X                       = load('../results/matfiles/baseline_merger_top');
muc_ij                  = X.mu_ij;
zc_ij                   = X.z_ij;

% 1. Compute baseline economy
options.mu_omega_only   = 'N';
options.solve_scale     = 'Y';
options.PE              = 'N';
out                     = Welfare_closed_form_SM_arnold(mu_ij0,z_ij0,Mj,[],[],param,glob,options);
mu0                     = out.mu;
omega0                  = out.omega;
N0                      = out.N;
C0                      = out.C;
glob.zbar               = out.Ztilde;
glob.varphibar          = out.varphibar;

% 2. Compute merger economy
options.mu_omega_only   = 'N';
options.solve_scale     = 'N';
options.PE              = 'N';
out                     = Welfare_closed_form_SM_arnold(muc_ij,zc_ij,Mj,[],[],param,glob,options);
mum                     = out.mu;
omegam                  = out.omega;
Nm                      = out.N;
Cm                      = out.C;

% Changes in mu and omega
dlogmu                  = 100*(log(mum)    - log(mu0));
dlogomega               = 100*(log(omegam) - log(omega0));

% Fraction due to mu and omega
options.mu_omega_only   = 'N';
options.solve_scale     = 'N';
options.PE              = 'N';
out                     = Welfare_closed_form_SM_arnold(mu_ij0,z_ij0,Mj,mum,omegam,param,glob,options);
lambdaH                 = out.lambdaE;
frac_mu                 = out.lambdaE_mu/out.lambdaE;
frac_omega              = out.lambdaE_omega/out.lambdaE;

% Welfare checks
xxx                     = 1+1/varphi;
lambdafun               = @(C0,C1,N0,N1) (1/C0)*(C1 - C0 - varphibar^(-1/varphi)*(N1^xxx - N0^xxx)/xxx);
lambdaE                 = lambdafun(C0,Cm,N0,Nm);
lambdaZ                 = lambdafun(X0.C,X.C,X0.N,X.N); 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

results_top        = [ factor*mean_dlogn_Y_w,    factor*mean_dlogn_Y_uw,     factor*mean_dlogwn_Y_w ,   ... % Table 3
                        factor*mean_dlogn_Y_w_small , factor*mean_dlogn_Y_w_big  ,          ... % Table 4
                        factor*mean_dlogw_Y    ,                                           ... % Table 5
                        factor*mean_dlogw_HI     , factor*mean_dlogw_MI,  factor*mean_dlogw_LI,                   ... % Table 6 
                        b_dlogw_dlogHHI     , b_dlogn_dlogHHI,                      ... % Table 7
                        b_hhi ,                                                     ... % Table 8
                        dlogHHI, ...
                        mean_n_pre_Y , median_n_pre_Y , mean_w_pre_Y/1000,               ... % Summary stats
                        lambda , dhhi_w , dhhi_uw, ...
                        100*lambdaH, frac_mu, frac_omega, dlogmu, dlogomega, ...
                        0,0,0,0, ...
                       ];
% return

%% MERGER - BOTTOM TWO
fprintf('Bottom 2\n');
n_ij0           = X0.n_ij;
w_ij0           = X0.w_ij/dcpi;
z_ij            = X0.z_ij;
mu_ij0          = X0.mu_ij;
X               = load('../results/matfiles/baseline_merger_bottom');
nc_ij           = X.n_ij;
wc_ij           = X.w_ij;
muc_ij          = X.mu_ij;

merge_ij        = X.merge_ij;
merge_j         = sum(merge_ij)>0;
n_ij0           = n_ij0(:,merge_j);
w_ij0           = w_ij0(:,merge_j);
mu_ij0          = mu_ij0(:,merge_j);
z_ij0           = z_ij(:,merge_j);
nc_ij           = nc_ij(:,merge_j);
wc_ij           = wc_ij(:,merge_j);
muc_ij          = muc_ij(:,merge_j);
zc_ij           = X0.z_ij(:,merge_j);
merge_ij        = merge_ij(:,merge_j);
swn_ij          = X.swn_ij(:,merge_j);

% return

Nc_bodies       = sum(sum(nc_ij));
Wc_mean         = sum(sum(wc_ij.*nc_ij))/Nc_bodies;

% Compute new HHIs
Y_ij            = logical(merge_ij) & (n_ij0>0);     % YES - merge
N_ij            = logical(~merge_ij) & (n_ij0>0);    % NO  - not merge
swn_ij_Y        = sum(Y_ij.*swn_ij);    % Combined share of merged firms
swn_ij_N        = N_ij.*swn_ij;         % Shares of non-merged firms
swn_ij_new      = [swn_ij_Y;swn_ij_N];  % Concatenate
hhi_j1          = sum(swn_ij_new.^2);

% Table 8 - Compute counterfactual initial HHIs with just firm removal
swn_ij0             = X0.swn_ij(:,merge_j);
swn_ij0_Y           = sum(Y_ij.*swn_ij0);                               % Combined share of merged firms
swn_ij0_N           = N_ij.*swn_ij0;                                    % Shares of non-merged firms
swn_ij0_new         = [swn_ij0_Y;swn_ij0_N];                            % Concatenate
swn_ij0_new         = bsxfun(@rdivide,swn_ij0_new,sum(swn_ij0_new));    % Renormalize
hhi_j_new           = sum(swn_ij0_new.^2);
dhhi_j_hat          = hhi_j_new - hhi_j0(merge_j);
dhhi_j_actual       = hhi_j1    - hhi_j0(merge_j);

Nmerge              = numel(dhhi_j_actual);
Ntotal              = floor(Nmerge/factor);
Zerovec             = zeros(Ntotal-Nmerge,1);
Onevec              = ones(Ntotal-Nmerge,1);
b_hhi               = regress([dhhi_j_actual';Zerovec],[[ones(size(dhhi_j_hat))';Onevec],[dhhi_j_hat';Zerovec]]);
b_hhi               = b_hhi(2);

n_ij0_Y_vec     = reshape(n_ij0(Y_ij),sum(sum(Y_ij)),1);
n_ij0_N_vec     = reshape(n_ij0(N_ij),sum(sum(N_ij)),1);
w_ij0_Y_vec     = reshape(w_ij0(Y_ij),sum(sum(Y_ij)),1);
w_ij0_N_vec     = reshape(w_ij0(N_ij),sum(sum(N_ij)),1);

n_ij1_Y_vec     = reshape(nc_ij(Y_ij),sum(sum(Y_ij)),1);
n_ij1_N_vec     = reshape(nc_ij(N_ij),sum(sum(N_ij)),1);
w_ij1_Y_vec     = reshape(wc_ij(Y_ij),sum(sum(Y_ij)),1);
w_ij1_N_vec     = reshape(wc_ij(N_ij),sum(sum(N_ij)),1);

mean_n_pre_Y    = mean(n_ij0_Y_vec);
mean_n_pre_N    = mean(n_ij0_N_vec);
median_n_pre_Y  = prctile(n_ij0_Y_vec,50);
median_n_pre_N  = prctile(n_ij0_N_vec,50);
mean_w_pre_Y    = mean(w_ij0_Y_vec);
mean_w_pre_N    = mean(w_ij0_N_vec);

dlogn_ij_Y      = log(n_ij1_Y_vec)-log(n_ij0_Y_vec);
dlogw_ij_Y      = log(w_ij1_Y_vec)-log(w_ij0_Y_vec);
dlogwn_ij_Y     = log(w_ij1_Y_vec.*n_ij1_Y_vec) - log(w_ij0_Y_vec.*n_ij0_Y_vec);

omega_ij_Y      = n_ij0_Y_vec/sum(n_ij0_Y_vec);
ones_ij_Y       = ones(size(n_ij0_Y_vec))/sum(ones(size(n_ij0_Y_vec)));

% Table 3
mean_dlogn_Y_w  = sum(omega_ij_Y.*dlogn_ij_Y);
mean_dlogn_Y_uw = sum(ones_ij_Y.*dlogn_ij_Y);
mean_dlogwn_Y_w = sum(omega_ij_Y.*dlogwn_ij_Y);

% Table 4 - Weighted
i_big                   = (n_ij0_Y_vec>prctile(n_ij0_Y_vec,50));
% i_big                   = (n_ij0_Y_vec>prctile(n_ij0(n_ij0>0),50));
omega_ij_Y_small        = n_ij0_Y_vec(~i_big)/sum(n_ij0_Y_vec(~i_big));
mean_dlogn_Y_w_small    = sum(omega_ij_Y_small.*dlogn_ij_Y(~i_big));
omega_ij_Y_big          = n_ij0_Y_vec(i_big)/sum(n_ij0_Y_vec(i_big));
mean_dlogn_Y_w_big      = sum(omega_ij_Y_big.*dlogn_ij_Y(i_big));

% Table 5 - Earnings 
% Here we weight by the number of workers, since its estimated by Arnold i
% worker level data. Therefore, e.g. if there are 100 workers in one firm
% and 200 workers in another firm, he would have 300 observations.
% Weighting by employment blows up the observations at the bigger firms,
% which is exactly what worker level data does.
mean_dlogw_Y            = sum(omega_ij_Y.*dlogw_ij_Y);  % (-0.008)

% Table 6 - By Initial concentration
dhhi_j              = log(hhi_j1) - log(hhi_j0(merge_j));
dhhi_cut            = prctile(dhhi_j_hat,75);
LI_j                = (dhhi_j_hat < dhhi_cut);  % Below top quantile for dhhi
hhi_j0_merge        = hhi_j0(merge_j);
hhi_cut             = prctile(hhi_j0_merge(~LI_j),50);
MI_j                = (~LI_j) & (hhi_j0_merge<hhi_cut);     % Above top quantile for dhhi & below median hhi in these
HI_j                = (~LI_j) & (hhi_j0_merge>hhi_cut);     % Above top quantile for dhhi & above median hhi in these

% (a) High impact
n_ij1_conc          = nc_ij(:,HI_j);
w_ij0_conc          = w_ij0(:,HI_j);
w_ij1_conc          = wc_ij(:,HI_j );
Y_ij_conc           = Y_ij(:,HI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_HI    = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_HI       = sum(omega_ij_conc.*dlogw_ij_conc_HI);    % (-0.031)

% (b) Medium impact
n_ij1_conc          = nc_ij(:,MI_j);
w_ij0_conc          = w_ij0(:,MI_j);
w_ij1_conc          = wc_ij(:,MI_j );
Y_ij_conc           = Y_ij(:,MI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_MI       = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_MI       = sum(omega_ij_conc.*dlogw_ij_conc_MI);    % (-0.031)

% (c) Low impact
n_ij1_conc          = nc_ij(:,LI_j);
w_ij0_conc          = w_ij0(:,LI_j);
w_ij1_conc          = wc_ij(:,LI_j );
Y_ij_conc           = Y_ij(:,LI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_LI    = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_LI       = sum(omega_ij_conc.*dlogw_ij_conc_LI);    % (-0.031)

% Table 7 - Correlations
dlogn_j             = log(sum(nc_ij))-log(sum(n_ij0));
wbar_j0             = sum(n_ij0.*w_ij0)./sum(n_ij0);
wbar_j1             = sum(nc_ij.*wc_ij)./sum(nc_ij);
dlogw_j             = log(wbar_j1) - log(wbar_j0); 
dlogHHI_j           = log(hhi_j1) - log(hhi_j0(merge_j));

ones_j              = ones(size(dlogw_j));

Nmerge              = numel(dlogw_j);
Ntotal              = floor(Nmerge/factor);
Zerovec             = zeros(Ntotal-Nmerge,1);
Onevec              = ones(Ntotal-Nmerge,1);

% Unweighted
b_dlogw_dlogHHI     = regress([dlogw_j';Zerovec],[[ones_j';Onevec],[dlogHHI_j';Zerovec]]);
b_dlogw_dlogHHI     = b_dlogw_dlogHHI(2);

b_dlogn_dlogHHI     = regress([dlogn_j';Zerovec],[[ones_j';Onevec],[dlogHHI_j';Zerovec]]);
b_dlogn_dlogHHI     = b_dlogn_dlogHHI(2);

% Figure 8A
dlogHHI             = mean(dlogHHI_j);
% dlogHHI             = mean(dlogHHI_j(dlogHHI_j>prctile(dlogHHI_j,95)));

% Table 10 - Again for top ventile
wn_j0          = sum(n_ij0.*w_ij0);
n_j0           = sum(n_ij0);
w_j0           = wn_j0./n_j0;
wn_j1          = sum(nc_ij.*wc_ij);
n_j1           = sum(nc_ij);
w_j1           = wn_j1./n_j1;
dlogw_j        = log(w_j1)-log(w_j0);
p95dloghhi     = prctile(dlogHHI_j,95);
i_j            = (dlogHHI_j>p95dloghhi);
hhi_j_merge    = hhi_j(merge_j);

Nmerge         = numel(dlogw_j(i_j));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j)';Zerovec],[[ones(size(dlogw_j(i_j)))';Onevec],[dlogHHI_j(i_j)';Zerovec]]);
b_elast1       = b(2);
i_j2           = i_j & (hhi_j_merge>prctile(hhi_j_merge(i_j),50));
Nmerge         = numel(dlogw_j(i_j2));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j2)';Zerovec],[[ones(size(dlogw_j(i_j2)))';Onevec],[dlogHHI_j(i_j2)';Zerovec]]);
b_elast2       = b(2);
i_j3           = i_j & (hhi_j_merge<prctile(hhi_j_merge(i_j),50));
Nmerge         = numel(dlogw_j(i_j3));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j3)';Zerovec],[[ones(size(dlogw_j(i_j3)))';Onevec],[dlogHHI_j(i_j3)';Zerovec]]);
b_elast3       = b(2);

fprintf('All = %4.3f \t,\t Above median = %4.3f \t,\t Below median = %4.3f\n',[b_elast1,b_elast2,b_elast3]);
%______________________________________________________
% Figure 7
% wn_j0N          = sum(n_ij0.*w_ij0.*N_ij);
% n_j0N           = sum(n_ij0.*N_ij);
% w_j0N           = wn_j0N./n_j0N;
% wn_j1N          = sum(nc_ij.*wc_ij.*N_ij);
% n_j1N           = sum(nc_ij.*N_ij);
% w_j1N           = wn_j1N./n_j1N;
% dlogw_j         = log(w_j1N)-log(w_j0N);
% 
% dhhibins        = [0:0.05:0.7];
% Nbins           = numel(dhhibins);
% mean_dlogw_bin  = zeros(Nbins-1,1);
% for bb = (1:Nbins-1)
%     id_j                = (dlogHHI_j>=dhhibins(bb)) & (dlogHHI_j<=dhhibins(bb+1));
%     mean_dlogw_bin(bb)  = mean(dlogw_j(id_j & ~isnan(dlogw_j)));
% end
% p95 = prctile(dlogHHI_j,95);
% plot(dhhibins(2:end),mean_dlogw_bin,'bo-');grid on;hold on;
% plot([-10,10],[0,0],'k:');hold on;
% plot([p95,p95],[-0.3,0.2],'k-');hold on;
% xlim([0,0.7]);ylim([-0.3,0.2]);
%______________________________________________________

% lambda calculation
r           = (1/beta) - 1;
nc_j        = sum(nc_ij.^((eta+1)/eta)).^(eta/(eta+1));
Nc          = sum(nc_j.^((theta+1)/theta)).^(theta/(theta+1));
agg_y       = sum(sum(zbar.*z_ij0.*(n_ij0.^(alpha))))/(1-Delta);   
agg_yc      = sum(sum(zbar.*z_ij0.*(nc_ij.^(alpha))))/(1-Delta);
agg_cc      = agg_yc*(1-(delta/(r+delta))*KS); 
agg_c       = agg_y*(1-(delta/(r+delta))*KS);

xxx         = 1+1/varphi;
lambdafun   = @(C0,C1,N0,N1) (1/C0)*(C1 - C0 - varphibar^(-1/varphi)*(N1^xxx - N0^xxx)/xxx);
lambda      = 100*lambdafun(agg_c,agg_cc,N0_index,Nc);

% change in concentration
% wn_j0       = sum(w_ij0.*n_ij0);
% wn_j1       = sum(wc_ij.*nc_ij);
% omega_j0    = wn_j0/sum(wn_j0);
% omega_j1    = wn_j1/sum(wn_j1);
% 
% hhi_0_w     = sum(omega_j0.*hhi_j_merge);
% hhi_1_w     = sum(omega_j1.*hhi_j1);
% hhi_0_uw    = mean(hhi_j_merge);
% hhi_1_uw    = mean(hhi_j1);
% 
% dhhi_w      = hhi_1_w   - hhi_0_w;
% dhhi_uw     = hhi_1_uw  - hhi_0_uw;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Welfare
param.sigma             = 0;
param.varphi            = glob.varphi;
X0                      = load('../results/matfiles/baseline_merger_none');
z_ij0                   = X0.z_ij;
mu_ij0                  = X0.mu_ij;
Mj                      = sum(z_ij0>0);
X                       = load('../results/matfiles/baseline_merger_bottom');
muc_ij                  = X.mu_ij;
zc_ij                   = X.z_ij;

% 1. Compute baseline economy
options.mu_omega_only   = 'N';
options.solve_scale     = 'Y';
options.PE              = 'N';
out                     = Welfare_closed_form_SM_arnold(mu_ij0,z_ij0,Mj,[],[],param,glob,options);
mu0                     = out.mu;
omega0                  = out.omega;
N0                      = out.N;
C0                      = out.C;
Y0                      = out.Y;
glob.zbar               = out.Ztilde;
glob.varphibar          = out.varphibar;

% 2. Compute merger economy
options.mu_omega_only   = 'N';
options.solve_scale     = 'N';
options.PE              = 'N';
out                     = Welfare_closed_form_SM_arnold(muc_ij,zc_ij,Mj,[],[],param,glob,options);
mum                     = out.mu;
omegam                  = out.omega;
Nm                      = out.N;
Cm                      = out.C;
Ym                      = out.Y;

% N0/N0_index
% Nm/Nc
% C0/agg_c
% Cm/agg_cc

% Changes in mu and omega
dlogmu                  = 100*(log(mum)    - log(mu0));
dlogomega               = 100*(log(omegam) - log(omega0));

% Fraction due to mu and omega
options.mu_omega_only   = 'N';
options.solve_scale     = 'N';
options.PE              = 'N';
out                     = Welfare_closed_form_SM_arnold(mu_ij0,z_ij0,Mj,mum,omegam,param,glob,options);
lambdaH                 = out.lambdaE;
frac_mu                 = out.lambdaE_mu/out.lambdaE;
frac_omega              = out.lambdaE_omega/out.lambdaE;

% Welfare checks
% xxx                     = 1+1/varphi;
% lambdafun               = @(C0,C1,N0,N1) (1/C0)*(C1 - C0 - varphibar^(-1/varphi)*(N1^xxx - N0^xxx)/xxx);
lambdaE                 = lambdafun(C0,Cm,N0,Nm);
lambdaZ                 = lambdafun(X0.C,X.C,X0.N,X.N); 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

results_bottom        = [ factor*mean_dlogn_Y_w,    factor*mean_dlogn_Y_uw,     factor*mean_dlogwn_Y_w ,   ... % Table 3
                        factor*mean_dlogn_Y_w_small , factor*mean_dlogn_Y_w_big  ,          ... % Table 4
                        factor*mean_dlogw_Y    ,                                           ... % Table 5
                        factor*mean_dlogw_HI     , factor*mean_dlogw_MI,  factor*mean_dlogw_LI,                   ... % Table 6 
                        b_dlogw_dlogHHI     , b_dlogn_dlogHHI,                      ... % Table 7
                        b_hhi ,                                                     ... % Table 8
                        dlogHHI, ...
                        mean_n_pre_Y , median_n_pre_Y , mean_w_pre_Y/1000,               ... % Summary stats
                        lambda , dhhi_w , dhhi_uw, ...
                        100*lambdaH, frac_mu, frac_omega, dlogmu, dlogomega, ...
                        0,0,0,0, ...
                       ];
% return
                   
%% MERGER - RANDOM
fprintf('Random 1\n');
n_ij0           = X0.n_ij;
X               = load('../results/matfiles/baseline_merger_random');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Cut off
merge_ij        = X.merge_ij;
Y_ij            = logical(merge_ij) & (n_ij0>0);
n_ij0_merge     = n_ij0.*Y_ij;

nbig_j          = max(n_ij0_merge);
nsmall_j        = sum(n_ij0_merge) - nbig_j;

% - SELECT ON SMALLEST SIZE OF MERGING FIRMS
nbar_j          = sum(n_ij0_merge)/2;
% nbar_j          = nsmall_j;

ncut_max        = 200;
ncut            = [1:1:ncut_max];
for i=1:length(ncut)
    n_ij_temp       = n_ij0_merge(:,nbar_j>ncut(i));
    n_ij_temp       = n_ij_temp(n_ij_temp>0);
    
    mean_size(i)    = mean(n_ij_temp);
    median_size(i)  = median(n_ij_temp);
    nmarkets(i)     = sum(nbar_j>ncut(i));
end

% subplot(1,2,1);
% plot(ncut,mean_size);hold on;
% plot(ncut,median_size,'r-');hold on;
% subplot(1,2,2);
% plot(ncut,nmarkets);

median_size_target  = 116;   % Arnold Table 1 - Pseudo median
% median_size_target  = 30;   % Arnold Table 1 - Pseudo median
[~,icut]            = min(abs(median_size - median_size_target));
jstar               = (nbar_j>ncut(icut));
fprintf('ncut = %5i\n',ncut(icut));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_ij0           = X0.n_ij(:,jstar);
w_ij0           = X0.w_ij(:,jstar)/dcpi;
z_ij            = X0.z_ij(:,jstar);

merge_ij        = X.merge_ij(:,jstar);
swn_ij          = X.swn_ij(:,jstar);
nc_ij           = X.n_ij(:,jstar);
wc_ij           = X.w_ij(:,jstar);

Nc_bodies       = sum(sum(nc_ij));
Wc_mean         = sum(sum(wc_ij.*nc_ij))/Nc_bodies;

% Compute new HHIs
Y_ij            = logical(merge_ij) & (n_ij0>0);        % YES - merge
N_ij            = logical(~merge_ij) & (n_ij0>0);       % NO  - not merge
swn_ij_Y        = sum(Y_ij.*swn_ij);                    % Combined share of merged firms
swn_ij_N        = N_ij.*swn_ij;                         % Shares of non-merged firms
swn_ij_new      = [swn_ij_Y;swn_ij_N];                  % Concatenate
hhi_j1          = sum(swn_ij_new.^2);

% Table 8 - Compute counterfactual initial HHIs with just firm removal
swn_ij0             = X0.swn_ij(:,jstar);
swn_ij0_Y           = sum(Y_ij.*swn_ij0);                               % Combined share of merged firms
swn_ij0_N           = N_ij.*swn_ij0;                                    % Shares of non-merged firms
swn_ij0_new         = [swn_ij0_Y;swn_ij0_N];                            % Concatenate
swn_ij0_new         = bsxfun(@rdivide,swn_ij0_new,sum(swn_ij0_new));    % Renormalize
hhi_j_new           = sum(swn_ij0_new.^2);
dhhi_j_hat          = hhi_j_new - hhi_j0(jstar);        % Simple predicted dhhi
dhhi_j_actual       = hhi_j1    - hhi_j0(jstar);

Nmerge              = numel(dhhi_j_actual);
Ntotal              = floor(Nmerge/factor);
Zerovec             = zeros(Ntotal-Nmerge,1);
Onevec              = ones(Ntotal-Nmerge,1);
b_hhi               = regress([dhhi_j_actual';Zerovec],[[ones(size(dhhi_j_hat))';Onevec],[dhhi_j_hat';Zerovec]]);
b_hhi               = b_hhi(2);

n_ij0_Y_vec     = reshape(n_ij0(Y_ij),sum(sum(Y_ij)),1);
n_ij0_N_vec     = reshape(n_ij0(N_ij),sum(sum(N_ij)),1);
w_ij0_Y_vec     = reshape(w_ij0(Y_ij),sum(sum(Y_ij)),1);
w_ij0_N_vec     = reshape(w_ij0(N_ij),sum(sum(N_ij)),1);
wn_ij0_Y_vec    = reshape(w_ij0(Y_ij).*n_ij0(Y_ij),sum(sum(Y_ij)),1);
wn_ij0_N_vec    = reshape(w_ij0(N_ij).*n_ij0(N_ij),sum(sum(N_ij)),1);

n_ij1_Y_vec     = reshape(nc_ij(Y_ij),sum(sum(Y_ij)),1);
n_ij1_N_vec     = reshape(nc_ij(N_ij),sum(sum(N_ij)),1);
w_ij1_Y_vec     = reshape(wc_ij(Y_ij),sum(sum(Y_ij)),1);
w_ij1_N_vec     = reshape(wc_ij(N_ij),sum(sum(N_ij)),1);

mean_n_pre_Y    = mean(n_ij0_Y_vec);
mean_n_pre_N    = mean(n_ij0_N_vec);
median_n_pre_Y  = prctile(n_ij0_Y_vec,50);
median_n_pre_N  = prctile(n_ij0_N_vec,50);
mean_w_pre_Y    = sum(wn_ij0_Y_vec)/sum(n_ij0_Y_vec);
mean_w_pre_N    = sum(wn_ij0_N_vec)/sum(n_ij0_N_vec);

dlogn_ij_Y      = log(n_ij1_Y_vec)-log(n_ij0_Y_vec);
dlogw_ij_Y      = log(w_ij1_Y_vec)-log(w_ij0_Y_vec);
dlogwn_ij_Y     = log(w_ij1_Y_vec.*n_ij1_Y_vec) - log(w_ij0_Y_vec.*n_ij0_Y_vec);

omega_ij_Y      = n_ij0_Y_vec/sum(n_ij0_Y_vec);
ones_ij_Y       = ones(size(n_ij0_Y_vec))/sum(ones(size(n_ij0_Y_vec)));

% Table 3
mean_dlogn_Y_w          = sum(omega_ij_Y.*dlogn_ij_Y);
mean_dlogn_Y_uw         = sum(ones_ij_Y.*dlogn_ij_Y);
mean_dlogwn_Y_w         = sum(omega_ij_Y.*dlogwn_ij_Y);

% Table 4 - Weighted
i_big                   = (n_ij0_Y_vec>prctile(n_ij0_Y_vec,50));
% i_big                   = (n_ij0_Y_vec>prctile(n_ij0(n_ij0>0),50));
omega_ij_Y_small        = n_ij0_Y_vec(~i_big)/sum(n_ij0_Y_vec(~i_big));
mean_dlogn_Y_w_small    = sum(omega_ij_Y_small.*dlogn_ij_Y(~i_big));
omega_ij_Y_big          = n_ij0_Y_vec(i_big)/sum(n_ij0_Y_vec(i_big));
mean_dlogn_Y_w_big      = sum(omega_ij_Y_big.*dlogn_ij_Y(i_big));

% Table 5 - Earnings 
% Here we weight by the number of workers remaining (as he uses stayers), since its estimated by Arnold i
% worker level data. Therefore, e.g. if there are 100 workers in one firm
% and 200 workers in another firm, he would have 300 observations.
% Weighting by employment blows up the observations at the bigger firms,
% which is exactly what worker level data does.
omega_ij1_Y         = n_ij1_Y_vec/sum(n_ij1_Y_vec);
mean_dlogw_Y        = sum(omega_ij1_Y.*dlogw_ij_Y);  % (-0.008)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Table 6 - By Initial concentration
dhhi_cut            = prctile(dhhi_j_hat,75);
LI_j                = (dhhi_j_hat < dhhi_cut);  % Below top quantile for dhhi
hhi_j0_merge        = hhi_j0(jstar);
hhi_cut             = prctile(hhi_j0_merge(~LI_j),50);
MI_j                = (~LI_j) & (hhi_j0_merge<hhi_cut);     % Above top quantile for dhhi & below median hhi in these
HI_j                = (~LI_j) & (hhi_j0_merge>hhi_cut);     % Above top quantile for dhhi & above median hhi in these

% (a) High impact
n_ij1_conc          = nc_ij(:,HI_j);
w_ij0_conc          = w_ij0(:,HI_j);
w_ij1_conc          = wc_ij(:,HI_j );
Y_ij_conc           = Y_ij(:,HI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_HI    = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_HI       = sum(omega_ij_conc.*dlogw_ij_conc_HI);  

% (b) Medium impact
n_ij1_conc          = nc_ij(:,MI_j);
w_ij0_conc          = w_ij0(:,MI_j);
w_ij1_conc          = wc_ij(:,MI_j );
Y_ij_conc           = Y_ij(:,MI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_MI       = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_MI       = sum(omega_ij_conc.*dlogw_ij_conc_MI);    

% (c) Low impact
n_ij1_conc          = nc_ij(:,LI_j);
w_ij0_conc          = w_ij0(:,LI_j);
w_ij1_conc          = wc_ij(:,LI_j );
Y_ij_conc           = Y_ij(:,LI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_LI    = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_LI       = sum(omega_ij_conc.*dlogw_ij_conc_LI);    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Table 7 - Correlations
dlogn_j             = log(sum(nc_ij))-log(sum(n_ij0));
wbar_j0             = sum(n_ij0.*w_ij0)./sum(n_ij0);
wbar_j1             = sum(nc_ij.*wc_ij)./sum(nc_ij);
dlogw_j             = log(wbar_j1) - log(wbar_j0); 
dlogHHI_j           = log(hhi_j1) - log(hhi_j0(jstar));

ones_j              = ones(size(dlogw_j));

Nmerge              = numel(dlogw_j);
Ntotal              = floor(Nmerge/factor);
Zerovec             = zeros(Ntotal-Nmerge,1);
Onevec              = ones(Ntotal-Nmerge,1);

b_dlogw_dlogHHI     = regress([dlogw_j';Zerovec],[[ones_j';Onevec],[dlogHHI_j';Zerovec]]);
b_dlogw_dlogHHI     = b_dlogw_dlogHHI(2);

b_dlogn_dlogHHI     = regress([dlogn_j';Zerovec],[[ones_j';Onevec],[dlogHHI_j';Zerovec]]);
b_dlogn_dlogHHI     = b_dlogn_dlogHHI(2);

% Figure 8A
dlogHHI             = mean(dlogHHI_j);
% dlogHHI             = mean(dlogHHI_j(dlogHHI_j>prctile(dlogHHI_j,95)));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Table 10 - Again for top ventile
%__________________________________________________________________________
% Taking his heading for Table 10 literally:
logw_ij0            = log(w_ij0);
logw_ij0(n_ij0==0)  = 0;
logwc_ij            = log(wc_ij);
logwc_ij(nc_ij==0)  = 0;

omega_ij0      = n_ij0./sum(n_ij0);
meanlogw_j0    = sum(omega_ij0.*logw_ij0);
n_j0           = sum(n_ij0);
omega_ij1      = nc_ij./sum(nc_ij);
meanlogw_j1    = sum(omega_ij1.*logwc_ij);
n_j1           = sum(nc_ij);
dlogw_j        = meanlogw_j1-meanlogw_j0;       % Change in average log worker earnings
dlogn_j        = log(n_j1)-log(n_j0);           % Log change in market employment
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This matches the dlogHHI from Figure 8A
dlogHHI_cut     = prctile(dlogHHI_j,16); 
fprintf('Mean dlogHHI_j =  %4.3f\n',mean(dlogHHI_j(dlogHHI_j>dlogHHI_cut)));
i_j             = (dlogHHI_j>dlogHHI_cut);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hhi_j_merge    = hhi_j(jstar);
HHI_cut        = prctile(hhi_j_merge(i_j),50); 
% Above median of selected markets *after* applying dlogHHI_cut
%__________________________________________________________________________
% Log C
i_j1           = i_j;
Nmerge         = numel(dlogw_j(i_j1));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j1)';Zerovec],[[ones(size(dlogw_j(i_j1)))';Onevec],[dlogHHI_j(i_j1)';Zerovec]]);
b_elast1_w     = b(2);
b              = regress([dlogn_j(i_j1)';Zerovec],[[ones(size(dlogn_j(i_j1)))';Onevec],[dlogHHI_j(i_j1)';Zerovec]]);
b_elast1_n     = b(2);
%__________________________________________________________________________
% Log C x Above Median C
i_j2           = i_j & (hhi_j_merge>HHI_cut);
Nmerge         = numel(dlogw_j(i_j2));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j2)';Zerovec],[[ones(size(dlogw_j(i_j2)))';Onevec],[dlogHHI_j(i_j2)';Zerovec]]);
b_elast2_w     = b(2); 
b              = regress([dlogn_j(i_j2)';Zerovec],[[ones(size(dlogn_j(i_j2)))';Onevec],[dlogHHI_j(i_j2)';Zerovec]]);
b_elast2_n       = b(2);
%__________________________________________________________________________
% Log C x Below Median C
i_j3           = i_j & (hhi_j_merge<HHI_cut);
Nmerge         = numel(dlogw_j(i_j3));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j3)';Zerovec],[[ones(size(dlogw_j(i_j3)))';Onevec],[dlogHHI_j(i_j3)';Zerovec]]);
b_elast3_w       = b(2);
b              = regress([dlogn_j(i_j3)';Zerovec],[[ones(size(dlogn_j(i_j3)))';Onevec],[dlogHHI_j(i_j3)';Zerovec]]);
b_elast3_n       = b(2);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('This is where I am playing with the Table 10 results\n');
fprintf('1. Wages\n');
fprintf('BHM: All = %4.3f \t,\t Above median = %4.3f \t,\t Below median = %4.3f\n',[b_elast1_w,b_elast2_w,b_elast3_w]);
fprintf('DA:  All = -0.147 \t,\t Above median = -0.176 \t,\t Below median = 0.058\n');
fprintf('1. Employment\n');
fprintf('BHM: All = %4.3f \t,\t Above median = %4.3f \t,\t Below median = %4.3f\n',[b_elast1_n,b_elast2_n,b_elast3_n]);
fprintf('DA:  All = -0.450 \t,\t Above median = N/A \t,\t Below median = N/A\n');
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
%______________________________________________________
% % Figure 7
% wn_j0N          = sum(n_ij0.*w_ij0.*N_ij);
% n_j0N           = sum(n_ij0.*N_ij);
% w_j0N           = wn_j0N./n_j0N;
% wn_j1N          = sum(nc_ij.*wc_ij.*N_ij);
% n_j1N           = sum(nc_ij.*N_ij);
% w_j1N           = wn_j1N./n_j1N;
% dlogw_j         = log(w_j1N)-log(w_j0N);
% 
% dhhibins        = [0:0.05:0.50];
% Nbins           = numel(dhhibins);
% mean_dlogw_bin  = zeros(Nbins-1,1);
% p5_dlogw_bin    = zeros(Nbins-1,1);
% p95_dlogw_bin   = zeros(Nbins-1,1);
% for bb = (1:Nbins-1)
%     id_j                = (dlogHHI_j>=dhhibins(bb)) & (dlogHHI_j<=dhhibins(bb+1));
%     mean_dlogw_bin(bb)  = mean(dlogw_j(id_j & ~isnan(dlogw_j)));
%     p5_dlogw_bin(bb)    = prctile(dlogw_j(id_j & ~isnan(dlogw_j)),95);
%     p95_dlogw_bin(bb)   = prctile(dlogw_j(id_j & ~isnan(dlogw_j)),5);
% end
% p95 = prctile(dlogHHI_j,95);
% plot(dhhibins(2:end),mean_dlogw_bin,'ro-');grid on;hold on;
% plot(dhhibins(2:end),p5_dlogw_bin,'r-.');grid on;hold on;
% plot(dhhibins(2:end),p95_dlogw_bin,'r-.');grid on;hold on;
% plot([-10,10],[0,0],'k:');hold on;
% plot([p95,p95],[-0.3,0.2],'r-');hold on;
% xlim([0,0.7]);ylim([-0.3,0.2]);
%______________________________________________________

% lambda calculation
r           = (1/beta) - 1;
nc_j        = sum(nc_ij.^((eta+1)/eta)).^(eta/(eta+1));
Nc          = sum(nc_j.^((theta+1)/theta)).^(theta/(theta+1));
agg_y       = sum(sum(zbar.*z_ij.*(n_ij0.^(alpha))))/(1-Delta); 
agg_yc      = sum(sum(zbar.*z_ij.*(nc_ij.^(alpha))))/(1-Delta);
agg_cc      = agg_yc*(1-(delta/(r+delta))*KS); 
agg_c       = agg_y*(1-(delta/(r+delta))*KS);
lambda      = (agg_cc - ((1/varphibar).^(1/varphi)*(Nc^(1+1/varphi))/(1+1/varphi)) + ((1/varphibar).^(1/varphi)*(N0_index.^(1+1/varphi))/(1+1/varphi)) )/agg_c ;
lambda      = 100*(lambda-1);

% change in concentration
wn_j0       = sum(w_ij0.*n_ij0);
wn_j1       = sum(wc_ij.*nc_ij);
omega_j0    = wn_j0/sum(wn_j0);
omega_j1    = wn_j1/sum(wn_j1);

hhi_0_w     = sum(omega_j0.*hhi_j_merge);
hhi_1_w     = sum(omega_j1.*hhi_j1);
hhi_0_uw    = mean(hhi_j_merge);
hhi_1_uw    = mean(hhi_j1);

dhhi_w      = hhi_1_w   - hhi_0_w;
dhhi_uw     = hhi_1_uw  - hhi_0_uw;
%______________________________________________________

results_random   = [ factor*mean_dlogn_Y_w,    factor*mean_dlogn_Y_uw,     factor*mean_dlogwn_Y_w ,   ... % Table 3
                        factor*mean_dlogn_Y_w_small , factor*mean_dlogn_Y_w_big  ,          ... % Table 4
                        factor*mean_dlogw_Y    ,                                           ... % Table 5
                        factor*mean_dlogw_HI     , factor*mean_dlogw_MI,  factor*mean_dlogw_LI,                   ... % Table 6 
                        b_dlogw_dlogHHI     , b_dlogn_dlogHHI,                      ... % Table 7
                        b_hhi ,                                                     ... % Table 8
                        mean(dlogHHI_j(dlogHHI_j>dlogHHI_cut)), ...
                        mean_n_pre_Y , median_n_pre_Y , mean_w_pre_Y/1000,               ... % Summary stats
                        lambda, dhhi_w, dhhi_uw, ...
                        0,0,0,0,0,...
                        b_elast1_w,b_elast2_w,b_elast3_w,b_elast1_n, ...
                       ];

%% MERGER - RANDOM 2
fprintf('Random 2\n');
n_ij0           = X0.n_ij;
w_ij0           = X0.w_ij/dcpi;
X               = load('../results/matfiles/baseline_merger_random');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Cut off
merge_ij        = X.merge_ij;
Y_ij            = logical(merge_ij) & (n_ij0>0);
n_ij0_merge     = n_ij0.*Y_ij;

nbig_j          = max(n_ij0_merge);
nsmall_j        = sum(n_ij0_merge) - nbig_j;

% - SELECT ON AVERAGE SIZE OF MERGING FIRMS
nbar_j          = sum(n_ij0_merge)/2;   
% nbar_j          = nsmall_j;

ncut_max        = 200;
ncut            = [1:1:ncut_max];
for i=1:length(ncut)
    n_ij_temp       = n_ij0_merge(:,nbar_j>ncut(i));
    n_ij_temp       = n_ij_temp(n_ij_temp>0);
    
    mean_size(i)    = mean(n_ij_temp);
    median_size(i)  = median(n_ij_temp);
    nmarkets(i)     = sum(nbar_j>ncut(i));
end

% subplot(1,2,1);
% plot(ncut,mean_size);hold on;
% plot(ncut,median_size,'r-');hold on;
% subplot(1,2,2);
% plot(ncut,nmarkets);

median_size_target  = 116;   % Arnold Table 1 - Pseudo median
% median_size_target  = 30;   % Arnold Table 1 - Pseudo median
[~,icut]            = min(abs(median_size - median_size_target));
jstar               = (nbar_j>ncut(icut));
fprintf('ncut = %5i\n',ncut(icut));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_ij0           = X0.n_ij(:,jstar);
w_ij0           = X0.w_ij(:,jstar)/dcpi;
z_ij            = X0.z_ij(:,jstar);

merge_ij        = X.merge_ij(:,jstar);
swn_ij          = X.swn_ij(:,jstar);
nc_ij           = X.n_ij(:,jstar);
wc_ij           = X.w_ij(:,jstar);

Nc_bodies       = sum(sum(nc_ij));
Wc_mean         = sum(sum(wc_ij.*nc_ij))/Nc_bodies;

% Compute new HHIs
Y_ij            = logical(merge_ij) & (n_ij0>0);     % YES - merge
N_ij            = logical(~merge_ij) & (n_ij0>0);    % NO  - not merge
swn_ij_Y        = sum(Y_ij.*swn_ij);    % Combined share of merged firms
swn_ij_N        = N_ij.*swn_ij;         % Shares of non-merged firms
swn_ij_new      = [swn_ij_Y;swn_ij_N];  % Concatenate
hhi_j1          = sum(swn_ij_new.^2);

% Table 8 - Compute counterfactual initial HHIs with just firm removal
swn_ij0             = X0.swn_ij(:,jstar);
swn_ij0_Y           = sum(Y_ij.*swn_ij0);                               % Combined share of merged firms
swn_ij0_N           = N_ij.*swn_ij0;                                    % Shares of non-merged firms
swn_ij0_new         = [swn_ij0_Y;swn_ij0_N];                            % Concatenate
swn_ij0_new         = bsxfun(@rdivide,swn_ij0_new,sum(swn_ij0_new));    % Renormalize
hhi_j_new           = sum(swn_ij0_new.^2);
dhhi_j_hat          = hhi_j_new - hhi_j0(jstar);
dhhi_j_actual       = hhi_j1    - hhi_j0(jstar);

Nmerge              = numel(dhhi_j_actual);
Ntotal              = floor(Nmerge/factor);
Zerovec             = zeros(Ntotal-Nmerge,1);
Onevec              = ones(Ntotal-Nmerge,1);
b_hhi               = regress([dhhi_j_actual';Zerovec],[[ones(size(dhhi_j_hat))';Onevec],[dhhi_j_hat';Zerovec]]);
b_hhi               = b_hhi(2);

n_ij0_Y_vec     = reshape(n_ij0(Y_ij),sum(sum(Y_ij)),1);
n_ij0_N_vec     = reshape(n_ij0(N_ij),sum(sum(N_ij)),1);
w_ij0_Y_vec     = reshape(w_ij0(Y_ij),sum(sum(Y_ij)),1);
w_ij0_N_vec     = reshape(w_ij0(N_ij),sum(sum(N_ij)),1);

n_ij1_Y_vec     = reshape(nc_ij(Y_ij),sum(sum(Y_ij)),1);
n_ij1_N_vec     = reshape(nc_ij(N_ij),sum(sum(N_ij)),1);
w_ij1_Y_vec     = reshape(wc_ij(Y_ij),sum(sum(Y_ij)),1);
w_ij1_N_vec     = reshape(wc_ij(N_ij),sum(sum(N_ij)),1);

mean_n_pre_Y    = mean(n_ij0_Y_vec);
mean_n_pre_N    = mean(n_ij0_N_vec);
median_n_pre_Y  = prctile(n_ij0_Y_vec,50);
median_n_pre_N  = prctile(n_ij0_N_vec,50);
mean_w_pre_Y    = mean(w_ij0_Y_vec);
mean_w_pre_N    = mean(w_ij0_N_vec);

dlogn_ij_Y      = log(n_ij1_Y_vec)-log(n_ij0_Y_vec);
dlogw_ij_Y      = log(w_ij1_Y_vec)-log(w_ij0_Y_vec);
dlogwn_ij_Y     = log(w_ij1_Y_vec.*n_ij1_Y_vec) - log(w_ij0_Y_vec.*n_ij0_Y_vec);

omega_ij_Y      = n_ij0_Y_vec/sum(n_ij0_Y_vec);
ones_ij_Y       = ones(size(n_ij0_Y_vec))/sum(ones(size(n_ij0_Y_vec)));

% Table 3
mean_dlogn_Y_w  = sum(omega_ij_Y.*dlogn_ij_Y);
mean_dlogn_Y_uw = sum(ones_ij_Y.*dlogn_ij_Y);
mean_dlogwn_Y_w = sum(omega_ij_Y.*dlogwn_ij_Y);

% Table 4 - Weighted
i_big                   = (n_ij0_Y_vec>prctile(n_ij0_Y_vec,50));
% i_big                   = (n_ij0_Y_vec>prctile(n_ij0(n_ij0>0),50));
omega_ij_Y_small        = n_ij0_Y_vec(~i_big)/sum(n_ij0_Y_vec(~i_big));
mean_dlogn_Y_w_small    = sum(omega_ij_Y_small.*dlogn_ij_Y(~i_big));
omega_ij_Y_big          = n_ij0_Y_vec(i_big)/sum(n_ij0_Y_vec(i_big));
mean_dlogn_Y_w_big      = sum(omega_ij_Y_big.*dlogn_ij_Y(i_big));

% Table 5 - Earnings 
% Here we weight by the number of workers remaining (as he uses stayers), since its estimated by Arnold i
% worker level data. Therefore, e.g. if there are 100 workers in one firm
% and 200 workers in another firm, he would have 300 observations.
% Weighting by employment blows up the observations at the bigger firms,
% which is exactly what worker level data does.
omega_ij1_Y         = n_ij1_Y_vec/sum(n_ij1_Y_vec);
mean_dlogw_Y        = sum(omega_ij1_Y.*dlogw_ij_Y);  % (-0.008)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Table 6 - By Initial concentration
dhhi_j              = log(hhi_j1) - log(hhi_j0(jstar));
dhhi_cut            = prctile(dhhi_j_hat,75);
LI_j                = (dhhi_j_hat < dhhi_cut);  % Below top quantile for dhhi
hhi_j0_merge        = hhi_j0(jstar);
hhi_cut             = prctile(hhi_j0_merge(~LI_j),50);
MI_j                = (~LI_j) & (hhi_j0_merge<hhi_cut);     % Above top quantile for dhhi & below median hhi in these
HI_j                = (~LI_j) & (hhi_j0_merge>hhi_cut);     % Above top quantile for dhhi & above median hhi in these

% (a) High impact
n_ij1_conc          = nc_ij(:,HI_j);
w_ij0_conc          = w_ij0(:,HI_j);
w_ij1_conc          = wc_ij(:,HI_j );
Y_ij_conc           = Y_ij(:,HI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_HI    = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_HI       = sum(omega_ij_conc.*dlogw_ij_conc_HI);    % (-0.031)

% (b) Medium impact
n_ij1_conc          = nc_ij(:,MI_j);
w_ij0_conc          = w_ij0(:,MI_j);
w_ij1_conc          = wc_ij(:,MI_j );
Y_ij_conc           = Y_ij(:,MI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_MI       = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_MI       = sum(omega_ij_conc.*dlogw_ij_conc_MI);    % (-0.031)

% (c) Low impact
n_ij1_conc          = nc_ij(:,LI_j);
w_ij0_conc          = w_ij0(:,LI_j);
w_ij1_conc          = wc_ij(:,LI_j );
Y_ij_conc           = Y_ij(:,LI_j);
n_ij1_conc_vec      = reshape(n_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij0_conc_vec      = reshape(w_ij0_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
w_ij1_conc_vec      = reshape(w_ij1_conc(Y_ij_conc),sum(sum(Y_ij_conc)),1);
omega_ij_conc       = n_ij1_conc_vec/sum(n_ij1_conc_vec);
dlogw_ij_conc_LI    = log(w_ij1_conc_vec) - log(w_ij0_conc_vec);
mean_dlogw_LI       = sum(omega_ij_conc.*dlogw_ij_conc_LI);    % (-0.031)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Table 7 - Correlations
dlogn_j             = log(sum(nc_ij))-log(sum(n_ij0));
wbar_j0             = sum(n_ij0.*w_ij0)./sum(n_ij0);
wbar_j1             = sum(nc_ij.*wc_ij)./sum(nc_ij);
dlogw_j             = log(wbar_j1) - log(wbar_j0); 
dlogHHI_j           = log(hhi_j1) - log(hhi_j0(jstar));

ones_j              = ones(size(dlogw_j));

Nmerge              = numel(dlogw_j);
Ntotal              = floor(Nmerge/factor);
Zerovec             = zeros(Ntotal-Nmerge,1);
Onevec              = ones(Ntotal-Nmerge,1);

b_dlogw_dlogHHI     = regress([dlogw_j';Zerovec],[[ones_j';Onevec],[dlogHHI_j';Zerovec]]);
b_dlogw_dlogHHI     = b_dlogw_dlogHHI(2);

b_dlogn_dlogHHI     = regress([dlogn_j';Zerovec],[[ones_j';Onevec],[dlogHHI_j';Zerovec]]);
b_dlogn_dlogHHI     = b_dlogn_dlogHHI(2);

% Figure 8A
dlogHHI             = mean(dlogHHI_j);
% dlogHHI             = mean(dlogHHI_j(dlogHHI_j>prctile(dlogHHI_j,95)));

arnoldFig7_plot = 'N';
if strcmp(arnoldFig7_plot,'Y');
    %______________________________________________________
    % Figure 7
    % (1) Average wage paid to workers
    wn_j0N          = sum(n_ij0.*w_ij0.*N_ij);
    n_j0N           = sum(n_ij0.*N_ij);
    w_j0N           = wn_j0N./n_j0N;
    wn_j1N          = sum(nc_ij.*wc_ij.*N_ij);
    n_j1N           = sum(nc_ij.*N_ij);
    w_j1N           = wn_j1N./n_j1N;
    dlogw_jN        = log(w_j1N)-log(w_j0N);

    wn_j0Y          = sum(n_ij0.*w_ij0.*Y_ij);
    n_j0Y           = sum(n_ij0.*Y_ij);
    w_j0Y           = wn_j0Y./n_j0Y;
    wn_j1Y          = sum(nc_ij.*wc_ij.*Y_ij);
    n_j1Y           = sum(nc_ij.*Y_ij);
    w_j1Y           = wn_j1Y./n_j1Y;
    dlogw_jY        = log(w_j1Y)-log(w_j0Y);

    % (2) Average wage paid by firms
    wn_j0N          = sum(w_ij0.*N_ij);
    n_j0N           = sum(N_ij);
    w_j0N           = wn_j0N./n_j0N;
    wn_j1N          = sum(wc_ij.*N_ij);
    n_j1N           = sum(N_ij);
    w_j1N           = wn_j1N./n_j1N;
    dlogw_jN2       = log(w_j1N)-log(w_j0N);

    wn_j0Y          = sum(w_ij0.*Y_ij);
    n_j0Y           = sum(Y_ij);
    w_j0Y           = wn_j0Y./n_j0Y;
    wn_j1Y          = sum(wc_ij.*Y_ij);
    n_j1Y           = sum(Y_ij);
    w_j1Y           = wn_j1Y./n_j1Y;
    dlogw_jY2       = log(w_j1Y)-log(w_j0Y);

    dhhibins        = [0:0.05:0.50];
    Nbins           = numel(dhhibins);
    % 1. Average wage paid to workers
    mean_dlogw_binN = zeros(Nbins-1,1);
    p5_dlogw_binN   = zeros(Nbins-1,1);
    p95_dlogw_binN  = zeros(Nbins-1,1);
    mean_dlogw_binY = zeros(Nbins-1,1);
    p5_dlogw_binY   = zeros(Nbins-1,1);
    p95_dlogw_binY  = zeros(Nbins-1,1);
    % 2. Average wage paid by firms
    mean_dlogw_binN2 = zeros(Nbins-1,1);
    p5_dlogw_binN2   = zeros(Nbins-1,1);
    p95_dlogw_binN2  = zeros(Nbins-1,1);
    mean_dlogw_binY2 = zeros(Nbins-1,1);
    p5_dlogw_binY2   = zeros(Nbins-1,1);
    p95_dlogw_binY2  = zeros(Nbins-1,1);
    for bb = (1:Nbins-1)
        id_j                 = (dlogHHI_j>=dhhibins(bb)) & (dlogHHI_j<=dhhibins(bb+1));
        % 1. Average wage paid to workers
        % Non merging firms
        mean_dlogw_binN(bb)  = mean(dlogw_jN(id_j & ~isnan(dlogw_jN)));
        p5_dlogw_binN(bb)    = prctile(dlogw_jN(id_j & ~isnan(dlogw_jN)),95);
        p95_dlogw_binN(bb)   = prctile(dlogw_jN(id_j & ~isnan(dlogw_jN)),5);
        % Merging firms
        mean_dlogw_binY(bb)  = mean(dlogw_jY(id_j & ~isnan(dlogw_jY)));
        p5_dlogw_binY(bb)    = prctile(dlogw_jY(id_j & ~isnan(dlogw_jY)),95);
        p95_dlogw_binY(bb)   = prctile(dlogw_jY(id_j & ~isnan(dlogw_jY)),5);
        % 2. Average wage paid by firms
        % Non merging firms
        mean_dlogw_binN2(bb)  = mean(dlogw_jN2(id_j & ~isnan(dlogw_jN2)));
        p5_dlogw_binN2(bb)    = prctile(dlogw_jN2(id_j & ~isnan(dlogw_jN2)),95);
        p95_dlogw_binN2(bb)   = prctile(dlogw_jN2(id_j & ~isnan(dlogw_jN2)),5);
        % Merging firms
        mean_dlogw_binY2(bb)  = mean(dlogw_jY2(id_j & ~isnan(dlogw_jY2)));
        p5_dlogw_binY2(bb)    = prctile(dlogw_jY2(id_j & ~isnan(dlogw_jY2)),95);
        p95_dlogw_binY2(bb)   = prctile(dlogw_jY2(id_j & ~isnan(dlogw_jY2)),5);
    end
    p95 = prctile(dlogHHI_j,95);
    p50 = prctile(dlogHHI_j,50);
    arnoldvec = -[0.01,0.01,0.005,linspace(0.0075,0.10,12)];
    arnoldbin = [0.025,0.050,0.10:0.05:0.70];

    F1 = figure(1);
    set(F1,'Pos',[170 558 1070 420]);

    subplot(1,2,1);
    plot(dhhibins(2:end),mean_dlogw_binN,'bo-','linewidth',2);grid on;hold on;
    plot(dhhibins(2:end),mean_dlogw_binY,'rx-','linewidth',2);grid on;hold on;
    plot(arnoldbin,arnoldvec,'gs-','linewidth',2);hold on;

    % plot(dhhibins(2:end),p5_dlogw_binN,'b-.');grid on;hold on;
    % plot(dhhibins(2:end),p95_dlogw_binN,'b-.');grid on;hold on;

    plot(dhhibins(2:end),mean_dlogw_binY,'rx-');grid on;hold on;
    % plot(dhhibins(2:end),p5_dlogw_binY,'r-.');grid on;hold on;
    % plot(dhhibins(2:end),p95_dlogw_binY,'r-.');grid on;hold on;

    plot([-10,10],[0,0],'k-');hold on;
    plot([p50,p50],[-0.3,0.2],'b-');hold on;
    xlim([0,0.5]);ylim([-0.3,0.1]);
    l = legend('Non-merging firms','Merging firms','Arnold for non-merging firms');
    set(l,'interpreter','latex','location','southwest');

    xlabel('Change in Log Local Labor Market Concentration','interpreter','latex');
    ylabel('Change in log average earnings','interpreter','latex');
    ytickformat('%3.2f');
    xtickformat('%3.1f');
    title('A. Change in average worker wages','interpreter','latex');
    set(gca,'ticklabelinterpreter','latex','fontsize',12);
    print '-dpng' Figure_Arnold_Figure7

    subplot(1,2,2);
    plot(dhhibins(2:end),mean_dlogw_binN2,'bo-','linewidth',2);grid on;hold on;
    plot(dhhibins(2:end),mean_dlogw_binY2,'rx-','linewidth',2);grid on;hold on;
    % plot(arnoldbin,arnoldvec,'gs-','linewidth',2);hold on;

    % plot(dhhibins(2:end),p5_dlogw_binN2,'b-.');grid on;hold on;
    % plot(dhhibins(2:end),p95_dlogw_binN2,'b-.');grid on;hold on;

    plot(dhhibins(2:end),mean_dlogw_binY2,'rx-');grid on;hold on;
    % plot(dhhibins(2:end),p5_dlogw_binY2,'r-.');grid on;hold on;
    % plot(dhhibins(2:end),p95_dlogw_binY2,'r-.');grid on;hold on;

    plot([-10,10],[0,0],'k-');hold on;
    plot([p50,p50],[-0.3,0.2],'b-');hold on;
    xlim([0,0.5]);ylim([-0.3,0.1]);
    l = legend('Non-merging firms','Merging firms');
    set(l,'interpreter','latex','location','southwest');

    xlabel('Change in Log Local Labor Market Concentration','interpreter','latex');
    ylabel('Change in log average earnings','interpreter','latex');
    ytickformat('%3.2f');
    xtickformat('%3.1f');
    title('B. Arithmetic average of firm wages','interpreter','latex');
    set(gca,'ticklabelinterpreter','latex','fontsize',12);
    print '-dpng' Figure_Arnold_Figure7
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Table 10 - Again for top ventile
%__________________________________________________________________________
% Taking his heading for Table 10 literally:
logw_ij0            = log(w_ij0);
logw_ij0(n_ij0==0)  = 0;
logwc_ij            = log(wc_ij);
logwc_ij(nc_ij==0)  = 0;

omega_ij0      = n_ij0./sum(n_ij0);
meanlogw_j0    = sum(omega_ij0.*logw_ij0);
n_j0           = sum(n_ij0);
omega_ij1      = nc_ij./sum(nc_ij);
meanlogw_j1    = sum(omega_ij1.*logwc_ij);
n_j1           = sum(nc_ij);
dlogw_j        = meanlogw_j1-meanlogw_j0;       % Change in average log worker earnings
dlogn_j        = log(n_j1)-log(n_j0);           % Log change in market employment
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This matches the dlogHHI from Figure 8A
dlogHHI_cut     = prctile(dlogHHI_j,35); 
fprintf('Mean dlogHHI_j =  %4.3f\n',mean(dlogHHI_j(dlogHHI_j>dlogHHI_cut)));
i_j             = (dlogHHI_j>dlogHHI_cut);
% return
% p50 looks great.
% This should be 95 if doing Arnold (column (3) and column (6)) as I
% haven't figured out how to do weighted regressions in MATLAB :(
% If this is set to p0, then just replicate moment 10
% Up to p75, get correct Above median more negative than Below median
% When go all the way up to p90, this flips though 
% His Figure 7 says "Go up to p95 before you start getting these effects"
% our equivalent Figure 7 says a number much lower, but we also get some
% positive effects when dlogHHI_j is pretty low
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hhi_j_merge    = hhi_j(jstar);
HHI_cut        = prctile(hhi_j_merge(i_j),50); 
% Above median of selected markets *after* applying dlogHHI_cut
%__________________________________________________________________________
% Log C
i_j1           = i_j;
Nmerge         = numel(dlogw_j(i_j1));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j1)';Zerovec],[[ones(size(dlogw_j(i_j1)))';Onevec],[dlogHHI_j(i_j1)';Zerovec]]);
b_elast1_w     = b(2);
b              = regress([dlogn_j(i_j1)';Zerovec],[[ones(size(dlogn_j(i_j1)))';Onevec],[dlogHHI_j(i_j1)';Zerovec]]);
b_elast1_n     = b(2);
%__________________________________________________________________________
% Log C x Above Median C
i_j2           = i_j & (hhi_j_merge>HHI_cut);
Nmerge         = numel(dlogw_j(i_j2));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j2)';Zerovec],[[ones(size(dlogw_j(i_j2)))';Onevec],[dlogHHI_j(i_j2)';Zerovec]]);
b_elast2_w     = b(2); 
b              = regress([dlogn_j(i_j2)';Zerovec],[[ones(size(dlogn_j(i_j2)))';Onevec],[dlogHHI_j(i_j2)';Zerovec]]);
b_elast2_n       = b(2);
%__________________________________________________________________________
% Log C x Below Median C
i_j3           = i_j & (hhi_j_merge<HHI_cut);
Nmerge         = numel(dlogw_j(i_j3));
Ntotal         = floor(Nmerge/factor);
Zerovec        = zeros(Ntotal-Nmerge,1);
Onevec         = ones(Ntotal-Nmerge,1);
b              = regress([dlogw_j(i_j3)';Zerovec],[[ones(size(dlogw_j(i_j3)))';Onevec],[dlogHHI_j(i_j3)';Zerovec]]);
b_elast3_w       = b(2);
b              = regress([dlogn_j(i_j3)';Zerovec],[[ones(size(dlogn_j(i_j3)))';Onevec],[dlogHHI_j(i_j3)';Zerovec]]);
b_elast3_n       = b(2);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('This is where I am playing with the Table 10 results\n');
fprintf('1. Wages\n');
fprintf('BHM: All = %4.3f \t,\t Above median = %4.3f \t,\t Below median = %4.3f\n',[b_elast1_w,b_elast2_w,b_elast3_w]);
fprintf('DA:  All = -0.147 \t,\t Above median = -0.176 \t,\t Below median = 0.058\n');
fprintf('1. Employment\n');
fprintf('BHM: All = %4.3f \t,\t Above median = %4.3f \t,\t Below median = %4.3f\n',[b_elast1_n,b_elast2_n,b_elast3_n]);
fprintf('DA:  All = -0.450 \t,\t Above median = N/A \t,\t Below median = N/A\n');
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');

% lambda calculation
r           = (1/beta) - 1;
nc_j        = sum(nc_ij.^((eta+1)/eta)).^(eta/(eta+1));
Nc          = sum(nc_j.^((theta+1)/theta)).^(theta/(theta+1));
agg_y       = sum(sum(zbar.*z_ij.*(n_ij0.^(alpha))))/(1-Delta); 
agg_yc      = sum(sum(zbar.*z_ij.*(nc_ij.^(alpha))))/(1-Delta);
agg_cc      = agg_yc*(1-(delta/(r+delta))*KS); 
agg_c       = agg_y*(1-(delta/(r+delta))*KS);
lambda      = (agg_cc - ((1/varphibar).^(1/varphi)*(Nc^(1+1/varphi))/(1+1/varphi)) + ((1/varphibar).^(1/varphi)*(N0_index.^(1+1/varphi))/(1+1/varphi)) )/agg_c ;
lambda      = 100*(lambda-1);

% change in concentration
wn_j0       = sum(w_ij0.*n_ij0);
wn_j1       = sum(wc_ij.*nc_ij);
omega_j0    = wn_j0/sum(wn_j0);
omega_j1    = wn_j1/sum(wn_j1);

hhi_0_w     = sum(omega_j0.*hhi_j_merge);
hhi_1_w     = sum(omega_j1.*hhi_j1);
hhi_0_uw    = mean(hhi_j_merge);
hhi_1_uw    = mean(hhi_j1);

dhhi_w      = hhi_1_w   - hhi_0_w;
dhhi_uw     = hhi_1_uw  - hhi_0_uw;
%______________________________________________________

results_random2      = [ factor*mean_dlogn_Y_w,    factor*mean_dlogn_Y_uw,     factor*mean_dlogwn_Y_w ,   ... % Table 3
                        factor*mean_dlogn_Y_w_small , factor*mean_dlogn_Y_w_big  ,          ... % Table 4
                        factor*mean_dlogw_Y    ,                                           ... % Table 5
                        factor*mean_dlogw_HI     , factor*mean_dlogw_MI,  factor*mean_dlogw_LI,                   ... % Table 6 
                        b_dlogw_dlogHHI     , b_dlogn_dlogHHI,                      ... % Table 7
                        b_hhi ,                                                     ... % Table 8
                        mean(dlogHHI_j(dlogHHI_j>dlogHHI_cut)), ...
                        mean_n_pre_Y , median_n_pre_Y , mean_w_pre_Y/1000,               ... % Summary stats
                        lambda, dhhi_w, dhhi_uw, ...
                        0,0,0,0,0,...
                        b_elast1_w,b_elast2_w,b_elast3_w,b_elast1_n, ...
                       ];
                   
%% OPEN TABLE
% Open the file for writing
fid = fopen('../results/tables/table2_arnold_replication.tex', 'w');

% Write the LaTeX table structure
fprintf(fid, '\\begin{table}[t!]\n');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\scalebox{0.85}{\n');
fprintf(fid, '\\begin{tabular}{lcccr}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '\\textbf{Moment} & \\multicolumn{2}{c}{\\textbf{A. Arnold (2020)}} & & \\textbf{B. Model} \\\\ \n');
fprintf(fid, '\\midrule\n');

% Define the correct order of moments matching the original table
printvec = [15, 1, 6, 3, 7, 8, 12];

% Loop through and format table sections according to the correct order
for i = 1:numel(printvec)
    mm = printvec(i);
    
    % Add section headers based on the row number
    if mm == 15
        fprintf(fid, '\\textbf{A. Targeted} & & & & \\\\ \n');
        fprintf(fid, 'Median employment pre-merger & Table 1 & %.0f & & %.0f \\\\ \n', ...
            results_arnold(mm), results_random(mm));
    elseif mm == 1
        fprintf(fid, '\\midrule\n');
        fprintf(fid, '\\textbf{B. Employment and wages} & & & & \\\\ \n');
    elseif mm == 7
        fprintf(fid, '\\midrule\n');
        fprintf(fid, '\\textbf{C. Interaction with concentration} & & & & \\\\ \n');
    end

    % Data rows formatting
    if mm == 12
        % Coefficients, displayed with 3 decimal places
        fprintf(fid, '$\\Delta HHI_j = \\alpha + \\widehat{\\beta}\\Delta{HHI}_j$, $\\widehat{\\beta}$ & Table 8 & %.3f & & %.3f \\\\ \n', ...
            results_arnold(mm), results_random(mm));
    elseif mm == 7 || mm == 8
        % Determine concentration level
        if mm == 7
            concentration = 'High';
        elseif mm == 8
            concentration = 'Medium';
        end
        
        % Values with phantom alignment for small numbers
        fprintf(fid, 'Change in log worker earnings (%s concentration) $(\\times 100)$ & Table 6 & \\phantom{0}%.1f & & \\phantom{0}%.1f \\\\ \n', ...
            concentration, 100 * results_arnold(mm), 100 * results_random(mm));
    elseif mm == 1 || mm == 3 || mm == 6
        % Regular log changes
        if mm == 1
            moment_name = 'Change in log employment $(\\times 100)$';
            table_ref = 3;
        elseif mm == 3
            moment_name = 'Change in log payroll $(\\times 100)$';
            table_ref = 3;
        elseif mm == 6
            moment_name = 'Change in log worker earnings $(\\times 100)$';
            table_ref = 5;
        end
        
        fprintf(fid, '%s & Table %d & %.1f & & %.1f \\\\ \n', ...
            moment_name, table_ref, 100 * results_arnold(mm), 100 * results_random(mm));
    end
end

% Close the table
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}}\n');
fprintf(fid, '\\caption{Mergers and replication of \\citet{Arnold} \\label{tab:arnold}}\n');
fprintf(fid, '\\vspace*{-.3cm}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\\end{table}\n');

% Close the file
fclose(fid);





