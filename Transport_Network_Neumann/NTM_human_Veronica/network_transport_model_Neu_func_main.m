function [N, M, F_in, F_out, M_tot,Mass_error,Flux_tot,i_zero] = network_transport_model_Neu_func_main(matdir,varargin) %(matdir,varargin)

if nargin < 4
    matdir = [cd filesep 'MatFiles'];
 %matdir = [cd filesep 'MatFiles'];
end
alpha_ = 0;
F_edge_0_ = 0;
beta_ = 1e-06;
gamma1_ =0.001;
gamma2_ = 0;
delta_ = 50;
epsilon_ = 50; 
lambda1_ = 1e-02;
%lambda2_ = 1e-02;
gamma1_dt_=0;
lambda1_dt_=0;
lambda_time_dip_roi_=[2 6 7 14 2+15 6+15 7+15 14+15];
gamma_time_dip_roi_=[2 6 7 14 2+15 6+15 7+15 14+15];
%study_ = 'Hurtado';
init_path_ = [];
init_rescale_ = 2e-2;
dt_ =0.01;
T_ = 0.1;

trange_ = [0:0.0005:0.002, 0.0025:0.0025:0.1, 0.11:0.01:6];

length(trange_)

frac_ = 0.7;
%L_int_ = 1000; % in micrometers - SET to Matrix from Human Dist. Template!!!
L1_ = 200;
%L2_ = 200; 
L_ais_ = 40;
L_syn_ = 40;
resmesh_ = 'fine';
plotting_ = 1;
reltol_ = 1e-6;
abstol_ = 1e-6;
fsolvetol_ = 1e-20;
connectome_subset_ = 'Human_All';
time_scale_ = 6*30*24*(60)^2;
len_scale_ = 1e-3;

mu_r_0_ = 1.5;
mu_r_L_ = mu_r_0_;
mu_u_0_ = 1.5;
mu_u_L_ = mu_u_0_;

ip = inputParser;
% validChar = @(x) ischar(x);
validScalar = @(x) isnumeric(x) && isscalar(x)&& (x>=0);
validLogical = @(x) validScalar(x) && (x == 0 || x == 1);
%validAxonDiv=@(x) strcmp(x,'r1') | strcmp(x,'half') | strcmp(x,'r2');
addParameter(ip, 'alpha', alpha_, validScalar);
addParameter(ip, 'F_edge_0', F_edge_0_, validScalar);
addParameter(ip, 'beta', beta_, validScalar);
addParameter(ip, 'gamma1', gamma1_, validScalar);
addParameter(ip, 'gamma2', gamma2_, validScalar);
addParameter(ip, 'delta', delta_, validScalar);
addParameter(ip, 'epsilon', epsilon_, validScalar);
addParameter(ip, 'frac', frac_, validScalar);
addParameter(ip, 'lambda1', lambda1_, validScalar);
%addParameter(ip, 'lambda2', lambda2_, validScalar);
addParameter(ip, 'gamma1_dt', gamma1_dt_, validScalar);
addParameter(ip, 'lambda1_dt', lambda1_dt_, validScalar);
addParameter(ip, 'lambda_time_dip_roi', lambda_time_dip_roi_);
addParameter(ip, 'gamma_time_dip_roi', gamma_time_dip_roi_);
addParameter(ip, 'L1', L1_, validScalar);
%addParameter(ip, 'L2', L2_, validScalar);
addParameter(ip, 'resmesh', resmesh_);
addParameter(ip, 'L_ais', L_ais_);
addParameter(ip, 'L_syn', L_syn_);
addParameter(ip, 'reltol', reltol_, validScalar);
addParameter(ip, 'abstol', abstol_, validScalar);
addParameter(ip, 'fsolvetol', fsolvetol_, validScalar);
addParameter(ip, 'connectome_subset', connectome_subset_);
addParameter(ip, 'len_scale', len_scale_, validScalar);
addParameter(ip, 'time_scale', time_scale_, validScalar);
addParameter(ip, 'mu_r_0',mu_r_0_);
addParameter(ip, 'mu_r_L',mu_r_L_);
addParameter(ip, 'mu_u_0',mu_u_0_);
addParameter(ip, 'mu_u_L',mu_u_L_);
%addParameter(ip, 'study', study_);
addParameter(ip, 'init_rescale', init_rescale_, validScalar);
addParameter(ip, 'dt', dt_);
addParameter(ip, 'T', T_);
addParameter(ip, 'trange', trange_);
addParameter(ip, 'init_path', init_path_);
addParameter(ip, 'plotting', plotting_, validLogical);
parse(ip, varargin{:});

%load([matdir filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct'); 
%load([matdir filesep 'DefaultAtlas.mat'],'DefaultAtlas'); % NEED TO ADD HUMAN VOLUMES
load([matdir filesep 'Human-DK-Vols_mm3.mat'], 'volumes');
load([matdir filesep 'Human-DK-dist.mat'], 'distances');

%distances = distances * 10^3;

%load([matdir filesep 'CCF_labels.mat'],'CCF_labels');
load([matdir filesep 'Connectome-Human-Directed-DK.mat'],'Consensus');
Conn = Consensus.SC;
Conn = Conn * 5.3522e+03;

Adj = logical(Conn);

% ADJUST TO REAL VOLUMES

for i=1:length(Conn)
    for j=1:length(Conn)
       if i==j
           Conn(i,j)=0;
       end
    end
end

thresh_C = 0.8 * mean(nonzeros(Conn(:)));
Conn(Conn < thresh_C) = 0;
Adj = logical(Conn);

% Conn = readmatrix([matdir filesep 'mouse_connectome_19_01.csv']);
% Adj = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);

% if ~isempty(ip.Results.init_path)
%     init_path = zeros(size(Conn,1),1);
%     for i = 1:length(ip.Results.init_path)
%         reghemstr = ip.Results.init_path{i};
%         reghemcell = split(reghemstr,'_');
%         reglog = ismember(CCF_labels(:,1),reghemcell{1});
%         if strcmp(reghemcell{2},'L')
%             hemlog = ismember(CCF_labels(:,4),'Left Hemisphere'); 
%         elseif strcmp(reghemcell{2},'R')
%             hemlog = ismember(CCF_labels(:,4),'Right Hemisphere'); 
%         else 
%             hemlog = ones(size(Conn,1),1);
%         end
%         init_path((reglog + hemlog) == 2) = 1;
%     end
% elseif isnan(mousedata_struct.(ip.Results.study).seed)
%     %init_path = logical(mousedata_struct.(ip.Results.study).data(:,1));
%     init_path = mousedata_struct.(ip.Results.study).data(:,1);
%     mean_nz = mean(nonzeros(init_path));
%     init_path = init_path ./ mean_nz; % For IbaP301S study where seed is not specified and data is nonbinary at t=1
%     init_path = DataToCCF(init_path,ip.Results.study,matdir);
% else
%     init_path = logical(mousedata_struct.(ip.Results.study).seed);
%     init_path = DataToCCF(init_path,ip.Results.study,matdir);
% end

% Initial seed pathology

init_path = zeros(size(Conn,1),1);
init_path(10,1) = 1;

beta_new = ip.Results.beta * ip.Results.time_scale;
gamma1_new = ip.Results.gamma1 * ip.Results.time_scale;
gamma2_new = ip.Results.gamma2 * ip.Results.time_scale;
%taufun = @(x) ip.Results.init_rescale - (x +...
    %(gamma1_new * x.^2)./(beta_new - gamma2_new * x));
mu_frac=ip.Results.mu_r_0/ip.Results.mu_u_0;
taufun = @(x) ip.Results.init_rescale - (x*(1+mu_frac) +...
    (gamma1_new *beta_new*(1+mu_frac^2)* x.^2-gamma1_new*gamma2_new*mu_frac*x.^3*(1+mu_frac))./((beta_new - gamma2_new * x)*(beta_new - gamma2_new *mu_frac* x)));
options_taufun = optimset('TolFun',ip.Results.fsolvetol,'Display','off');
init_rescale_n =fsolve(taufun,0,options_taufun); %2e-2;


id_nan=isnan(init_path);
init_path(id_nan)=0;

init_tau = init_rescale_n * init_path;

% connectome subset needs to be adjusted for human

% switch ip.Results.connectome_subset
%     case 'Hippocampus'
%         inds = ismember(CCF_labels(:,3),'Hippocampus');
%     case 'Hippocampus+PC+RSP'
%         inds_hipp = ismember(CCF_labels(:,3),'Hippocampus');
%         inds_pc = ismember(CCF_labels(:,1),'Piriform area');
%         inds_rsp = ismember(CCF_labels(:,3),'Retrosplenial Area');
%         inds = logical(inds_hipp + inds_pc + inds_rsp);
%     case 'RH'
%         inds = ismember(CCF_labels(:,4),'Right Hemisphere');
%     case 'LH'
%         inds = ismember(CCF_labels(:,4),'Left Hemisphere');
%     otherwise
%         inds = logical(ones(size(Conn,1),1)); %#ok<LOGL> 
% end

inds = logical(ones(size(Conn,1),1));

Adj = Adj(inds,inds);
Conn = Conn(inds,inds);

%Vol = DefaultAtlas.volumes(inds);

Vol = volumes(inds)*8; % Volumes in mm^3 - MAKE SURE UNITS ARE CORRECT

init_tau = init_tau(inds);
nroi = size(Adj,1);

% dont know what these were used for

% regnamecell = CCF_labels(inds,:);
% regnames = cell(size(regnamecell,1),1);
% for i = 1:length(regnames)
%     regname = regnamecell{i,1};
%     reghem = regnamecell{i,4};
%     if strcmp(reghem,'Right Hemisphere')
%         regnames{i} = [regname ' RH'];
%     else
%         regnames{i} = [regname ' LH'];
%     end
% end

if isempty(ip.Results.trange)
    t = 0:ip.Results.dt:ip.Results.T;
else
    t = ip.Results.trange;
end
   nt=length(t);
  N=zeros(nroi,nt);
 
  N(:,1)=init_tau(:);

  % for i=1:nroi
  %     if N(i,1)==0
  %         N(i,1)=1e-4;
  %     end
  % end
    seedregions_ = (N(:,1) > 0);  %    (N(:,1) > 1e-4); % 
 
 i_zero=N(:,1)==0;   %init_tau==0;

%Regional_edge_mass = zeros([nroi, nroi, nt, 4]);
 Fun_N=zeros(nroi,nt);

 netw_flux_0= zeros([nroi,size(N)]);
  netw_flux_L= zeros([nroi,size(N)]);
  Mass_edge=zeros([nroi,nroi,nt]);
  n_ss0=zeros([nroi,nroi,nt]);
 
F_source_edge=zeros([nroi,nroi,nt]);
 %gamma_time_dip_roi=[2 6 7 14 2+15 6+15 7+15 14+15];
 roi_id=zeros(nroi,1);
 roi_id( ip.Results.gamma_time_dip_roi,1)=1;
%lambda_time_dip_roi=[2 6 7 14 2+15 6+15 7+15 14+15];
roi_id_lambda=zeros(nroi,1);
 roi_id_lambda( ip.Results.lambda_time_dip_roi,1)=1;
 Gamma1_der=@(t_)ip.Results.gamma1_dt*ip.Results.time_scale*roi_id;
 gamma1_fun = @(t_) gamma1_new+ t_.*Gamma1_der(t_);
 Lambda1_der=@(t_)ip.Results.lambda1_dt*roi_id_lambda;
 lambda1_fun=@(t_) ip.Results.lambda1 + t_.*Lambda1_der(t_);


Idx_netw_x0=logical(init_tau).*Adj;
Idx_netw_xL= logical(init_tau);  



 
    function [Fun_N_app,netw_flux_0_app,netw_flux_L_app,F_source_edge_app,n_ss0_app,Mass_edge_app]=NTM_nodes(N_adj_app,N_app,n0_guess,Gamma1_x0_app,Lambda1_x0_app,t_step)
       
      [netw_flux_0_app,netw_flux_L_app,F_source_edge_app,n_ss0_app,Mass_edge_app]=NetworkFluxCalculator_Neu(N_adj_app,N_app,n0_guess,Adj,'beta',ip.Results.beta,...
                                    'delta',ip.Results.delta,'F_edge_0',ip.Results.F_edge_0,...
                                    'epsilon',ip.Results.epsilon,...
                                    'frac',ip.Results.frac,...
                                    'lambda1_x0',Lambda1_x0_app,...
                                    'lambda1_xL',(lambda1_fun(t_step)),...
                                    'gamma1_x0',Gamma1_x0_app,...
                                    'gamma1_xL',gamma1_fun(t_step),...
                                    'idx_netw_x0',Idx_netw_x0,'idx_netw_xL',Idx_netw_xL,'frac',ip.Results.frac,...
                                    'L1',ip.Results.L1,...
                                    'L_ais',ip.Results.L_ais,...
                                    'L_syn',ip.Results.L_syn,...
                                    'resmesh',ip.Results.resmesh,...
                                    'reltol',ip.Results.reltol,...
                                    'abstol',ip.Results.abstol,...
                                    'fsolvetol',ip.Results.fsolvetol, ...
                                    'time_scale', ip.Results.time_scale, ...
                                    'mu_r_0',ip.Results.mu_r_0,'mu_r_L',ip.Results.mu_r_L, 'mu_u_0',ip.Results.mu_u_0, 'mu_u_L',ip.Results.mu_u_L ); %compute the steady state  network flux at time t0


         m_t_app= (gamma1_fun(t_step)).*N_app.*((2*beta_new-ip.Results.gamma2.*N_app)./(beta_new-ip.Results.gamma2.*N_app).^2);
         F_in=netw_flux_L_app;
     F_out=netw_flux_0_app;
     F_out=F_out.';
     Fun_N_app=(1./(Vol.*(1+m_t_app))).*( ( diag((Conn.'*F_in)) - diag((Conn*F_out))) );
    end
 
%der=zeros(nroi,nt);
%fprintf('Flux calculation at initial time \n') 
% tic
n0_guess=0.0001*ones(nroi);
for h=1:  (nt-1)
  
  fprintf('Time step %d/%d\n',h,nt-1) 
   N_adj_in=N(:,h).*Adj;

Gamma1_x0_h=gamma1_fun(t(h)).*Adj;
Lambda1_x0_h=lambda1_fun(t(h)).*Adj;

  [Fun_N(:,h),netw_flux_0(:,:,h),netw_flux_L(:,:,h),F_source_edge(:,:,h),n_ss0(:,:,h),Mass_edge(:,:,h)]=NTM_nodes(N_adj_in,N(:,h),n0_guess,Gamma1_x0_h,Lambda1_x0_h,t(h));

n0_guess=n_ss0(:,:,h);

 
 

 if h<=5000 % previously 5
     N_k1=N(:,h)+Fun_N(:,h).*((t(h+1)-t(h))/2);
    N_adj_k1=N_k1.*Adj;

  Gamma1_x0_h_k1=(gamma1_fun((t(h+1)+t(h))/2)).*Adj;
 
  Lambda1_x0_h_k1=(lambda1_fun((t(h+1)+t(h))/2)).*Adj;


[Fun_N_k1,~,~]=NTM_nodes(N_adj_k1,N_k1,n0_guess,Gamma1_x0_h_k1, Lambda1_x0_h_k1,(t(h+1)+t(h))/2);

 N(:,h+1)=N(:,h)+(Fun_N_k1)*(t(h+1)-t(h));


  else
  N(:,h+1)=N(:,h)+Fun_N(:,h).*(t(h+1)-t(h));
  end
   


 
 if h==nt-1

      N_adj_h1=N(:,h+1).*Adj;
  
  Gamma1_x0_h1=gamma1_fun(t(h+1)).*Adj;
 
  Lambda1_x0_h1=lambda1_fun(t(h+1)).*Adj;


[netw_flux_0(:,:,h+1),netw_flux_L(:,:,h+1),F_source_edge(:,:,h+1),n_ss0(:,:,h+1),Mass_edge(:,:,h+1)]=NetworkFluxCalculator_Neu(N_adj_h1,N(:,h+1),n_ss0(:,:,h),Adj,'beta',ip.Results.beta,...
                                    'delta',ip.Results.delta,'F_edge_0',ip.Results.F_edge_0,...
                                    'epsilon',ip.Results.epsilon,...
                                    'frac',ip.Results.frac,...
                                    'lambda1_x0',Lambda1_x0_h1,...
                                    'lambda1_xL',lambda1_fun(t(h+1)),...
                                    'gamma1_x0',Gamma1_x0_h1,...
                                    'gamma1_xL',gamma1_fun(t(h+1)),...
                                    'idx_netw_x0',Idx_netw_x0,'idx_netw_xL',Idx_netw_xL,'frac',ip.Results.frac,...
                                    'L1',ip.Results.L1,...
                                    'L_ais',ip.Results.L_ais,...
                                    'L_syn',ip.Results.L_syn,...
                                    'resmesh',ip.Results.resmesh,...
                                    'reltol',ip.Results.reltol,...
                                    'abstol',ip.Results.abstol,...
                                    'fsolvetol',ip.Results.fsolvetol, ...
                                    'time_scale', ip.Results.time_scale, 'connectome_subset',ip.Results.connectome_subset, ...
                                    'mu_r_0',ip.Results.mu_r_0,'mu_r_L',ip.Results.mu_r_L, 'mu_u_0',ip.Results.mu_u_0, 'mu_u_L',ip.Results.mu_u_L ); %compute the steady state  network flux at time t0

 end

N(:,h+1) = N(:,h+1) + ip.Results.alpha*(t(h+1)-t(h))*N(:,h);
 toc

end

 M=(gamma1_fun(t).* N.^2)./(beta_new-ip.Results.gamma2 * N);
  N_int=mu_frac*N;
  M_int=(gamma1_fun(t).* N_int.^2)./(beta_new-ip.Results.gamma2 * N_int);

 

 Flux_tot=zeros(1,nt);
 for h=1:nt
      F_in=netw_flux_L(:,:,h);
    F_out=netw_flux_0(:,:,h);
    F_out=F_out.';
  Flux_tot(1,h)=sum(( diag((Conn.'*F_in)) - diag((Conn*F_out))));
 end
  fprintf('Total Flux = %d\n',Flux_tot)
      figure(1)
    plot(t,Flux_tot,'r')
    
    ylabel('t')
    xlabel('Flux_tot(t)')
    title('Total flux')
      %subtitle(txt,'Interpreter','latex');

  Mass_node = sum(Vol.*N,1)+sum(Vol.*M,1)+sum(Vol.*N_int,1)+sum(Vol.*M_int,1);
    % size(Mass_node)
     fprintf('Total Node Mass = %d\n',Mass_node)
    Mass_tot_edge=zeros(1,nt);
    for h=1:(nt)
        Mass_tot_edge(1,h) = sum(Conn.*Mass_edge(:,:,h),'all');
    end
    fprintf('Total Edge Mass = %d\n',Mass_tot_edge)
    for h=1:nt
        F=sum(F_source_edge(:,:,h),'all');
    end
     F_source_ed=zeros(1,nt);
    for h=1:nt
      F_source_ed(1,h)=sum(Conn.*F_source_edge(:,:,h),'all').*t(h);
    end
    F_source_nodes=zeros(1,nt);
    source_reg=ones(nroi,1);
    for h=1:nt
       F_source_nodes(1,h)=sum(Vol(seedregions_).*(ip.Results.alpha*source_reg(seedregions_,1)) ,'all').*t(h);
    end
     M_tot=Mass_node+Mass_tot_edge ;
     F_source=F_source_nodes+F_source_ed;
     M_tot_tru=M_tot(1,1)+F_source;
     Mass_error=(M_tot-(M_tot(1,1)+F_source))./(M_tot(1,1)+F_source);

    txt = ['$\mathbf{\gamma}_{1,t}$' num2str(ip.Results.gamma1_dt),',', '$\mathbf{\lambda}_{t}$',',' num2str(ip.Results.lambda1_dt),',' '$\mathbf{\lambda_1} = $' num2str(lambda1_)  ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
    figure(2)
    plot(t,M_tot,'r')
    hold on
    plot(t,M_tot_tru,'b')
    ylabel('t')
    xlabel('Mass(t)')
    title('Total mass vs total mass true')
      subtitle(txt,'Interpreter','latex');
      %saveas(figure(1),[ cd '/plot_tau/' 'Tot_mass' '_' 'wrong' '.png' ])
%     M_tot(1,1)
%   reldiffs = M_tot - M_tot(1,1);
% reldiffs = reldiffs ./ M_tot(1,1);
% 
 figure (3); 
% %hold on;
% %for i = 1:size(masstots,1)
     plot(t,Mass_error); 
     ylabel('t')
     xlabel('E(t)')
     title('Relative error Total mass')
  subtitle(txt,'Interpreter','latex');
%     % saveas(figure(2),[ cd '/plot_tau/' 'Rel_err' '_' 'wrong' '.png' ])
% %end
% 
% 
% 
% 
figure(4)
%figure
subplot(2,1,1)
plot(t,N);


xlabel('t');
ylabel('N(t)');
 title("N,M distributions on the network",'Fontsize',12);
  %txt = ['$\mathbf{\gamma}_{1,t} = $' num2str(gamma1_dt), ',','$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
 %subtitle(txt,'Interpreter','latex');
subplot(2,1,2)
plot(t,M);
xlabel('t');
ylabel('M(t)');
  %saveas(figure(3),[ cd '/plot_tau/' 'N_M' '_' 'wrong' '.png' ])
% 

figure(5)
%figure
subplot(2,1,1)
plot(t,N(i_zero,:));
 title("N,M distributions on the network",'Fontsize',12);
 % txt = ['$\mathbf{\gamma}_{1,t}$' num2str(gamma1_dt),',', '$\mathbf{\lambda}_{t}$',',' num2str(lambda_dt),',' '$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
 %subtitle(txt,'Interpreter','latex');
subplot(2,1,2)
plot(t,M(i_zero,:));
  %saveas(figure(4),[ cd '/plot_tau/' 'N_M_no_seed' '_' 'wrong' '.png' ])

% figure(5)
% %figure
% 
% plot(t,res_max);
% 
% 
% xlabel('t');
% ylabel('res_max');
%  title("max residual shooting parameter calculation",'Fontsize',12);
%    saveas(figure(4),[ cd '/plot_tau/' 'max_shooting' 'wrong' '.png' ])




end