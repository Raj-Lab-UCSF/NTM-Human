%% 2. Parameter definitions
simstr = 'sim_tau_human_1_';

%inputparams = cell(2,8);
inputparams = cell(1,10);

paramnames = {'gamma1', 'delta','epsilon','lambda1',...
    'gamma1_dt','lambda1_dt','mu_r_0','mu_u_0','F_edge_0','alpha'};
inputparams(1,:) = paramnames;
% inputparams{2,1} = 1e-6; % beta
% inputparams{2,2} = [0.002]; % gamma1
% inputparams{2,3} = [ 50]; % delta
% inputparams{2,4} = [25 ]; % epsilon
% inputparams{2,5} = [0.005 0.1]; % lambda1
% inputparams{2,6} = [0.005 0.1]; % lambda2
% inputparams{2,7} = [0]; % gamma1_dt
% inputparams{2,8} = [0]; % lambda_dt

% 2b. Create parameter array to grid search using allcomb()
% paramgrid = combvec(inputparams{2,1},...
%                     inputparams{2,2},...
%                     inputparams{2,3},...
%                     inputparams{2,4},...
%                     inputparams{2,5},...
%                     inputparams{2,6},...
%                     inputparams{2,7},...
%                     inputparams{2,8});


 %[1e-6, 0.001,1,0.1, 0.025,0.025,0,0,;
 paramgrid=   [ 0.001,50,50, 0.001, 0,0,1.5,1.5,0,0;
     0.002,50,50, 0.001,0,0,1.5,1.5,0,0;] 
     %0.001,1,0.1, 0.005,0.005,0,0,1.5,1.5,1.5,1.5,0,0;
     %0.002,0.1,1, 0.005,0.005,0,0,1.5,1.5,1.5,1.5,0,0];

%paramgrid=paramgrid.';

paramnamescell = repmat(paramnames,size(paramgrid,1),1);



 beta = 1e-06;
% gamma1 = 1e-05;%0.001; %2e-03;
gamma2 = 0;
 gamma_time_dip_roi=[2 6 7 14 2+15 6+15 7+15 14+15];
lambda_time_dip_roi=[2 6 7 14 2+15 6+15 7+15 14+15];
%study = 'Hurtado';
init_path = [];
init_rescale = 2e-2;
dt =0.01;
T = 0.1;
trange =[0:0.0005: 0.002, 0.0025:0.0025:0.1, 0.11:0.01:6];

frac = 0.7; % Average fraction of n diffusing 
% in micrometers
L1 = 200;

L_ais = 40;
L_syn= 40;
resmesh = 'fine';
plotting= 1;
reltol= 1e-8;
abstol= 1e-8;
fsolvetol= 1e-12;
connectome_subset= 'Human_All';
time_scale= 6*30*24*(60)^2;
len_scale= 1e-3;
matdir = [cd filesep 'MatFiles'];
ncores=2;
parpool(ncores)
tic
parfor i = 1:size(paramgrid,1)
    fprintf('Simulation %d/%d \n',i,size(paramgrid,1))
    paramlist = paramgrid(i,:);
    paramnames_i = paramnamescell(i,:); 
 [N, M, F_in, F_out,M_tot,Mass_error,Flux_tot,i_zero ]= network_transport_model_Neu_func_main(matdir, paramnames_i{1},paramlist(1),...
                                paramnames_i{2},paramlist(2),...
                                paramnames_i{3},paramlist(3),...
                                paramnames_i{4},paramlist(4),...
                                paramnames_i{5},paramlist(5),...
                                paramnames_i{6},paramlist(6),...
                                paramnames_i{7},paramlist(7),...
                                paramnames_i{8},paramlist(8),...    %'beta',beta,'gamma1',gamma1,'gamma2',gamma2,'delta',delta,'epsilon',epsilon,'lambda1',lambda1,'lambda2',lambda2,'gamma1_dt',gamma1_dt,'lambda_dt',lambda_dt, ...
                                 paramnames_i{9},paramlist(9),paramnames_i{10},paramlist(10), ...
                                 'L1',L1,...
                              'L_ais',L_ais,...
                                'L_syn',L_syn,...
                                'T',T,...
                                'dt',dt,...
                                'trange',trange,...
                                'resmesh', resmesh,...
                                'plotting',plotting,...
                                'reltol',reltol,...
                                'abstol',abstol,...
                                'fsolvetol',fsolvetol,...
                                'init_rescale',init_rescale,...
                                'init_path',init_path,...
                                'connectome_subset',connectome_subset, ...
                                'lambda_time_dip_roi',lambda_time_dip_roi,'gamma_time_dip_roi',gamma_time_dip_roi);


%M_tot=Mass_node+Mass_tot_edge;
t=trange;
   sim_no=num2str(i);
    %txt = ['$\mathbf{\gamma}_{1,t}$' num2str(gamma1_dt),',','$\mathbf{\gamma}_{1}$' num2str(gamma1),',', '$\mathbf{\lambda}_{t}$',',' num2str(lambda_dt),',' '$\mathbf{\lambda_1} = $' num2str(lambda1) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2) ',' , '$\mathbf{\beta}=$' num2str(beta) ','  '$\mathbf{\epsilon}=$' num2str(epsilon) ',', '$\mathbf{\delta}=$' num2str(delta) ];
    txt = ['$\mathbf{\gamma}_{1,t}$' num2str(paramlist(5)),',','$\mathbf{\gamma}_{1}$' num2str(paramlist(1)),',', '$\mathbf{\lambda}_{t}$' num2str(paramlist(6)),',' ,'$\mathbf{\lambda_1} = $' num2str(paramlist(4)), ',' ,'$\mathbf{\epsilon}=$' num2str(paramlist(3)) ,',', '$\mathbf{\delta}=$' num2str(paramlist(2)),',', '$\mathbf{\mu}_{r}$' num2str(paramlist(7)),',' ,'$\mathbf{\mu}_{r}$' num2str(paramlist(8)),',',  '$\mathbf{F_edge}_{0}$' num2str(paramlist(9))];
    
    figure(1)
    plot(t,M_tot)
    ylabel('t')
    xlabel('Mass(t)')
    title('Total mass')
     subtitle(txt,'Interpreter','latex');
   
     saveas(figure(1),[ cd '/plot/' simstr   sim_no  'Tot_mass' '_' 'Neu'  '.png' ])
   % M_tot(1,1)
  %reldiffs = M_tot - M_tot(1,1);
%reldiffs = reldiffs ./ M_tot(1,1);

figure(2); 
%hold on;

    plot(t,Mass_error); 
    ylabel('t')
    xlabel('E(t)')
    title('Relative error Total mass')
 subtitle(txt,'Interpreter','latex');
    saveas(figure(2),[ cd '/plot/' simstr sim_no 'Rel_err' '_' 'Neu'  '.png' ])

figure(3); 
%hold on;

    plot(t,Flux_tot); 
    ylabel('t')
    xlabel('Sum J')
    title('Sum Fluxes')
 subtitle(txt,'Interpreter','latex');
    saveas(figure(3),[ cd '/plot/' simstr sim_no 'Sum_fluxes' '_' 'Neu'  '.png' ])




figure(4)
%figure
subplot(2,1,1)
plot(t,N);


xlabel('t');
ylabel('N(t)');
 title("N,M distributions on the network",'Fontsize',12);
  %txt = ['$\mathbf{\gamma}_{1,t} = $' num2str(gamma1_dt), ',','$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
 subtitle(txt,'Interpreter','latex');
subplot(2,1,2)
plot(t,M);
xlabel('t');
ylabel('M(t)');
 saveas(figure(4),[ cd '/plot/' simstr  sim_no 'N_M' '_' 'Neu'  '.png' ])

figure(5)
%figure
subplot(2,1,1)
plot(t,N(i_zero,:));
 title("N,M distributions on the network",'Fontsize',12);
 % txt = ['$\mathbf{\gamma}_{1,t}$' num2str(gamma1_dt),',', '$\mathbf{\lambda}_{t}$',',' num2str(lambda_dt),',' '$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
 subtitle(txt,'Interpreter','latex');
subplot(2,1,2)
plot(t,M(i_zero,:));
 saveas(figure(5),[ cd '/plot/' simstr sim_no 'N_M_no_seed' '_' 'Neu'  '.png' ])

figure(6)
%figure
plot(t,N+M);


xlabel('t');
ylabel('N(t)+M(t)');
 title("N,M distributions on the network",'Fontsize',12);
  %txt = ['$\mathbf{\gamma}_{1,t}$' num2str(gamma1_dt),',','$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
 subtitle(txt,'Interpreter','latex');
 saveas(figure(6),[ cd '/plot/'  simstr sim_no 'Tot_tau' '_' 'Neu'  '.png' ])

 
figure(7)
plot(t,N(i_zero,:)+M(i_zero,:));


xlabel('t');
ylabel('N(t)+M(t)');
 title("Total Tau on the network",'Fontsize',12);
  %txt = ['$\mathbf{\gamma}_{1,t}$' num2str(gamma1_dt),',','$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
 subtitle(txt,'Interpreter','latex');
 saveas(figure(7),[ cd '/plot/' simstr   sim_no 'Tot_tau_no_seed' '_' 'Neu'  '.png' ])
%nroi_idx=1:nroi;
  %figure(6)
figure(8)
%heatmap(N(i_zero,1:end)+M(i_zero,1:end));
heatmap(N(i_zero,1:end));

xlabel('t');
ylabel('N(t)');
 title("Soluble Tau on the network");
% subtitle(txt,'Interpreter','latex');
 saveas(figure(8),[ cd '/plot/'  simstr sim_no 'Sol_tau_heatmap' '_' 'Neu'  '.png' ])
 %  txt = ['$\mathbf{\gamma}_{1,t}$' num2str(gamma1_dt),',','$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
 %subtitle(txt,'Interpreter','latex');


end

delete(gcp('nocreate'));




% figure(8)
% for h=1:length(t)
% 
%     subplot(2,1,1)
%     %heatmap(1:nroi,1:nroi, B(:,:,h))
%     plot(t(h),max(abs(B(:,:,h)),[],'all'),'o',LineWidth=2)
%     hold on
%     title('B network')
%     %subtitle(txt,'Interpreter','latex')
%     subplot(2,1,2)
%     %heatmap(1:nroi,1:nroi,6*30*24*(60)^2*netw_flux(:,:,h))
% 
%     plot(t(h),max(abs(6*30*24*(60)^2*netw_flux(:,:,h)),[],'all'),'*',LineWidth=2)
%    hold on
%  % xlabel('nroi')
%  % ylabel('J,B (o-* line)' )
%  title('J network')
% 
% end
% % saveas(figure(8),[ cd '/plot/' 'B_J' '_' 'lambda_var' '.png' ])
% saveas(figure(8),[ cd '/plot/' 'BJ' '_' 'gamma_lambda_var_new' '.png' ])
