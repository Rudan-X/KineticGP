function varargout = optim_obj_MCMC(Variable,accession,KE_type,estimated_ind)
% This function calculate the chi-square between simulated and measured
% profiles
% INPUTS:
%    Variable:              Kinetic parameters
%    acc:                   Accession to be optimized
%
%
% OUTPUT:
%    varargout:             varargout{1} is the summation of chi-square of
%                           all fitted curves

init_sol=load_initial_solution();
param_name=load_parameter_name();
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));


% parameters22=Variable(1:nvar,1);
% parameters23=parameters22;
% parameters23(vmaxind)=Variable((nvar+1):end);

parameters0=zeros((length(init_sol)+length(vmaxind)),1);
parameters0(estimated_ind)=Variable;
parameters22=parameters0(1:length(init_sol));
parameters23=parameters22;
parameters23(vmaxind)=parameters0((length(init_sol)+1):end);
%% Transform the variable from logarithmic to original scale

ratio=10.^parameters22;
parameters22=ratio.*init_sol;
ratio=10.^parameters23;
parameters23=ratio.*init_sol;


%% Load experimental data
% Load measured ACI curves

% Measured data are generated using /parameterization/store_measured_curves.m
load("data/processed_data/measured_curves22_23.mat")
data=load("data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;

% Load the profiles for the corresponding accession
[~,argind]=ismember(accession,final_acc);


%%
[simA22,simgs22]=simulate_ACI(parameters22,mean(ACa22.Tair(:,argind)),KE_type);
[simA23,simgs23]=simulate_ACI(parameters23,mean(ACa23.Tair(:,argind)),KE_type);

[simAQ22]=simulate_AQ(parameters22,25,KE_type);
[simAQ23]=simulate_AQ(parameters23,25,KE_type);
% minAstd=min([ACa22.sd(:,argind);ACa23.sd(:,argind)]);
% minGsstd=min([GsCa22.sd(:,argind);GsCa23.sd(:,argind)]);
% minAQstd=min([AQ22.sd(:,argind);AQ23.sd(:,argind)]);

maxAstd=max([ACa22.sd(:,argind);ACa23.sd(:,argind)]);
maxGsstd=max([GsCa22.sd(:,argind);GsCa23.sd(:,argind)]);
maxAQstd=max([AQ22.sd(:,argind);AQ23.sd(:,argind)]);

zACI22=sum((ACa22.meas(:,argind)-simA22).^2./maxAstd^2,'omitnan');
zgs22=sum((GsCa22.meas(:,argind)-simgs22).^2./maxGsstd^2,'omitnan');

zACI23=sum((ACa23.meas(:,argind)-simA23).^2./maxAstd^2,'omitnan');
zgs23=sum((GsCa23.meas(:,argind)-simgs23).^2./maxGsstd^2,'omitnan');

zAQ22=sum((AQ22.meas(:,argind)-simAQ22).^2./maxAstd^2,'omitnan');
zAQ23=sum((AQ23.meas(:,argind)-simAQ23).^2./maxAstd^2,'omitnan');

varargout{1}=zACI22+zgs22+zACI23+zgs23+zAQ22+zAQ23;

varargout{2} = [];
varargout{3} = [];



