function [init_sol,acc_final,param_name,vmaxind]=optim_initialization_parameters()
params=[5.2000    0.0576   0.7200    0.6700    0.0001    2.2820];

[km,vmax,rate_act,permeab]=load_parameter_values_Wang();


init_sol=[params';km';vmax;rate_act;permeab];

data=load('../data/final_accessions.mat');

acc_final=data.final_acc;
acc_final([17,34,69])=[];
param_name=load_parameter_name();

param_name=replace(param_name,'_','-');
param_name(strcmp(param_name,'Vpmax'))='Vpmax';
param_name(strcmp(param_name,'Vmax'))='Vcmax';
param_name(strcmp(param_name,'Vmax-OAA-MC-MCchl[35]'))='Vmax-OAA-MC-MCchl';
param_name(strcmp(param_name,'Vmax-PYR-BC[35]'))='Vmax-PYR-BCchl-BC';
param_name(strcmp(param_name,'Vmax-PYR-MC[35]'))='Vmax-PYR-MC-MCchl';
param_name(strcmp(param_name,'Vmax-PEP-MC[35]'))='Vmax-PEP-MCchl-MC';

vmaxind=find(contains(param_name,'Vmax[1]')):find(contains(param_name,'Vmax-Hep'));
vmaxind=setdiff(vmaxind,find(strcmp(param_name,'ARratio')));
