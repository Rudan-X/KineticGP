function sensitivity_sampling_eQ(arg_ind,varind)
% INPUTS:
%    arg_ind:    index of accession to be optimized
%
% OUTPUT:
%    Optimized parameters saved in filen
%    Intermediate solutions saved in saving_path

KE_type='equilibrator';
%% Add path of PESTO optimizer and C4 model
% userpath='C:\Users\Rudan\Documents\MATLAB_Drive\';

userpath='/work/xu2/';
addpath(strcat(userpath,'PESTO-master/'),'-begin')
addpath(strcat(userpath,'KineticGP/C4_dynamic_model/'))
addpath(strcat(userpath,'KineticGP/parameterization/'))
addpath(strcat(userpath,'KineticGP/data/'))

%% Load all accessions
data=load("../data/processed_data/final_acc22_23.mat");
[final_acc,~,~,~,~]=load_common_accessions22_23();
param_name=load_parameter_name();

vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));

accession=final_acc(arg_ind);


% Intermediate solutions could be stored by adding the global variable in the code PESTO-master/private/performPT.m
% performPT.m is called by getParameterSamples.m because 'PT' is the chosen MCMC algorithm


% Seed random number generator
% rng(0);


%% Generation of the structs and options for PESTO
% The structs and the PestoOptions object, which are necessary for the 
% PESTO routines to work are created and set to convenient values

% Prepare parameters structure
display(' Prepare structs and options...')

param_name=[param_name;strcat(param_name(vmaxind),"y23")]; 
parameters.name   = param_name(varind);
parameters.min    = -1*ones(length(parameters.name), 1);
parameters.max    = 1*ones(length(parameters.name), 1);

parameters.number = length(parameters.name);

% PestoOptions
load('../data/optionsPesto.mat') % this data is only readable in linux
% optionsPesto           = PestoOptions();
% save('optionsPesto.mat',"../data/optionsPesto")

optionsPesto.comp_type = 'sequential'; 
optionsPesto.mode      = 'text';
optionsPesto.save      = logical(1);
optionsPesto.tempsave  = logical(1);
PestoOptions.trace     = logical(1);
optionsPesto.objOutNumber=1;

optionsPesto.localOptimizerOptions.Hessian='off';


% objective function

objectiveFunction = @(x)sensitivity_obj(x,accession,KE_type,varind);
optionsPesto.obj_type  = 'negative log-posterior'; % if minimization is



%% Parameter Sampling
% Covering all sampling options in one struct
display(' Sampling without prior information...');
optionsPesto.MCMC.nIterations  = 500;

% PT (with only 1 chain -> AM) specific options:
optionsPesto.MCMC.samplingAlgorithm = 'PT';
optionsPesto.MCMC.PT.nTemps         = 4;
% optionsPesto.MCMC.PT.exponentT      = 6;    
% optionsPesto.MCMC.PT.regFactor      = 1e-8;

% Initialize the chains by choosing a random initial point and a 'large'
% covariance matrix

    
%%

optionsPesto.MCMC.theta0 = zeros(parameters.number,1);
optionsPesto.MCMC.sigma0 = 1000 * eye(parameters.number);

optionsPesto.MCMC.saveEach=1;
optionsPesto.MCMC.saveFileName=strcat('../SensitivityAnalysis/sensitivity_results/intermediate_',char(accession),char(parameters.name));

% Run the sampling
tic
parameters = getParameterSamples(parameters, objectiveFunction, optionsPesto);
toc


filen=strcat('../SensitivityAnalysis/sensitivity_results/MCMCres_',accession,'_',char(parameters.name),'.mat');
save(filen,'parameters')


