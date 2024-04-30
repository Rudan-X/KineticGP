function optim_MCMC_sampling(arg_ind)
TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

% Seed random number generator
rng(0);

addpath('/work/xu2/PESTO-master/','-begin')
% addpath('/home/rudan/Documents/MATLAB/PESTO-master/','-begin')
savepath
% addpath(genpath('/home/rudan/Documents/MATLAB/PESTO-master2/'))
%% Generation of the structs and options for PESTO
% The structs and the PestoOptions object, which are necessary for the 
% PESTO routines to work are created and set to convenient values

% Prepare parameters structure
display(' Prepare structs and options...')

[~,final_acc,param_name,vmaxind]=optim_initialization_parameters();

% global accession and global saving_path will be used by
% PESTO-master/private/performPT.m

accession=final_acc(arg_ind);

global saving_path
saving_path=strcat("/work/xu2/kineticGP/results/optimized_parameters/history_MCMCres_",accession, ".mat");

% saving_path=strcat("../results/optimized_parameters/history_MCMCres_",accession, ".mat");
param_name=[param_name;strcat(param_name(vmaxind),"y22")]; 
parameters.name   = param_name;
parameters.min    = -1*ones(length(parameters.name), 1);
parameters.max    = 1*ones(length(parameters.name), 1);
parameters.number = length(parameters.name);

% PestoOptions
load('../data/optionsPesto.mat')
% optionsPesto           = PestoOptions();
% save('optionsPesto.mat',"optionsPesto")
% optionsPesto.obj_type  = 'negative log-posterior';
optionsPesto.comp_type = 'sequential'; 
optionsPesto.mode      = 'text';
optionsPesto.save      = logical(1);
optionsPesto.tempsave  = logical(1);
PestoOptions.trace     = logical(1);
optionsPesto.objOutNumber=1;

optionsPesto.localOptimizerOptions.Hessian='off';


% objective function

objectiveFunction = @(x)optim_obj_MCMC(x,accession);
optionsPesto.obj_type  = 'negative log-posterior';
display("Chi-square function")


%% Parameter Sampling
% Covering all sampling options in one struct
display(' Sampling without prior information...');
optionsPesto.MCMC.nIterations  = 2000;

% PT (with only 1 chain -> AM) specific options:
optionsPesto.MCMC.samplingAlgorithm = 'PT';
% optionsPesto.MCMC.PT.nTemps         = 6;
% optionsPesto.MCMC.PT.exponentT      = 6;    
% optionsPesto.MCMC.PT.regFactor      = 1e-8;

% Initialize the chains by choosing a random initial point and a 'large'
% covariance matrix


optionsPesto.MCMC.theta0 = zeros(parameters.number,1);
optionsPesto.MCMC.sigma0 = 1000 * eye(parameters.number);

% Run the sampling
tic
parameters = getParameterSamples(parameters, objectiveFunction, optionsPesto);
toc


filen=strcat("../results/optimized_parameters/MCMC_result_",accession,".mat");
save(filen,'parameters')


