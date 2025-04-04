function optim_MCMC_sampling(arg_ind,topX)
% INPUTS:
%    arg_ind:    index of accession to be optimized
%
% OUTPUT:
%    Optimized parameters saved in filen
%    Intermediate solutions saved in saving_path

KE_type='equilibrator';

%% Load all accessions
data=load("data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;

param_name=load_parameter_name();

vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
accession=final_acc(arg_ind);

% Seed random number generator
% rng(0);

%% Load parameters to be optimied

optim_ind=sort(optimized_var_ind(topX));

%% Generation of the structs and options for PESTO
% The structs and the PestoOptions object, which are necessary for the 
% PESTO routines to work are created and set to convenient values

% Prepare parameters structure
display(' Prepare structs and options...')

param_name=[param_name;strcat(param_name(vmaxind),"y23")]; 
parameters.name   = param_name(optim_ind);
parameters.min    = -1*ones(length(parameters.name), 1);
parameters.max    = 1*ones(length(parameters.name), 1);

parameters.number = length(parameters.name);

% PestoOptions
load('data/optionsPesto.mat') % this data is only readable in linux
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

objectiveFunction = @(x)optim_obj_MCMC(x,accession,KE_type,optim_ind);
optionsPesto.obj_type  = 'negative log-posterior'; % if minimization is



%% Parameter Sampling
% Covering all sampling options in one struct
display(' Sampling without prior information...');
optionsPesto.MCMC.nIterations  = 1000;

% PT (with only 1 chain -> AM) specific options:
optionsPesto.MCMC.samplingAlgorithm = 'PT';
optionsPesto.MCMC.PT.nTemps         = 15;
% optionsPesto.MCMC.PT.exponentT      = 6;    
% optionsPesto.MCMC.PT.regFactor      = 1e-8;

% Initialize the chains by choosing a random initial point and a 'large'
% covariance matrix

    
%%

optionsPesto.MCMC.theta0 = zeros(parameters.number,1);
optionsPesto.MCMC.sigma0 = 1000 * eye(parameters.number);

% Intermediate solutions could be stored:
optionsPesto.MCMC.saveEach=1;
optionsPesto.MCMC.saveFileName=strcat('results/parameterization/param_top',char(string(topX)),'/intermediate_',char(accession));

% Run the sampling
tic
parameters = getParameterSamples(parameters, objectiveFunction, optionsPesto);
toc


filen=strcat('results/parameterization/param_top',char(string(topX)),'/MCMCres_',char(accession),'.mat');


save(filen,'parameters')


