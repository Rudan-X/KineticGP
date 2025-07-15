function optim_MCMC_sampling_CI(arg_ind,t)
% INPUTS:
%    arg_ind:    index of accession to be optimized
%
% OUTPUT:
%    Optimized parameters saved in filen
%    Intermediate solutions saved in saving_path

KE_type='equilibrator';
top_x=[10,20,30,40];
topX=top_x(t);
%% Load all accessions
data=load("data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;

param_name=load_parameter_name();

vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
accession=final_acc(arg_ind);

% Seed random number generator
% rng(0);
%% Load parameters to be optimied

optim_ind=optimized_var_ind(topX);



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
optionsPesto.MCMC.nIterations  = 300;

% PT (with only 1 chain -> AM) specific options:
optionsPesto.MCMC.samplingAlgorithm = 'PT';
optionsPesto.MCMC.PT.nTemps         = 15;
% optionsPesto.MCMC.PT.exponentT      = 6;    
% optionsPesto.MCMC.PT.regFactor      = 1e-8;

% Initialize the chains by choosing a random initial point and a 'large'
% covariance matrix


%%
filen=strcat('results/parameterization/param_top',char(string(topX)),'/MCMCres_',char(accession),'.mat');
if exist(filen,'file')==2
    MCMCres=load(filen);
else
    filen=strcat('results/parameterization/param_top',char(string(topX)),'/intermediate_',char(accession),'.mat');
    MCMCres=load(filen);
    MCMCres.parameters.S=MCMCres.res;
end

MCMCres.parameters.S.logPost=abs(MCMCres.parameters.S.logPost);

deg=56-length(optim_ind);
red_chi2=MCMCres.parameters.S.logPost/deg;
valid_indices = find(red_chi2 < 1);

if ~isempty(valid_indices)
    ideal_fval=max(red_chi2(valid_indices));
    index=find(red_chi2==ideal_fval);
    optindex=index(end);
else
    ind=find(MCMCres.parameters.S.logPost==min(MCMCres.parameters.S.logPost));
    optindex=ind(end);
end
pars=MCMCres.parameters.S.par;
optimized_S=pars(:,optindex);
pars=pars(:,~isnan(pars(1,:)));
%%

optionsPesto.MCMC.theta0 = optimized_S;
optionsPesto.MCMC.sigma0 = cov(pars');

% Intermediate solutions could be stored:
optionsPesto.MCMC.saveEach=1;
optionsPesto.MCMC.saveFileName=strcat('results/parameterization/param_top',char(string(topX)),'/intermediateCI_',char(accession));

% Run the sampling
tic
parameters = getParameterSamples(parameters, objectiveFunction, optionsPesto);
toc


filen=strcat('results/parameterization/param_top',char(string(topX)),'/MCMCresCI_',char(accession),'.mat');


save(filen,'parameters')


