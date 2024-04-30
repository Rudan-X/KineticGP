function varargout = optim_obj_MCMC(Variable,acc)

[init_sol,~,~,vmaxind]=optim_initialization_parameters();

nvar=length(init_sol);

parameters21=Variable(1:nvar);
parameters22=Variable(1:nvar);
parameters22(vmaxind)=Variable((nvar+1):end);

%% Transform the variable from logarithmic to original scale
ratio=10.^parameters21;
parameters21=ratio.*init_sol;
ratio=10.^parameters22;
parameters22=ratio.*init_sol;

%%
[~,zACI21,zgs21]=simulate_photosynthesis(parameters21,acc,'ACI');
[~,zAQ22,~]=simulate_photosynthesis(parameters22,acc,'AQ');

varargout{1}=zACI21+zgs21+zAQ22;

varargout{2} = [];
varargout{3} = [];



