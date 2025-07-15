function SFsec41_calculate_CIs()

userpath='/work/xu2/';
addpath(genpath(strcat(userpath,'PESTO-master/')))

addpath(genpath(strcat(userpath,'KineticGP/code')))
addpath(genpath(strcat(userpath,'KineticGP/data')))

%%
top_x=[10,20,30,40];
for t=1:4
    topX=top_x(t);
    optim_ind=optimized_var_ind(topX);

    data=load("data/processed_data/final_acc22_23.mat");
    final_acc=data.final_acc;
    param_name=load_parameter_name();
    
    vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
    
    initial_solution=load_initial_solution();
    ini=[initial_solution;initial_solution(vmaxind)];
    ini=ini(optim_ind);
    
    load('data/optionsPesto.mat')
    optionsPesto.mode      = 'text';
    confLevels =[0.8, 0.9, 1];
    
    %%
    for i=1:length(final_acc)
        acc=final_acc(i);
        filen=strcat('results/parameterization/param_top',char(string(topX)),'/MCMCresCI_',char(acc),'.mat');
        MCMCres=load(filen);
        parameters=MCMCres.parameters;
        parameters.number=length(optim_ind);
        
        parameters = getParameterConfidenceIntervals(parameters, confLevels, optionsPesto);
        
        parameters.CI.S=10.^parameters.CI.S.*ini;
        
        CI=parameters.CI;
        filen=strcat('results/confidence_interval/param_top',char(string(topX)),'_CI300_',char(acc),'.mat');
        save(filen,"CI")
    end
end







