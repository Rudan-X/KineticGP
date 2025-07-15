userpath='C:/Users/Rudan/Documents/GitHub/KineticGP/';

addpath(strcat(userpath,'code/matlab/C4_dynamic_model/'))
addpath(strcat(userpath,'data/'))
addpath(strcat(userpath,'code/matlab/utils/'))

data=load(strcat(userpath,"data/processed_data/final_acc22_23.mat"));
final_acc=data.final_acc;

param_name=load_parameter_name();
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
param_name=[param_name;strcat(param_name(vmaxind),"y23")]; 


init_sol=load_initial_solution();
init_sol=[init_sol;init_sol(vmaxind)];

np=length(param_name);
%%

for acci=1:68
    parameters=NaN(500,np);
    parameters_log=NaN(500,np);
    chi2=NaN(500,np);

    for v=1:length(param_name)
        filen=strcat(userpath,'results/sensitivity_results/MCMCres_',final_acc(acci),'_',char(param_name(v)),'.mat');

        MCMCres=load(filen);
        MCMCres.parameters.S.logPost=abs(MCMCres.parameters.S.logPost);

        ratio=10.^MCMCres.parameters.S.par;
        optimized_S=ratio.*init_sol(v);
        parameters(:,v)=optimized_S;
        parameters_log(:,v)=MCMCres.parameters.S.par;

        chi2(:,v)=MCMCres.parameters.S.logPost;
    end

    data=array2table(parameters,"VariableNames",param_name');
    writetable(data,strcat(userpath,"/results/sensitivity_results/MCMCsampled_parameters_",final_acc(acci),".xlsx"), "Sheet", "parameters","WriteVariableNames",true,"WriteRowNames",false)
    
    data=array2table(chi2,"VariableNames",param_name');
    writetable(data,strcat(userpath,"/results/sensitivity_results/MCMCsampled_parameters_",final_acc(acci),".xlsx"), "Sheet", "chi-square","WriteVariableNames",true,"WriteRowNames",false)
    
    data=array2table(parameters_log,"VariableNames",param_name');
    writetable(data,strcat(userpath,"/results/sensitivity_results/MCMCsampled_parameters_log_",final_acc(acci),".xlsx"), "Sheet", "parameters","WriteVariableNames",true,"WriteRowNames",false)
end
