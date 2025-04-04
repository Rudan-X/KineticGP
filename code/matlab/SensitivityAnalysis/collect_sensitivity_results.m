userpath='C:\Users\Rudan\Documents\GitHub\KineticGP\';

addpath(strcat(userpath,'code\matlab\C4_dynamic_model/'))
addpath(strcat(userpath,'data\'))
addpath(strcat(userpath,'code\matlab\utils/'))

data=load(strcat(userpath,"data\processed_data/final_acc22_23.mat"));
final_acc=data.final_acc;

param_name=load_parameter_name();
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
param_name=[param_name;strcat(param_name(vmaxind),"y23")]; 
param_name(vmaxind)=strcat(param_name(vmaxind),"y22");

init_sol=load_initial_solution();
init_sol=[init_sol;init_sol(vmaxind)];

np=length(param_name);

parameters=NaN(68,np);

for acci=1:68
    for v=1:length(param_name)
        filen=strcat(userpath,'results/sensitivity_results/MCMCres_',final_acc(acci),'_',char(param_name(v)),'.mat');

        MCMCres=load(filen);
        MCMCres.parameters.S.logPost=abs(MCMCres.parameters.S.logPost);
        ind=find(MCMCres.parameters.S.logPost==min(MCMCres.parameters.S.logPost));
        
        ratio=10.^MCMCres.parameters.S.par(ind(end));
        optimized_S=ratio.*init_sol(v);
        parameters(acci,v)=optimized_S;
    end
end

data=array2table(parameters,"VariableNames",param_name',"RowNames",final_acc);
writetable(data,strcat(userpath,"/results/sensitivity_results/fitted_parameters.csv"),"WriteVariableNames",true,"WriteRowNames",true)
