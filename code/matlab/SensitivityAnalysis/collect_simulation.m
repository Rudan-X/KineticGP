userpath='C:/Users/Rudan/Documents/GitHub/';

addpath(genpath(strcat(userpath,'KineticGP/code')))
addpath(genpath(strcat(userpath,'KineticGP/data')))

cd(strcat(userpath,'KineticGP/'))

%% Load data
KE_type="equilibrator";
data=load("data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;

ranked_params=readtable("results/sensitivity_results/ranked_parameters_median.csv");
[~,ind]=ismember(["X32","Y32","F32","Q32","D32","E33","G34"],ranked_params{:,"Ranked_parameters"});
ranked_params(ind,:)=[];

top_params=ranked_params{1:40,"Ranked_parameters"};

optimized_param=readtable("results/sensitivity_results/fitted_parameters.csv");
param_name=optimized_param.Properties.VariableNames;
param_name(1)=[];
optimized_param=table2array(optimized_param(:,2:end));

[~,estimated_ind]=ismember(top_params,param_name);


%%
for k_i= 1:40
    filen=strcat("results/sensitivity_results/fitted_A/optim_simulation_",param_name(estimated_ind(k_i)),".mat");
    load(filen)
    mat=[ ACa22.sim',AQ22.sim',ACa23.sim',AQ23.sim'];
    
    sim_var=["Accession"; strcat("2022_CO2_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
        strcat("2022_PAR_",string([1800,1100,500,300,150,50]))';
        strcat("2023_CO2_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
        strcat("2023_PAR_",string([1800,1100,500,300,150,50]))'];
    
    T = array2table([final_acc,mat],"VariableNames",sim_var);
    
    writetable(T,strcat("results/sensitivity_results/fitted_A/fitted_simulation_",param_name(estimated_ind(k_i)),".csv"),'WriteVariableNames',true,'WriteRowNames',true);
end

%%

reduced_chi=zeros(68,40);
for k_i= 1:40
    filen=strcat("results/sensitivity_results/fitted_A/optim_simulation_",param_name(estimated_ind(k_i)),".mat");
    load(filen)

    deg=56-1;
    temp=ACa22.z + GsCa22.z + ACa23.z + GsCa23.z + AQ22.z + AQ23.z;
    reduced_chi(:,k_i)=round(temp',2)/deg;
end

T = array2table(reduced_chi,'VariableNames',top_params);
writetable(T,strcat("results/sensitivity_results/fitted_A/reduced_chi_square_top40.csv"),'WriteVariableNames',true);
