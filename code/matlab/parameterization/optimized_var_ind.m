function estimated_ind=optimized_var_ind(topX)
param_name=load_parameter_name();

vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));

ranked_params=readtable("results/sensitivity_results/ranked_parameters.csv");
ranked_params=ranked_params{1:topX,"Ranked_parameters"};
ind=find(contains(ranked_params,"Vm"));
Kms=string(ranked_params(setdiff(1:length(ranked_params),ind)));
Vmaxs=ranked_params(ind);
Vmaxs=unique(replace(string(Vmaxs),'y23',''));

Vmaxs=strcat(repmat(Vmaxs,2,1),[repmat("y22"',length(Vmaxs),1);repmat("y23"',length(Vmaxs),1)]);

final_params=[Kms;Vmaxs];

param_name=[param_name;strcat(param_name(vmaxind),"y23")]; 
param_name(vmaxind)=strcat(param_name(vmaxind),"y22");

[~,estimated_ind]=ismember(final_params,param_name);
