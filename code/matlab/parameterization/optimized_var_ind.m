function estimated_ind=optimized_var_ind(topX)
param_name=load_parameter_name();

vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));

ranked_params=readtable("results/sensitivity_results/ranked_parameters_median.csv");
ranked_params=ranked_params{:,"Ranked_parameters"};

excluded=["X32";"Y32";"F32";"Q32";"D32";"E33";"G34"];
[~,ind0]=ismember(excluded,ranked_params);
ranked_params(ind0)=[];
ranked_params=ranked_params(1:topX);
ind1=find(contains(ranked_params,"Vm"));
ind2=find(contains(ranked_params,"Jmax"));
ind=sort([ind1;ind2]);
Kms=string(ranked_params(setdiff(1:length(ranked_params),ind)));
Vmaxs=ranked_params(ind);
Vmaxs=unique(replace(string(Vmaxs),'y23',''));

Vmaxs=strcat(repmat(Vmaxs,2,1),[repmat(""',length(Vmaxs),1);repmat("y23"',length(Vmaxs),1)]);

final_params=[Kms;Vmaxs];

param_name=[param_name;strcat(param_name(vmaxind),"y23")]; 
[~,estimated_ind]=ismember(final_params,param_name);
estimated_ind=sort(estimated_ind);