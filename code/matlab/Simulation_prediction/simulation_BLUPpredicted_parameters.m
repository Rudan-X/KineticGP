function simulation_BLUPpredicted_parameters(topX,method,scale)
param_name=load_parameter_name();
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
init_sol=load_initial_solution();
np=length(init_sol);
init_sol=[init_sol;init_sol(vmaxind)];
estimated_ind=optimized_var_ind(topX);

%%

KE_type="equilibrator";

estim_params=readtable(strcat("results/parameterization/param_top",string(topX),"/optimized_BLUPparameters.csv"));
final_acc=estim_params{:,1};
estim_params=table2array(estim_params(:,2:end));
% method="rrBLUP";
% scale="original";

%% Load BLUPs predicted parameters from 2022 and 2023

predicted_params=readtable(strcat("results/KineticGP_Kprediction/predictedK_BLUP202223_top",string(topX),"_",scale,"_predicted_parameters_",method,".csv"));

predicted_lines22=predicted_params{:,1};

predicted_lines23=predicted_lines22;


predicted_params=predicted_params(:,2:end);
npred=size(predicted_params,2);
estimated_ind=estimated_ind(1:npred);
predicted_params=table2array(predicted_params);

recon=repmat(init_sol',size(predicted_params,1),1);
recon(:,estimated_ind)=predicted_params;
predicted_params=recon;

%%
AQ=readtable("data/processed_data/Testing_AQcurves_years22&23_accession.csv");
testing22=AQ{AQ{:,"Year"}==2022,"Accession"};
testing23=AQ{AQ{:,"Year"}==2023,"Accession"};



simA22=zeros(length(testing22),6);

sim_var=["Year";"Accession";strcat("PAR_",string([1800,1100,500,300,150,50]))'];

%%

for k=1:length(testing22)
    fprintf("Accession %d\n",k)
    [~,ind]=ismember(testing22(k),predicted_lines22);

    if sum(isnan(predicted_params(ind,:)))<1
         parameters22=predicted_params(ind,1:np);
         simA22(k,:)=simulate_AQ(parameters22,25,KE_type);
    else
        simA22(k,:)=NaN;
    end
end


simA=simA22;
accs=string([testing22]);
years=repmat("BLUP",length(testing22),1);


% simA=[simA22;simA23];
% accs=string([testing22;testing23]);
% years=[2022*ones(length(testing22),1);2023*ones(length(testing23),1)];

filen=strcat("results/KineticGP_Asimulation/simulatedA_BLUP2223_top",string(topX),"_",scale,"_",method,".csv");
writetable(array2table([years,accs,simA],"VariableNames",sim_var),filen,'WriteVariableNames',true,'WriteRowNames',true);