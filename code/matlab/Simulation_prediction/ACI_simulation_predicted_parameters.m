function ACI_simulation_predicted_parameters(topX,method,scale)
param_name=load_parameter_name();
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
init_sol=load_initial_solution();
np=length(init_sol);
init_sol=[init_sol;init_sol(vmaxind)];
estimated_ind=optimized_var_ind(topX);

%%

KE_type="equilibrator";
env="control";
% method="rrBLUP";
% scale="original";


%% Load predicted parameters with BLUPs of Vmaxs

predicted_params21=readtable(strcat("results/KineticGP_Kprediction/predictedK_ACI_BLUP2021_top",string(topX),"_",scale,"_predicted_parameters_",method,".csv"),"TreatAsMissing","NA");
predicted_lines21=predicted_params21{:,1};
predicted_params21=predicted_params21(:,2:end);
predicted_params21=table2array(predicted_params21);

recon=repmat(init_sol(1:np)',size(predicted_params21,1),1);
recon(:,estimated_ind(1:size(predicted_params21,2)))=predicted_params21;
predicted_params21=recon;

simA21=zeros(length(predicted_lines21),11);


for k=1:length(predicted_lines21)
    fprintf("Accession %d\n",k)
    parameters21=predicted_params21(k,:);
    simA21(k,:)=simulate_ACI(parameters21,25,KE_type);
end
% AQ=readtable("data/processed_data/Testing_Asat21_accession.csv");
% predicted_lines21=AQ{:,"Accession"};
%% Load predicted parameters from 2022 

predicted_params22=readtable(strcat("results/KineticGP_Kprediction/predictedK_ACI2022_top",string(topX),"_",scale,"_predicted_parameters_",method,".csv"));

predicted_lines22=predicted_params22{:,1};

predicted_params22=predicted_params22(:,2:end);
predicted_params22=table2array(predicted_params22(:,1:topX));


recon=repmat(init_sol',size(predicted_params22,1),1);
recon(:,estimated_ind(1:size(predicted_params22,2)))=predicted_params22;
predicted_params22=recon;

% recon=repmat(init_sol(1:np)',size(predicted_params22,1),1);
% recon(:,estimated_ind(1:size(predicted_params22,2)))=predicted_params22;

simA22=zeros(length(predicted_lines22),11);

for k=1:length(predicted_lines22)
    fprintf("Accession %d\n",k)
    parameters22=predicted_params22(k,1:np);
    simA22(k,:)=simulate_ACI(parameters22,25,KE_type);
end

%%
predicted_params23=readtable(strcat("results/KineticGP_Kprediction/predictedK_ACI2023_top",string(topX),"_",scale,"_predicted_parameters_",method,".csv"));
predicted_lines23=predicted_params23{:,1};

predicted_params23=predicted_params23(:,2:end);
predicted_params23=table2array(predicted_params23);

predicted_lines23(isnan(predicted_params23(:,1)))=[];
predicted_params23(isnan(predicted_params23(:,1)),:)=[];

recon=repmat(init_sol',size(predicted_params23,1),1);
recon(:,estimated_ind(1:size(predicted_params23,2)))=predicted_params23;
predicted_params23=recon;

%%

simA23=zeros(length(predicted_lines23),11);
% sim_var=["Year";"Accession";strcat("PAR_",string([1800,1100,500,300,150,50]))'];

for k=1:length(predicted_lines23)
    fprintf("Accession %d\n",k)
    parameters23=predicted_params23(k,1:np);
    parameters23(vmaxind)=predicted_params23(k,(np+1):end);

    simA23(k,:)=simulate_ACI(parameters23,25,KE_type);
end

simA=[simA21;simA22;simA23];
accs=string([predicted_lines21;predicted_lines22;predicted_lines23]);
years=[2021*ones(length(predicted_lines21),1);2022*ones(length(predicted_lines22),1);2023*ones(length(predicted_lines23),1)];

sim_var=["Year";"Accession";strcat("CO2_",string([400, 600, 800, 1000, 1250, 300, 250, 200, 100, 75, 25]))'];
% simA=[simA22;simA23];
% accs=string([predicted_lines22;predicted_lines23]);
% years=[2022*ones(length(predicted_lines22),1);2023*ones(length(predicted_lines23),1)];

filen=strcat("results/KineticGP_Asimulation/simulatedACI_top",string(topX),"_",scale,"_",method,"_",env,".csv");
writetable(array2table([years,accs,simA],"VariableNames",sim_var),filen,'WriteVariableNames',true,'WriteRowNames',true);

