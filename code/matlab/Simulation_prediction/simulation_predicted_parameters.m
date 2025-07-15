function simulation_predicted_parameters(topX,method,scale,env)
param_name=load_parameter_name();
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
init_sol=load_initial_solution();
np=length(init_sol);
init_sol=[init_sol;init_sol(vmaxind)];
estimated_ind=optimized_var_ind(topX);

%%

fieldcond=readtable("data/processed_data/Testing_Asat3years_fieldcond_accession.csv",'Delimiter',',','VariableNamingRule','preserve');
KE_type="equilibrator";
% env="field";
% method="rrBLUP";
% scale="original";


%% Load predicted parameters with BLUPs of Vmaxs

predicted_params21=readtable(strcat("results/KineticGP_Kprediction/predictedK_BLUP2021_top",string(topX),"_",scale,"_predicted_parameters_",method,".csv"),"TreatAsMissing","NA");
predicted_lines21=predicted_params21{:,1};
predicted_params21=predicted_params21(:,2:end);
predicted_params21=table2array(predicted_params21);

recon=repmat(init_sol(1:np)',size(predicted_params21,1),1);
recon(:,estimated_ind(1:size(predicted_params21,2)))=predicted_params21;
predicted_params21=recon;

AQ=readtable("data/processed_data/Testing_Asat21_accession.csv");
testing21=AQ{:,"Accession"};
%% Load predicted parameters from 2022 and 2023

predicted_params=readtable(strcat("results/KineticGP_Kprediction/predictedK202223_top",string(topX),"_",scale,"_predicted_parameters_",method,".csv"));

predicted_lines22=predicted_params{:,1};
predicted_lines23=predicted_lines22;
predicted_params=predicted_params(:,2:end);
predicted_params=table2array(predicted_params);

recon=repmat(init_sol',size(predicted_params,1),1);
recon(:,estimated_ind)=predicted_params;
predicted_params=recon;

%%
AQ=readtable("data/processed_data/Testing_AQcurves_years22&23_accession.csv");
testing22=AQ{AQ{:,"Year"}==2022,"Accession"};
testing23=AQ{AQ{:,"Year"}==2023,"Accession"};
if env=="field"
    simA21=zeros(length(testing21),7);
    simA22=zeros(length(testing22),7);
    simA23=zeros(length(testing23),7);
    sim_var=["Year";"Accession";strcat("PAR_",string([301,1800,1100,500,300,150,50]))'];
elseif env=="control"
    simA21=zeros(length(testing21),6);
    simA22=zeros(length(testing22),6);
    simA23=zeros(length(testing23),6);
    sim_var=["Year";"Accession";strcat("PAR_",string([1800,1100,500,300,150,50]))'];
end
%%

for k=1:length(testing21)
    fprintf("Accession %d\n",k)
    [~,ind]=ismember(testing21(k),predicted_lines21);

    if sum(isnan(predicted_params21(ind,:)))<1
         parameters21=predicted_params21(ind,:);
        if env=="field"
            indT=find(strcmp(fieldcond{:,"Accession"},testing21(k)));
            indT=indT(fieldcond{indT,"Year"}==2021);
            if isempty(indT)
                simA21(k,:)=NaN;
            else
                Tfield=fieldcond{indT,"meanTemperature"};
                simA21(k,:)=simulate_AQ_field(parameters21,Tfield,KE_type);
            end
        elseif env=="control"
            simA21(k,:)=simulate_AQ(parameters21,25,KE_type);
        end
    else
        simA21(k,:)=NaN;
    end
end

for k=1:length(testing22)
    fprintf("Accession %d\n",k)
    [~,ind]=ismember(testing22(k),predicted_lines22);

    if sum(isnan(predicted_params(ind,:)))<1
         parameters22=predicted_params(ind,1:np);
        if env=="field"
            indT=find(strcmp(fieldcond{:,"Accession"},testing22(k)));
            indT=indT(fieldcond{indT,"Year"}==2022);
            if isempty(indT)
                simA22(k,:)=NaN;
            else
                Tfield=fieldcond{indT,"meanTemperature"};
                simA22(k,:)=simulate_AQ_field(parameters22,Tfield,KE_type);
            end
        elseif env=="control"
            simA22(k,:)=simulate_AQ(parameters22,25,KE_type);
        end
    else
        simA22(k,:)=NaN;
    end
end

for k=1:length(testing23)
    fprintf("Accession %d\n",k)
    [~,ind]=ismember(testing23(k),predicted_lines23);

    if sum(isnan(predicted_params(ind,:)))<1
        parameters23=predicted_params(ind,1:np);
        parameters23(vmaxind)=predicted_params(ind,(np+1):end);
        if env=="field"
            indT=find(strcmp(fieldcond{:,"Accession"},testing23(k)));
            indT=indT(fieldcond{indT,"Year"}==2023);
            if isempty(indT)
                simA23(k,:)=NaN;
            else
                Tfield=fieldcond{indT,"meanTemperature"};
                simA23(k,:)=simulate_AQ_field(parameters23,Tfield,KE_type);
            end
        elseif env=="control"
            simA23(k,:)=simulate_AQ(parameters23,25,KE_type);
        end
    else
        simA23(k,:)=NaN;
    end
end

simA=[simA21;simA22;simA23];
accs=string([testing21;testing22;testing23]);
years=[2021*ones(length(testing21),1);2022*ones(length(testing22),1);2023*ones(length(testing23),1)];


% simA=[simA22;simA23];
% accs=string([testing22;testing23]);
% years=[2022*ones(length(testing22),1);2023*ones(length(testing23),1)];

filen=strcat("results/KineticGP_Asimulation/simulatedA_top",string(topX),"_",scale,"_",method,"_",env,".csv");
writetable(array2table([years,accs,simA],"VariableNames",sim_var),filen,'WriteVariableNames',true,'WriteRowNames',true);