function simulation_training21(topX,method,scale,env)

init_sol=load_initial_solution();
np=length(init_sol);

estimated_ind=optimized_var_ind(topX);
estimated_ind(estimated_ind>np)=[];
%%

fieldcond=readtable("data/processed_data/Testing_Asat3years_fieldcond_accession.csv",'Delimiter',',','VariableNamingRule','preserve');
KE_type="equilibrator";
%%
measBLUP=readtable("data/processed_data/Training_Asat21_accession.csv");

predicted_params=readtable(strcat("results/KineticGP_Kprediction/predictedK_top",string(topX),"_",scale,"_predicted_training_parameters_",method,"_BLUP.csv"),"TreatAsMissing","NA");

predicted_lines21=predicted_params{:,1};

[~,ind]=ismember(measBLUP{:,"Accession"},predicted_lines21);
predicted_params=table2array(predicted_params(:,2:end));
predicted_params=predicted_params(ind,:);
predicted_lines21=predicted_lines21(ind);
testing21=predicted_lines21;

recon=repmat(init_sol',size(predicted_params,1),1);
recon(:,estimated_ind)=predicted_params;
predicted_params=recon;

if env=="field"
    simA21=zeros(length(testing21),7);
    sim_var=["Year";"Accession";strcat("PAR_",string([301,1800,1100,500,300,150,50]))'];

elseif env=="control"
    simA21=zeros(length(testing21),6);
    sim_var=["Year";"Accession";strcat("PAR_",string([1800,1100,500,300,150,50]))'];

end
%%
for k=1:length(testing21)
    fprintf("Accession %d\n",k)
    [~,ind]=ismember(testing21(k),predicted_lines21);

    if ~isnan(predicted_params(ind,1)) % no SNP data for the accessions with NaN in predicted params
        
        parameters21=predicted_params(ind,:);
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

accs=string(testing21);
years=2021*ones(length(testing21),1);

filen=strcat("results/KineticGP_Asimulation/simulatedA_training21_top",string(topX),"_",scale,"_",method,"_",env,".csv");

writetable(array2table([years,accs,simA21],"VariableNames",sim_var),filen,'WriteVariableNames',true,'WriteRowNames',true);

