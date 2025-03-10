function simulation_from_predicted_parameters_CVbased_training21(t)
param_type="original";
env="control";
thres=[12.5,10,7.5,6,5,4,3,2,1,0];
thre=thres(t);
KE_type="equilibrator";
method="BGLR";
folder="equilibrator_parameters_1round";
% folder="equilibrator_11parameters";
%%
folderdir='C:/Users/Rudan/Documents/MATLAB_Drive/KineticGP/';
if ~exist(folderdir, 'dir')
    folderdir='/work/xu2/KineticGP/';
end
addpath(strcat(folderdir,'analysis_publication/'))
addpath(strcat(folderdir,'C4_dynamic_model/'))
addpath(strcat(folderdir,'parameterization/'))

init_sol0=load_initial_solution();
nvar=length(init_sol0);

param_name=load_parameter_name();
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));

param11=["Ki57","Kd57","MRd","BBslope","BBintercept","tao_ActRubisco","Vm2","Vm6"];
[~,ind11]=ismember(param11,param_name);


fieldcond=readtable("../data/processed_data/Testing_Asat3years_fieldcond_accession.csv",'Delimiter',',','VariableNamingRule','preserve');

%%
all_traits=readtable(strcat("../results/",folder,"/optimized_parameters_",KE_type,"_BLUP.csv"),'Delimiter',',','VariableNamingRule','preserve');
all_traits(:,1)=[];
training_params=table2array(all_traits)';
np=size(all_traits,2);

measBLUP=readtable("../data/processed_data/Training_Asat21_accession.csv");
predicted_params=readtable(strcat("../results/",folder,"/",param_type,"_trained_parameters_BLUP_",KE_type,"_",method,".csv"),"TreatAsMissing",["NA"]);
predicted_lines21=predicted_params{:,1};

[~,ind]=ismember(measBLUP{:,"Accession"},predicted_lines21);
predicted_params=predicted_params(ind,2:end);
predicted_lines21=predicted_lines21(ind);

predicted_params=table2array(predicted_params);
ind2=find(isnan(predicted_params(1,:)));
if ~isempty(ind2)
    predicted_params(:,ind2)=repmat(mean(training_params(ind2,:),2)',size(predicted_params,1),1);
end

CVvec=zeros(np,1);
for k=1:np
    CVvec(k,1)=std(predicted_params(:,k),"omitnan")/mean(predicted_params(:,k),"omitnan")*100;
end

goodind=find(CVvec>=thre);

new_params=repmat(mean(training_params,2),1,length(predicted_lines21));
new_params(goodind,:)=predicted_params(:,goodind)';

var21=new_params;

testing21=predicted_lines21;

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
        if folder=="equilibrator_11parameters"
            parameters21=init_sol0;
            parameters21(ind11)=var21(:,ind);
        else
            parameters21=var21(:,ind);
        end

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

if method=="rrBLUP"
    ml="";
elseif method=="BGLR"
    ml="BGLR_";
end

if folder=="equilibrator_11parameters"
    filen=strcat("../GenomicPrediction/testing/",ml,param_type,"_11param_",KE_type,"_",env,"_CV",string(thre),"_year21_training.csv");
elseif folder=="equilibrator_parameters_1round"
    filen=strcat("../GenomicPrediction/testing/",ml,param_type,"_",KE_type,"_",env,"_CV",string(thre),"_year21_training.csv");
end

writetable(array2table([years,accs,simA21],"VariableNames",sim_var),filen,'WriteVariableNames',true,'WriteRowNames',true);

