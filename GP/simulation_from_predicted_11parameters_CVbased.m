function simulation_from_predicted_11parameters_CVbased(t)
% Change the following variables:
param_type="original";
env="control";
thres=[12.5,10,7.5,6,5,4,3,2,1,0];
thre=thres(t);
KE_type="equilibrator";
% method="rrBLUP";
method="BGLR";
folder="equilibrator_11parameters";
%%
folderdir='C:/Users/Rudan/Documents/MATLAB_Drive/KineticGP/';
if ~exist(folderdir, 'dir')
    folderdir='/work/xu2/KineticGP/';
end
addpath(strcat(folderdir,'analysis_publication/'))
addpath(strcat(folderdir,'C4_dynamic_model/'))
addpath(strcat(folderdir,'parameterization/'))

init_sol0=load_initial_solution();
nvar=8;

param_name=load_parameter_name();
param11=["Ki57","Kd57","MRd","BBslope","BBintercept","tao_ActRubisco","Vm2","Vm6"];
[~,ind11]=ismember(param11,param_name);


%%

fieldcond=readtable("../data/processed_data/Testing_Asat3years_fieldcond_accession.csv",'Delimiter',',','VariableNamingRule','preserve');

%%
all_traits=readtable(strcat("../results/",folder,"/optimized_parameters_",KE_type,"_BLUP.csv"),'Delimiter',',','VariableNamingRule','preserve');
all_traits(:,1)=[];
training_params=table2array(all_traits)';
np=size(all_traits,2);

predicted_params=readtable(strcat("../results/",folder,"/",param_type,"_predicted_parameters_BLUP_",KE_type,"_",method,".csv"),"TreatAsMissing","NA");
predicted_lines21=predicted_params{:,1};
predicted_params=predicted_params(:,2:end);
predicted_params=table2array(predicted_params);

if method=="BGLR"
    ind2=find(isnan(predicted_params(1,:)));
    if ~isempty(ind2)
        predicted_params(:,ind2)=repmat(mean(training_params(ind2,:),2)',size(predicted_params,1),1);
    end
end

CVvec=zeros(np,1);
for k=1:np
    CVvec(k,1)=std(predicted_params(:,k),"omitnan")/mean(predicted_params(:,k),"omitnan")*100;
end

goodind=find(CVvec>=thre);

new_params=repmat(mean(training_params,2),1,length(predicted_lines21));
new_params(goodind,:)=predicted_params(:,goodind)';

var21=new_params;

AQ=readtable("../data/processed_data/Testing_Asat21_accession.csv");
testing21=AQ{:,"Accession"};

%%
all_traits=readtable(strcat("../results/",folder,"/optimized_parameters_",KE_type,".csv"),'Delimiter',',','VariableNamingRule','preserve');
all_traits(:,1)=[];
training_params=table2array(all_traits)';
np=size(all_traits,2);

predicted_params=readtable(strcat("../results/",folder,"/",param_type,"_predicted_parameters_",KE_type,"_",method,".csv"));

predicted_lines22=predicted_params{:,1};
predicted_lines23=predicted_lines22;
predicted_params=predicted_params(:,2:end);
predicted_params=table2array(predicted_params);


CVvec=zeros(np,1);
for k=1:np
    CVvec(k,1)=std(predicted_params(:,k),"omitnan")/mean(predicted_params(:,k),"omitnan")*100;
end

goodind=find(CVvec>=thre);

new_params=repmat(mean(training_params,2),1,length(predicted_lines22));
new_params(goodind,:)=predicted_params(:,goodind)';

var22=new_params(1:nvar,:);
var23=var22;
var23(7:8,:)=new_params((nvar+1):end,:);


%%


AQ=readtable("../data/processed_data/Testing_AQcurves_years22&23_accession.csv");
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

    if sum(isnan(var21(:,ind)))<1
         parameters21=init_sol0;
         parameters21(ind11)=var21(:,ind);
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

    if sum(isnan(var22(:,ind)))<1
         parameters22=init_sol0;
         parameters22(ind11)=var22(:,ind);
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

    if sum(isnan(var23(:,ind)))<1
        parameters23=init_sol0;
        parameters23(ind11)=var23(:,ind);
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

%%

if method=="rrBLUP"
    ml="";
elseif method=="BGLR"
    ml="BGLR_";
end

simA=[simA21;simA22;simA23];
accs=string([testing21;testing22;testing23]);
years=[2021*ones(length(testing21),1);2022*ones(length(testing22),1);2023*ones(length(testing23),1)];

filen=strcat("../GenomicPrediction/testing/",ml,param_type,"_11param_",KE_type,"_",env,"_CV",string(thre),".csv");
writetable(array2table([years,accs,simA],"VariableNames",sim_var),filen,'WriteVariableNames',true,'WriteRowNames',true);

