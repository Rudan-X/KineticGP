%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART3: Store A at saturating light from 2021

[acc68,~,~,~,~]=load_common_accessions22_23();

acc_info21=readtable("../data/curves/2021_maize_reference.csv",'VariableNamingRule','preserve');
acc_info21{:,"Accession"}=replace(acc_info21{:,"Accession"},'_','_00');

all_acc=unique(string(acc_info21{:,"Accession"}));
testing_acc=setdiff(all_acc,acc68);

na=length(testing_acc);
nq=1;

AQ21=struct();
AQ21.meas=zeros(nq,na);


measA21_all=[];
plots21=[];
acc21=[];

% Only for testing accessions
for argind=1: length(testing_acc)
    genotype=testing_acc(argind);
    % ACa21.meas(:,argind)=load_ACIdata21(genotype);
    [A,plots]=load_Asatdata21(genotype);
    if ~isempty(A)
        measA21_all=[measA21_all;A];
        plots21=[plots21;plots];
        acc21=[acc21;repmat(string(genotype),length(plots),1)];
    end
    AQ21.meas(:,argind)=mean(A); 
end

testing_acc=testing_acc(~isnan(AQ21.meas));
AQ21.meas=AQ21.meas(~isnan(AQ21.meas));


AQ21.x=repmat(1800',1,na);
writetable(array2table([repmat(2021,length(testing_acc),1),testing_acc,AQ21.meas'],"VariableNames",["Year","Accession","A_sat"]),"../data/processed_data/Testing_Asat21_accession.csv",'WriteVariableNames',true,"WriteRowNames",false);


sim_var=["Year";"Accession";"Plot_rep";"A_sat"];
writetable(array2table([repmat(2021,length(acc21),1),acc21,plots21,measA21_all],"VariableNames",sim_var),"../data/processed_data/Testing_Asat21_plot.csv",'WriteVariableNames',true,'WriteRowNames',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART3B: Store A at saturating light from 2021 training set

[acc68,~,~,~,~]=load_common_accessions22_23();

na=length(acc68);
nq=1;

AQ21=struct();
AQ21.meas=zeros(nq,na);


measA21_all=[];
plots21=[];
acc21=[];

% Only for testing accessions
for argind=1: length(acc68)
    genotype=acc68(argind);
    [A,plots]=load_Asatdata21(genotype);
    if ~isempty(A)
        measA21_all=[measA21_all;A];
        plots21=[plots21;plots];
        acc21=[acc21;repmat(string(genotype),length(plots),1)];
    end
    AQ21.meas(:,argind)=mean(A); 
end

acc68=acc68(~isnan(AQ21.meas));
AQ21.meas=AQ21.meas(~isnan(AQ21.meas));


AQ21.x=repmat(1800',1,na);
writetable(array2table([repmat(2021,length(acc68),1),acc68,AQ21.meas'],"VariableNames",["Year","Accession","A_sat"]),"../data/processed_data/Training_Asat21_accession.csv",'WriteVariableNames',true,"WriteRowNames",false);


sim_var=["Year";"Accession";"Plot_rep";"A_sat"];
writetable(array2table([repmat(2021,length(acc21),1),acc21,plots21,measA21_all],"VariableNames",sim_var),"../data/processed_data/Training_Asat21_plot.csv",'WriteVariableNames',true,'WriteRowNames',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART4: Store testing ACI from 2021
[acc68,~,~,~,~]=load_common_accessions22_23();

acc_info21=readtable("../data/curves/2021_maize_reference.csv",'VariableNamingRule','preserve');
acc_info21{:,"Accession"}=replace(acc_info21{:,"Accession"},'_','_00');

all_acc=unique(string(acc_info21{:,"Accession"}));
testing_acc=setdiff(all_acc,acc68);

na=length(testing_acc);
nc=11;

plots21=[];
acc21=[];
measA21_all=[];

ACa21=struct();
ACa21.meas=zeros(nc,na);
for argind=1: length(testing_acc)
    genotype=testing_acc(argind);
    [A,plots]=load_ACIdata21(genotype);
    ACa21.meas(:,argind)=mean(A,2); 

    if ~isnan(A(1,:))
        measA21_all=[measA21_all;A'];
        plots21=[plots21;plots];
        acc21=[acc21;repmat(string(genotype),length(plots),1)];
    end
end
ind=find(~isnan(ACa21.meas(1,:)));
testing_acc=testing_acc(ind);
ACa21.meas=ACa21.meas(:,ind);


sim_var=["Year";"Accession";"Plot_rep";strcat("CO2_",string([400, 600, 800, 1000, 1250, 300, 250, 200, 100, 75, 25]'))];
writetable(array2table([repmat(2021,length(acc21),1),acc21,plots21,measA21_all],"VariableNames",sim_var),"../data/processed_data/Testing_ACI21_plot.csv",'WriteVariableNames',true,'WriteRowNames',true);

sim_var=["Year";"Accession";strcat("CO2_",string([400, 600, 800, 1000, 1250, 300, 250, 200, 100, 75, 25]'))];
writetable(array2table([repmat(2021,length(testing_acc),1),testing_acc,ACa21.meas'],"VariableNames",sim_var),"../data/processed_data/Testing_ACI21_accession.csv",'WriteVariableNames',true,"WriteRowNames",false);


%% PART5: Store training ACI from 2021
[acc68,~,~,~,~]=load_common_accessions22_23();

na=length(acc68);
nc=11;

plots21=[];
acc21=[];
measA21_all=[];

ACa21=struct();
ACa21.meas=zeros(nc,na);
for argind=1: length(acc68)
    genotype=acc68(argind);
    [A,plots]=load_ACIdata21(genotype);
    ACa21.meas(:,argind)=mean(A,2); 

    if ~isnan(A(1,:))
        measA21_all=[measA21_all;A'];
        plots21=[plots21;plots];
        acc21=[acc21;repmat(string(genotype),length(plots),1)];
    end
end
ind=find(~isnan(ACa21.meas(1,:)));
acc68=acc68(ind);
ACa21.meas=ACa21.meas(:,ind);


sim_var=["Year";"Accession";"Plot_rep";strcat("CO2_",string([400, 600, 800, 1000, 1250, 300, 250, 200, 100, 75, 25]'))];
writetable(array2table([repmat(2021,length(acc21),1),acc21,plots21,measA21_all],"VariableNames",sim_var),"../data/processed_data/Training_ACI21_plot.csv",'WriteVariableNames',true,'WriteRowNames',true);

sim_var=["Year";"Accession";strcat("CO2_",string([400, 600, 800, 1000, 1250, 300, 250, 200, 100, 75, 25]'))];
writetable(array2table([repmat(2021,length(acc68),1),acc68,ACa21.meas'],"VariableNames",sim_var),"../data/processed_data/Training_ACI21_accession.csv",'WriteVariableNames',true,"WriteRowNames",false);

