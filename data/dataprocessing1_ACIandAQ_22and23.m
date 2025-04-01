%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART1: Store 68 training accessions with both ACI and AQ curves from 2022 and 2023
clear
if exist('/work/xu2/KineticGP/', 'dir')
    folderdir='/work/xu2/KineticGP/';
else
    folderdir='C:/Users/Rudan/Documents/MATLAB_Drive/KineticGP/';
end
% 
addpath(strcat(folderdir,'C4_dynamic_model/'))
addpath(strcat(folderdir,'parameterization/'))
addpath(strcat(folderdir,'data/'))

%%
[final_acc,~,~,~,~]=load_common_accessions22_23();
save('../data/processed_data/final_acc22_23.mat','final_acc')
writetable(array2table(final_acc,"VariableNames","Accession"),strcat("../data/processed_data/Training68genotypes.csv"),'WriteVariableNames',true,'WriteRowNames',false);

param_name=load_parameter_name();
param11=["Ki57","Kd57","MRd","BBslope","BBintercept","tao_ActRubisco","Vm2","Vm6"];
[~,ind11]=ismember(param11,param_name);
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));

init_sol=load_initial_solution();


writetable(array2table([init_sol(ind11);init_sol(ind11(7:8))],"VariableNames","Initial_value"),strcat("../data/processed_data/parameters11Wang.csv"),'WriteVariableNames',true,'WriteRowNames',false);

writetable(array2table([init_sol;init_sol(vmaxind)],"VariableNames","Initial_value"),strcat("../data/processed_data/parametersWang.csv"),'WriteVariableNames',true,'WriteRowNames',false);

%%
na=length(final_acc);
nc=11;
nq=6;
ACa22=struct();
GsCa22=struct();
ACa23=struct();
GsCa23=struct();
AQ22=struct();
AQ23=struct();

ACa22.meas=zeros(nc,na);
ACa22.Tair=zeros(nc,na);
ACa22.sd=zeros(nc,na);
ACa22.samples=zeros(1,na);

GsCa22.meas=zeros(nc,na);
GsCa22.sd=zeros(nc,na);


ACa23.meas=zeros(nc,na);
ACa23.Tair=zeros(nc,na);
ACa23.sd=zeros(nc,na);
ACa23.samples=zeros(1,na);

GsCa23.meas=zeros(nc,na);
GsCa23.sd=zeros(nc,na);


%%

for argind=1: length(final_acc)
    genotype=final_acc(argind);
    [A,gs,Tair,plots]=load_ACIdata22(genotype);
    ACa22.meas(:,argind)=mean(A,2,"omitnan");
    ACa22.sd(:,argind)=std(A,0,2,"omitnan");
    ACa22.samples(:,argind)=size(A,2);
    GsCa22.meas(:,argind)=mean(gs,2,"omitnan");
    GsCa22.sd(:,argind)=std(gs,0,2,"omitnan");
    ACa22.Tair(:,argind)=mean(Tair,2,"omitnan");
    ACa22.Tsd(:,argind)=std(Tair,0,2,"omitnan");

    [A,gs,Tair,plot]=load_ACIdata23(genotype);
    ACa23.meas(:,argind)=mean(A,2,"omitnan");
    ACa23.sd(:,argind)=std(A,0,2,"omitnan");
    ACa23.samples(:,argind)=size(A,2);

    GsCa23.meas(:,argind)=mean(gs,2,"omitnan");
    GsCa23.sd(:,argind)=std(gs,0,2,"omitnan");
    ACa23.Tair(:,argind)=mean(Tair,2,"omitnan");
    ACa23.Tsd(:,argind)=std(Tair,0,2,"omitnan");
end
%%
AQ22.meas=zeros(nq,na);
AQ22.sd=zeros(nq,na);
AQ22.samples=zeros(1,na);

AQ23.meas=zeros(nq,na);
AQ23.sd=zeros(nq,na);
AQ23.samples=zeros(1,na);

measA22_all=[];
measA23_all=[];
plots22=[];
plots23=[];
acc22=[];
acc23=[];
for argind=1: length(final_acc)
    genotype=final_acc(argind);
    [A,plots]=load_AQdata22(genotype,"training");
    AQ22.meas(:,argind)=mean(A,2,"omitnan");
    AQ22.sd(:,argind)=std(A,0,2,"omitnan");
    AQ22.samples(:,argind)=size(A,2);
    measA22_all=[measA22_all;A'];
    plots22=[plots22;plots'];
    acc22=[acc22;repmat(string(genotype),length(plots),1)];

    [A,plots]=load_AQdata23(genotype,"training");
    AQ23.meas(:,argind)=mean(A,2,"omitnan");
    AQ23.sd(:,argind)=std(A,0,2,"omitnan");
    AQ23.samples(:,argind)=size(A,2);
    measA23_all=[measA23_all;A'];
    plots23=[plots23;plots'];
    acc23=[acc23;repmat(string(genotype),length(plots),1)];
end


%%
ACa22.x=repmat([400, 600, 800, 1000, 1250, 300, 250, 200, 100, 75, 25]',1,na);
GsCa22.x=ACa22.x;
ACa23.x=ACa22.x;
GsCa23.x=ACa22.x;

AQ22.x=repmat([1800,1100,500,300,150,50]',1,na);
AQ23.x=AQ22.x;
save("../data/measured_curves22_23.mat",'ACa22','GsCa22','AQ22','ACa23','GsCa23','AQ23');
%%
measA_all=[measA22_all;measA23_all];
plotrep=[plots22;plots23];
accs=[acc22;acc23];
years=[2022*ones(length(plots22),1);2023*ones(length(plots23),1)];

sim_var=["Year";"Accession";"Plot_rep";strcat("PAR_",string([1800,1100,500,300,150,50]))'];

writetable(array2table([years,accs,plotrep,measA_all],"VariableNames",sim_var),"../data/processed_data/Training_AQcurves_years22&23_plot.csv",'WriteVariableNames',true,'WriteRowNames',true);

%%
save("../data/processed_data/measured_curves22_23.mat",'ACa22','GsCa22','AQ22','ACa23','GsCa23','AQ23');


measA_all=[AQ22.meas';AQ23.meas'];
accs=[final_acc;final_acc];
years=[2022*ones(length(final_acc),1);2023*ones(length(final_acc),1)];

sim_var=["Year";"Accession";strcat("PAR_",string([1800,1100,500,300,150,50]))'];

writetable(array2table([years,accs,measA_all],"VariableNames",sim_var),"../data/processed_data/Training_AQcurves_years22&23_accession.csv",'WriteVariableNames',true,'WriteRowNames',true);
