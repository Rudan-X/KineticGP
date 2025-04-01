userpath='C:\Users\Rudan\Documents\MATLAB_Drive\';

addpath(strcat(userpath,'KineticGP/C4_dynamic_model/'))
addpath(strcat(userpath,'KineticGP/parameterization/'))

data=load("../data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;

param_name=load_parameter_name();
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
param_name=[param_name;strcat(param_name(vmaxind),"y23")]; 

np=length(param_name);

gainfit=NaN(68,np);

for acci=52:68
    for v=1:length(param_name)
        filen=strcat('../SensitivityAnalysis/sensitivity_results/MCMCres_',final_acc(acci),'_',char(param_name(v)),'.mat');
    
        if exist(filen,'file')==2
            MCMCres=load(filen);
            MCMCres.parameters.S.logPost=abs(MCMCres.parameters.S.logPost);
            diff=MCMCres.parameters.S.logPost(1)-min(MCMCres.parameters.S.logPost);
            gainfit(acci,v)=diff/MCMCres.parameters.S.logPost(1)*100;
        else
            filen=strcat('../SensitivityAnalysis/sensitivity_results/intermediate_',final_acc(acci),char(param_name(v)),'.mat');
            if exist(filen,'file')==2
                MCMCres=load(filen);
                MCMCres.parameters.S=MCMCres.res;
                MCMCres.parameters.S.logPost=abs(MCMCres.parameters.S.logPost);
                diff=MCMCres.parameters.S.logPost(1)-min(MCMCres.parameters.S.logPost);
                gainfit(acci,v)=diff/MCMCres.parameters.S.logPost(1)*100;
            else
                gainfit(acci,v)=-10000;
            end
        end
        
    end
end


param_sum=nansum(gainfit,1)/2;
[sorted_gain,ind] = sort(param_sum,'descend');
topX=[param_name(ind(1:30)),sorted_gain(1:30)'];

[sorted_gain,ind] = sort(param_sum,'ascend');
worseX=[param_name(ind(1:100)),sorted_gain(1:100)'];


param11=["Ki57","Kd57","MRd","BBslope","BBintercept","tao_ActRubisco","Vm2","Vm6"];
[~,ind11]=ismember(param11,param_name);
check=[param_name(ind11),param_sum(ind11)'];

subplot(1,2,1)
plot(MCMCres.parameters.S.logPost,'.')
ylabel("posterior")

subplot(1,2,2)
plot(MCMCres.parameters.S.par,'.')
ylabel("Parameter value")


bestfitACCU=bestfit;
maxiterACCU=maxiter;


bestfit=zeros(68,1);
maxiter=zeros(68,1);
for acci=1:68
    [~,~,iter,fval]=load_result11(final_acc(acci),"equilibrator_11parameters");
    bestfit(acci)=fval;
    maxiter(acci)=iter;
end

