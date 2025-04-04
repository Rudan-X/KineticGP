userpath='C:/Users/Rudan/Documents/GitHub/KineticGP/';

addpath(genpath(strcat(userpath,'code/matlab/')))

data=load("data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;


param_name=load_parameter_name();
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
param_name=[param_name;strcat(param_name(vmaxind),"y23")]; 

init_sol=load_initial_solution();
init_sol=[init_sol;init_sol(vmaxind)];

% y0=NaN(68,1);
% x0=zeros(length(init_sol),1);
% for acci=4:68
%     display(acci)
%     y0(acci)=sensitivity_obj(x0,final_acc(acci),"equilibrator",1:length(init_sol));
% end
% 
% save('results/sensitivity_results/initial_chi2.mat','y0')

load('results/sensitivity_results/initial_chi2.mat')
np=length(param_name);

Ci_deri=NaN(68,np);
Ci_log=NaN(68,np);
dlny=NaN(68,np);
for acci=1:68
    for v=1:length(param_name)
        filen=strcat(userpath,'results/sensitivity_results/MCMCres_',final_acc(acci),'_',char(param_name(v)),'.mat');

        MCMCres=load(filen);
        MCMCres.parameters.S.logPost=abs(MCMCres.parameters.S.logPost);

        ind=find(MCMCres.parameters.S.logPost==min(MCMCres.parameters.S.logPost));
        xf=10.^MCMCres.parameters.S.par(ind(end))*init_sol(v);
        x0=init_sol(v);

        yf=min(MCMCres.parameters.S.logPost);
        dlny(acci,v)=log(yf/y0(acci));

        dlnx=log(xf/x0);
        if abs(dlnx)<1e-10
            dlnx=1e10;
        end
        
        dy=yf-y0(acci);
        dx=xf-x0;
        if abs(dx)<1e-10
            dx=1e10;
        end


        Ci_log(acci,v)=dlny(acci,v)/dlnx;
        Ci_deri(acci,v)=dy/y0(acci)*x0/dx;
    end
end

% since it is a minimization problem, better use absolute values

Ci_log=abs(Ci_log);
Ci_deri=abs(Ci_deri);

data=array2table(Ci_log,"VariableNames",param_name',"RowNames",final_acc);
writetable(data,"results/sensitivity_results/control_coefficient_log.csv","WriteVariableNames",true,"WriteRowNames",true)


maxn=30;
param_sum=nansum(Ci_deri,1);
[sorted,ind] = sort(param_sum,'descend');
topX2=[param_name(ind(1:maxn)),sorted(1:maxn)',init_sol(ind(1:maxn))];

param_sum=nansum(Ci_log,1);
[sorted,ind] = sort(param_sum,'descend');
topX=[param_name(ind(1:maxn)),sorted(1:maxn)',init_sol(ind(1:maxn))];

params=array2table([param_name(ind),sorted'],"VariableNames",["Ranked_parameters","Averaged_coefficient"]);
writetable(params,"results/sensitivity_results/ranked_parameters.csv","WriteVariableNames",true,"WriteRowNames",false)
