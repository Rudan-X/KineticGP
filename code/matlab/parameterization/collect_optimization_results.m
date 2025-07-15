userpath='C:/Users/Rudan/Documents/GitHub/';

addpath(genpath(strcat(userpath,'KineticGP/code')))
addpath(genpath(strcat(userpath,'KineticGP/data')))

cd(strcat(userpath,'KineticGP/'))

data=load("data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;

param_name=load_parameter_name();
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
param_name=[param_name;strcat(param_name(vmaxind),"y23")]; 
param_name(vmaxind)=strcat(param_name(vmaxind),"y22");

init_sol=load_initial_solution();
init_sol=[init_sol;init_sol(vmaxind)];

np=length(param_name);

topX=10:10:40;
%%


bestfit=zeros(68,4);
reduced_bestfit=zeros(68,4);
optindex=zeros(68,4);
complete=zeros(68,4);

for x=1:4
    estimated_ind=optimized_var_ind(topX(x));
    np=length(estimated_ind);
    parameters=NaN(68,np);
    deg=56-np;

    for arg_ind=1:68
        filen=strcat('results/parameterization/param_top',char(string(topX(x))),'/MCMCres_',final_acc(arg_ind),'.mat');
        if exist(filen,'file')==2
            MCMCres=load(filen);
            complete(arg_ind,x)=1;
        else
            filen=strcat('results/parameterization/param_top',char(string(topX(x))),'/intermediate_',final_acc(arg_ind),'.mat');
            if exist(filen,'file')==2
                MCMCres=load(filen);
                MCMCres.parameters.S=MCMCres.res;
            else
                parameters(arg_ind,:)=NaN;
            end
        end

        MCMCres.parameters.S.logPost=abs(MCMCres.parameters.S.logPost);

        red_chi2=MCMCres.parameters.S.logPost/deg;
        valid_indices = find(red_chi2 < 1);

        if ~isempty(valid_indices)
            ideal_fval=max(red_chi2(valid_indices));
            index=find(red_chi2==ideal_fval);
            optindex(arg_ind,x)=index(end);
            bestfit(arg_ind,x)=MCMCres.parameters.S.logPost(index(end));
            reduced_bestfit(arg_ind,x)=MCMCres.parameters.S.logPost(index(end))/deg;
        else
            ind=find(MCMCres.parameters.S.logPost==min(MCMCres.parameters.S.logPost));
            optindex(arg_ind,x)=ind(end);
            bestfit(arg_ind,x)=min(MCMCres.parameters.S.logPost);
            reduced_bestfit(arg_ind,x)=min(MCMCres.parameters.S.logPost)/deg;
        end

        ratio=10.^MCMCres.parameters.S.par(:,optindex(arg_ind,x));
        optimized_S=ratio.*init_sol(estimated_ind);
        parameters(arg_ind,:)=optimized_S;
        optindex(arg_ind,x)=optindex(arg_ind,x)+100;
        
    end
    data=array2table(parameters,"VariableNames",param_name(estimated_ind)',"RowNames",final_acc);
    writetable(data,strcat("results/parameterization/param_top",string(topX(x)),"/optimized_parameters.csv"),"WriteVariableNames",true,"WriteRowNames",true)
end
%%


writetable(array2table(optindex(:,1)),strcat('results/parameterization/param_top',char(string(topX(x))),'/Optimal_index_burnin.xlsx'),'WriteVariableNames',false);
%% Store all genotypes per parameter


estim_params=struct();
for p=1:11
    estim_params(p).parameters=NaN(68,1000);
end

for x=1:1
    estimated_ind=optimized_var_ind(topX(x));
    np=length(estimated_ind);
    parameters=NaN(68,np);
    deg=56-np;

    for arg_ind=1:68
        filen=strcat('results/parameterization/param_top',char(string(topX(x))),'/MCMCres_',final_acc(arg_ind),'.mat');
        if exist(filen,'file')==2
            MCMCres=load(filen);
            complete(arg_ind,x)=1;
        else
            filen=strcat('results/parameterization/param_top',char(string(topX(x))),'/intermediate_',final_acc(arg_ind),'.mat');
            if exist(filen,'file')==2
                MCMCres=load(filen);
                MCMCres.parameters.S=MCMCres.res;
            else
                parameters(arg_ind,:)=NaN;
            end
        end

        MCMCres.parameters.S.logPost=abs(MCMCres.parameters.S.logPost);

        % Discard the first 100 iterations
        % MCMCres.parameters.S.logPost=MCMCres.parameters.S.logPost(101:end);
        % MCMCres.parameters.S.par=MCMCres.parameters.S.par(:,101:end);
        % 
        % 
        ratio=10.^MCMCres.parameters.S.par;
        optimized_S=ratio.*init_sol(estimated_ind);

        for p=1:11
            estim_params(p).parameters(arg_ind,:)=optimized_S(p,:);
        end
    end
end

% for p=1:11
%     tab=array2table(estim_params(p).parameters);
%     writetable(tab,"results/supp_tables/Sampled_parameters.xlsx",'Sheet',param_name(estimated_ind(p)),'WriteVariableNames',false);
% end


for p=1:11
    tab=array2table(estim_params(p).parameters);
    writetable(tab,strcat("results/supp_tables/Sampled_parameters_",param_name(estimated_ind(p)),".xlsx"),'WriteVariableNames',true);
end


%%
sum_fit=sum(reduced_bestfit,2);
% best
[sorted,ind] = sort(sum_fit,'ascend');
final_acc(ind(1:5))

%worst
[sorted,ind] = sort(sum_fit,'descend');
final_acc(ind(1:5))

sum_fit(ind(1:5))
%%
thre=zeros(4,1);
degs=zeros(4,1);
nvars=zeros(4,1);
for x=1:4
    estim_ind=optimized_var_ind(topX(x));

    % thre(x)=chi2inv(0.05,56-length(estim_ind));
    degs(x)=56-length(estim_ind);
    nvars(x)=length(estim_ind);
end


