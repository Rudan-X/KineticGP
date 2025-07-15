function SFX2_prepare_figures(t)

top_x=[10,20,30,40];
topX=top_x(t);
%% Load data
KE_type="equilibrator";
param_name=load_parameter_name();
vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
init_sol=load_initial_solution();
np=length(init_sol);
init_sol=[init_sol;init_sol(vmaxind)];


load("data/processed_data/measured_curves22_23.mat");
data=load("data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;

metnames=load_metabolite_names();


%%

estim_params=readtable(strcat("results/parameterization/param_top",string(topX),"/optimized_parameters.csv"));

estim_params=estim_params(:,2:end);
estim_params=table2array(estim_params);

recon=repmat(init_sol',size(estim_params,1),1);

estimated_ind=optimized_var_ind(topX);
recon(:,estimated_ind)=estim_params;
estim_params=recon;
%%
simulated_conc22=struct();
simulated_conc23=struct();

display("Simulation started")
for argind=1:68
    parameters22=estim_params(argind,1:np);
    parameters23=parameters22;
    parameters23(vmaxind)=estim_params(argind,(np+1):end);

    sim_met22=SFX2_simulate_ACI_met(parameters22,mean(ACa22.Tair(:,argind)),KE_type);
    simulated_conc22(argind).met=sim_met22;
    sim_met23=SFX2_simulate_ACI_met(parameters23,mean(ACa23.Tair(:,argind)),KE_type);
    simulated_conc23(argind).met=sim_met23;
end

%%
filen=strcat("results/metabolite_profiles/top",string(topX),"_optim_simulation_metabolites.mat");
save(filen,'simulated_conc22','simulated_conc23')

fprintf("Simulation completed, data stored as %s \n", filen)

%%
filen=strcat(strcat("results/metabolite_profiles/top",string(topX),"_optim_simulation_metabolites.mat"));

load(filen)

mets=struct();
mets(1).names=["PGA[MC]","PGA[Mchl]","PGA[Bchl]"];
mets(length(mets)+1).names=["FBP[MC]","FBP[Bchl]"];
mets(length(mets)+1).names="Starch[Bchl]";
mets(length(mets)+1).names="SUC[MC]";
mets(length(mets)+1).names=["GCA[Bchl]","GCA[Bper]"];

mets(length(mets)+1).names=["Malate[MC]","Malate[Mchl]","Malate[BC]","Malate[Bchl]"];
mets(length(mets)+1).names=["Pyruvate[MC]","Pyruvate[Mchl]","Pyruvate[BC]""Pyruvate[Bchl]"];
mets(length(mets)+1).names=["PEP[MC]","PEP[Mchl]","PEP[BC]","PEP[Bchl]"];

plotmets=["PGA","FBP","Starch","Sucrose","GCA","Malate","Pyruvate","PEP"];

%%
for year=22:23
    metcon=zeros(68,length(mets));
    for m=1:length(mets)
        ind=find(contains(metnames,mets(m).names));
        for argind=1:68
            if year==22
                metcon(argind,m)=sum(simulated_conc22(argind).met(1,ind));
            else
                metcon(argind,m)=sum(simulated_conc23(argind).met(1,ind));
            end
        end
    end
    T=array2table(real(metcon),'VariableNames',plotmets);
    writetable(T,strcat("results/metabolite_profiles/top",string(topX),"_metabolites_across_genotypes",string(year),".csv"),'WriteVariableNames',true)
end
