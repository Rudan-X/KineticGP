function simulate_individual_params(k_i)
    userpath='/work/xu2/';
    addpath(genpath(strcat(userpath,'KineticGP/code')))
    addpath(genpath(strcat(userpath,'KineticGP/data')))
    
    cd(strcat(userpath,'KineticGP/'))

    %% Load data
    load("data/processed_data/measured_curves22_23.mat");
    KE_type="equilibrator";
    
    ranked_params=readtable("results/sensitivity_results/ranked_parameters_median.csv");

    top_params_old=ranked_params{1:40,"Ranked_parameters"};

    [~,ind]=ismember(["X32","Y32","F32","Q32","D32","E33","G34"],ranked_params{:,"Ranked_parameters"});
    ranked_params(ind,:)=[];
    
    top_params=ranked_params{1:40,"Ranked_parameters"};
    top_params=setdiff(top_params,top_params_old);
    optimized_param=readtable("results/sensitivity_results/fitted_parameters.csv");
    param_name=optimized_param.Properties.VariableNames;
    param_name(1)=[];
    optimized_param=table2array(optimized_param(:,2:end));

    [~,estimated_ind]=ismember(top_params,param_name);

    vmaxind=find(contains(param_name,'Vm1'),1):find(contains(param_name,'Vm35_Hep'));
    nvar=length(param_name);
    
    init_sol=load_initial_solution();
    np=length(init_sol);
    init_sol=[init_sol;init_sol(vmaxind)];
    
    
    recon=repmat(init_sol',size(optimized_param,1),1);
    recon(:,estimated_ind(k_i))=optimized_param(:,estimated_ind(k_i));
    optimized_param=recon;
    deg=56-1;
    %%
    
    
    display("Simulation started")
    for argind=1:68
        
        parameters22=optimized_param(argind,1:np);
        parameters23=optimized_param(argind,1:np);
        parameters23(vmaxind)=optimized_param(argind,(np+1):end);
    
        [simA22,simgs22]=simulate_ACI(parameters22,mean(ACa22.Tair(:,argind)),KE_type);
        [simA23,simgs23]=simulate_ACI(parameters23,mean(ACa23.Tair(:,argind)),KE_type);
        [simAQ22]=simulate_AQ(parameters22,25,KE_type);
        [simAQ23]=simulate_AQ(parameters23,25,KE_type);
    
        
        zACI22=sum((ACa22.meas(:,argind)-simA22).^2./ACa22.sd(:,argind).^2,'omitnan');
        zgs22=sum((GsCa22.meas(:,argind)-simgs22).^2./GsCa22.sd(:,argind).^2,'omitnan');
        
        zACI23=sum((ACa23.meas(:,argind)-simA23).^2./ACa23.sd(:,argind).^2,'omitnan');
        zgs23=sum((GsCa23.meas(:,argind)-simgs23).^2./GsCa23.sd(:,argind).^2,'omitnan');
        
        zAQ22=sum((AQ22.meas(:,argind)-simAQ22).^2./AQ22.sd(:,argind).^2,'omitnan');
        zAQ23=sum((AQ23.meas(:,argind)-simAQ23).^2./AQ23.sd(:,argind).^2,'omitnan');
    
    
        ACa22.sim(:,argind)=simA22;
        ACa22.z(argind)=zACI22;
        GsCa22.sim(:,argind)=simgs22;
        GsCa22.z(argind)=zgs22;
    
        ACa23.sim(:,argind)=simA23;
        ACa23.z(argind)=zACI23;
        GsCa23.sim(:,argind)=simgs23;
        GsCa23.z(argind)=zgs23;
    
    
        AQ22.sim(:,argind)=simAQ22;
        AQ22.z(argind)=zAQ22;
    
        AQ23.sim(:,argind)=simAQ23;
        AQ23.z(argind)=zAQ23;
    
        temp=round(zACI22+zgs22+zACI23+zgs23+zAQ22+zAQ23,2);
        fprintf("Simulated accession index: %d, reduced chi2: %4.2f \n",argind,temp/deg)
    
    end
    
    %% Store simulation in matlab .mat
     
    
    filen=strcat("results/sensitivity_results/fitted_A/optim_simulation_",param_name(estimated_ind(k_i)),".mat");
    save(filen,'ACa22','GsCa22','AQ22','ACa23','GsCa23','AQ23')
    
    fprintf("Simulation completed, data stored as %s \n", filen)
end

