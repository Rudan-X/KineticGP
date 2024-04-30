% function simulate_initial_result()
%% Load data
[init_sol,final_acc,~,~]=optim_initialization_parameters();

%%
ACa21=struct();
GsCa21=struct();
AQ22=struct();
nca=12;
ACa21.sim0=zeros(nca,length(final_acc));
ACa21.meas=zeros(nca,length(final_acc));
ACa21.sd=zeros(nca,length(final_acc));
ACa21.x=zeros(nca,length(final_acc));
ACa21.z0=5e4*ones(length(final_acc),1);

GsCa21.sim0=zeros(nca,length(final_acc));
GsCa21.meas=zeros(nca,length(final_acc));
GsCa21.sd=zeros(nca,length(final_acc));
GsCa21.x=zeros(nca,length(final_acc));
GsCa21.z0=5e4*ones(length(final_acc),1);

nQ=6;
AQ22.sim0=zeros(nQ,length(final_acc));
AQ22.meas=zeros(nQ,length(final_acc));
AQ22.sd=zeros(nQ,length(final_acc));
AQ22.x=zeros(nQ,length(final_acc));
AQ22.z0=5e4*ones(length(final_acc),1);

nvar=length(init_sol);
var=zeros(nvar,1);
ratio=10.^var;
parameters=ratio.*init_sol;

%%
for argind=32: length(final_acc)
    
    display(argind)
    
    
    [res,zA21,zgs21] = simulate_photosynthesis(parameters,final_acc(argind),'ACI');
    if ~isempty(res)
        ACa21.sim0(:,argind)=res(:,1);
        ACa21.z0(argind)=zA21;
        ACa21.meas(:,argind)=res(:,2);
        ACa21.sd(:,argind)=res(:,3);
        ACa21.x(:,argind)=res(:,end);

        GsCa21.sim0(:,argind)=res(:,4);
        GsCa21.z0(argind)=zgs21;
        GsCa21.meas(:,argind)=res(:,5);
        GsCa21.sd(:,argind)=res(:,6);
        GsCa21.x(:,argind)=res(:,end);
    end

    
    [res,zA22,~] = simulate_photosynthesis(parameters,final_acc(argind),'AQ');
    if ~isempty(res)
        AQ22.sim0(:,argind)=res(:,1);
        AQ22.meas(:,argind)=res(:,2);
        AQ22.sd(:,argind)=res(:,3);
        AQ22.x(:,argind)=res(:,4);
        AQ22.z0(argind)=zA22;
    end

end

%%
save("../results/simulations/simulation_initial.mat",'ACa21','GsCa21','AQ22')
