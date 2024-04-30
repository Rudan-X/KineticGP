%%


[~,final_acc,param_name,vmaxind]=optim_initialization_parameters();

for i=1:10%length(final_acc)
    accession=final_acc(i);
    load(strcat("../results/optimized_parameters/history_MCMCres_",accession,".mat"));
    ind=find(isnan(res.par(1,:)),1);
    Variable=res.par(:,(ind-1));
    [init_sol,~,~,vmaxind]=optim_initialization_parameters();

    nvar=length(init_sol);
    
    parameters21=Variable(1:nvar);
    parameters22=Variable(1:nvar);
    parameters22(vmaxind)=Variable((nvar+1):end);
    
    % Transform the variable from logarithmic to original scale
    ratio=10.^parameters21;
    parameters21=ratio.*init_sol;
    ratio=10.^parameters22;
    parameters22=ratio.*init_sol;

    Variable=parameters21;
    
    curve='ACI';
    envFactor=struct();

    if strcmp(curve,'ACI')
        % Load ACI curves
        [Ca_t,Ci,A_t,A_t_sd,gs_t,gs_t_sd]=load_ACIdata(accession);
        air_CO2_t0=Ca_t(1);
        Ci_t0=Ci(1);
        
        % envFactor.t_inter_ca=t_ca;
        envFactor.Ca_t=Ca_t;
        % Set empty vector for PAR changes
        envFactor.Q_t=[];
        % set constant PAR:
        PAR_t0=1800;
    
    elseif strcmp(curve,'AQ')
        % Load AQ curves
        [A_t,A_t_sd,Q_t,Ci]=load_AQdata(accession);
    
        envFactor.Q_t=Q_t;
        % Set empty vector for Ca changes
        envFactor.Ca_t=[];
        % set constant Air_CO2:
        air_CO2_t0=400;
        
        global Ci_t
        Ci_t=Ci;
        
        if Ci_t(1)>0
            Ci_t0=Ci_t(1);
        else
            aveCI=63.12;
            Ci_t0=aveCI;
        end
    
        PAR_t0=Q_t(1);
    end
    
    % Set experimental condition and load global variables for C4 model
    air_RH=0.65;
    air_temp=25;
    
    optim_initialization_global_env_variables(air_temp,air_RH,air_CO2_t0,Ci_t0,PAR_t0) 
    
    
    % Load metabolite names
    mets_name=load_metabolite_names();
    [~,ind_gs]=ismember("Gs[Leaf]",mets_name);
    
    
    % Rearrange the variables into kintic parameters
    params=Variable(1:6);
    
    global KVlen
    pend=6;
    km1=pend+1;
    kmend=pend+sum(KVlen);
    vm1=kmend+1;
    vmend=kmend+53; %(vmax=53)
    act1=vmend+1;
    actend=vmend+10;
    perm1=actend+1;
    
    kvalues=Variable(km1:kmend);
    vmaxs=Variable(vm1:vmend);
    act_rate=Variable(act1:actend);
    permeab=Variable(perm1:end);
    
    [Ini,max_vel,kms]=RAC4leafMetaIni(params,kvalues,vmaxs);% Initial values
    
    % OBJECTIVE FUNCTION: Get simulated A and gsw
    % Set ODE options
    nm=length(Ini);
    options=odeset('NonNegative',1:nm, 'RelTol', 1e-04,'Events', @reaching_steadyA); %  
    
    tspan=[0 60*60]; % maximum simulation time, 60minutes
    
    % Initialize simulation vectors
    Gs_VEL_tot=[];
    simA=zeros(length(A_t),1);
    simGs=zeros(length(A_t),1);
    
    % global Gs_VEL initialized in optim_initialization...
    if strcmp(curve,'ACI')
        nchange=length(Ca_t);
    elseif strcmp(curve,'AQ')
        nchange=length(Q_t);
    end
    
    for i=1:nchange
        % ode_sol=[];
        initialize_reaction_rates();
        global Gs_VEL
        if i==1
            xt0=Ini;
        else
            xt0=ode_sol.y(:,end);
        end
        if strcmp(curve,'ACI')
            envFactor.Ca_t=Ca_t(i);
        elseif strcmp(curve,'AQ')
            envFactor.Q_t=Q_t(i);
        end
        ode_sol=ode15s(@(t,x)RAC4leafMetaMB(t,x,params,max_vel,kms,act_rate,permeab,envFactor),tspan,xt0,options); 
        
        if ~isempty(ode_sol) && max(ode_sol.x)~=tspan(2)
            ind_A=2;
            simA(i)=Gs_VEL(end,ind_A);
        else
            simA(i)=1e10;
        end
    
        if i==1
            Gs_VEL_tot=Gs_VEL;
        else
            Gs_VEL(:,1)=Gs_VEL_tot(end,1)+Gs_VEL(:,1);
            Gs_VEL_tot=[Gs_VEL_tot;Gs_VEL];
        end
    
        simGs(i)=ode_sol.y(ind_gs,end);
    end
    figure
    plot(Gs_VEL_tot(200:end,1),Gs_VEL_tot(200:end,2))
end

%%
