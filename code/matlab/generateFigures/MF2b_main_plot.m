%% Set path
userpath='C:/Users/Rudan/Documents/GitHub/';
addpath(genpath(strcat(userpath,'KineticGP/code')))
addpath(genpath(strcat(userpath,'KineticGP/data')))

top=10:10:40;
years=[2022,2023];


tops1=[11, 20, 34, 44];
degs=56-tops1;

colors = [
    245/255, 212/255, 145/255;  % #F5D491
     41/255, 114/255,  59/255;  % #29723B
    166/255,  46/255,  56/255;  % #A62E38
     59/255,  91/255, 165/255   % #3B5BA5
];

pos=zeros(4,4,2);
pos(:,:,1)=[.20 .62 .3 .3; 
            .20 .57 .3 .3; 
            .20 .52 .3 .3; 
            .20 .47 .3 .3];

pos(:,:,2)=[.20 .20 .3 .3;
            .20 .15 .3 .3;
            .20 .10 .3 .3;
            .20 .05 .3 .3;] ;

%%

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(2, 1, "TileSpacing", "tight");
plot_acc=17;

for i = 1:2
    nexttile
    for x=1:4
        topX=top(x);
        estim_ind=optimized_var_ind(topX);
        sim_data=load(strcat("results/parameterization/param_top",string(topX),"/optim_simulation.mat"));
    
        if i==1
           ACa=sim_data.ACa22;
        else
           ACa=sim_data.ACa23;
        end

        if x==1
            simulate_measured=true;
        else
            simulate_measured=false;
        end
        
        MF2b_plot_curves(ACa,"ACI",plot_acc,degs(x),colors(x,:),'*',14,3,simulate_measured,pos(x,:,i))
    end

    if i==1
       % title("A-vs-Ca (2021)",FontWeight='bold')
       xticks([])
    end
    ax=gca;
    ax.FontSize = 16;
end


xlabel(t,"C_a (µmol/mol)",'FontSize', 16,'FontWeight','bold');
ylabel(t,"Photosynthetic rate (µmol/m^2/s)",'FontSize', 18,'FontWeight','bold');

filen="figures/tempFigures/MF2b_ACa22and23.png";
set(gcf, 'PaperPosition', [0 0 4.5 8]); %
saveas(gcf,filen)

%%

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(2, 1, "TileSpacing", "tight");
plot_acc=17;


for i = 1:2
   nexttile
   for x=1:4
        topX=top(x);
        estim_ind=optimized_var_ind(topX);
        sim_data=load(strcat("results/parameterization/param_top",string(topX),"/optim_simulation.mat"));
        if i==1
           GsCa=sim_data.GsCa22;
        else
           GsCa=sim_data.GsCa23;
        end

        if x==1
            simulate_measured=true;
        else
            simulate_measured=false;
        end
        MF2b_plot_curves(GsCa,"ACI",plot_acc,degs(x),colors(x,:),'*',14,3,simulate_measured,pos(x,:,i))

   end
   if i==1
       % title("A-vs-Ca (2021)",FontWeight='bold')
       xticks([])
   end
   ax=gca;
   ax.FontSize = 16;
end

xlabel(t,"C_a (µmol/mol)",'FontSize', 16,'FontWeight','bold');
ylabel(t,"Stomatal conductance (mol/m^2/s)",'FontSize', 18,'FontWeight','bold');


filen="figures/tempFigures/MF2b_GsCa22and23.png";
set(gcf, 'PaperPosition', [0 0 4.5 8]); %
saveas(gcf,filen)

%%

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(2, 1, "TileSpacing", "tight");
plot_acc=17;

for i = 1:2
   nexttile
   for x=1:4
        topX=top(x);
        estim_ind=optimized_var_ind(topX);
        sim_data=load(strcat("results/parameterization/param_top",string(topX),"/optim_simulation.mat"));
        if i==1
           AQ=sim_data.AQ22;
        else
           AQ=sim_data.AQ23;
        end

        if x==1
            simulate_measured=true;
        else
            simulate_measured=false;
        end

        MF2b_plot_curves(AQ,"AQ",plot_acc,degs(x),colors(x,:),'*',14,3,simulate_measured,pos(x,:,i))
   end
   if i==1
       % title("A-vs-Ca (2021)",FontWeight='bold')
       xticks([])
   end
   ax=gca;
   ax.FontSize = 16;
end

xlabel(t,"PAR (µmol/m^2/s)",'FontSize', 16,'FontWeight','bold');
ylabel(t,"Photosynthetic rate (µmol/m^2/s)",'FontSize', 18,'FontWeight','bold');

% leg=legend({"Measurement 2022","Measurement +/- std 2022","Measurement 2023","Measurement +/- std 2023","Wang et al., + updated Ke","Optimized parameters"},'Location','northeastoutside');

filen="figures/tempFigures/MF2b_AQ22and23.png";
set(gcf, 'PaperPosition', [0 0 4.5 8]); %
saveas(gcf,filen)
% 
% 
% %%
% 
% MF2b_plot_curves(years(i),AQ,"AQ",plot_acc,colors(x,:),'*',14,3,simulate_measured,pos(x,:,i))
% %%
% figure
% plot(1:10,1:10)
% leg=legend({"Measurement (purple), Simulation with N=11(yellow), N=21(green), N=34(red),N=46(blue)"});
% 
% filen="figures/tempFigures/MF2b_legend.png";
% set(gcf, 'PaperPosition', [0 0 6 3]); %
% saveas(gcf,filen)
