%% Set path
userpath='C:/Users/Rudan/Documents/GitHub/KineticGP/';

addpath(genpath(strcat(userpath,'code/matlab/')))


filen=strcat("results/parameterization/param_top10/optim_simulation.mat");
load(filen);

[~,ind]=sort(ACa22.x(:,1),'ascend');
[~,ind2]=sort(AQ22.x(:,1),'ascend');

%%
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(2, 1, "TileSpacing", "tight");

nexttile
plot(ACa22.x(ind,1),ACa22.meas(ind,:),'-o')
xlim([0,1300])
ax = gca;
ax.FontSize = 14; 

nexttile
[~,ind]=sort(ACa23.x(:,1),'ascend');
plot(ACa23.x(ind,1),ACa23.meas(ind,:),'-o')
xlim([0,1300])
ax = gca;
ax.FontSize = 14; 

xlabel(t,"C_a (µmol/mol)",'FontSize', 16,'FontWeight','bold');
ylabel(t,"Photosynthetic rate (µmol/m^2/s)",'FontSize', 16,'FontWeight','bold');
title(t,"A-C_a curve",'FontSize', 18,'FontWeight','bold')

filen="figures/tempFigures/SF1_measACa22and23.png";
set(gcf, 'PaperPosition', [0 0 4.5 8]); %
saveas(gcf,filen)
%%
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(2, 1, "TileSpacing", "tight");
nexttile
plot(GsCa22.x(ind,1),GsCa22.meas(ind,:),'-o')
xlim([0,1300])
ax = gca;
ax.FontSize = 14; 

nexttile
plot(GsCa23.x(ind,1),GsCa23.meas(ind,:),'-o')
xlim([0,1300])
ax = gca;
ax.FontSize = 14; 

xlabel(t,"C_a (µmol/mol)",'FontSize', 16,'FontWeight','bold');
ylabel(t,"Stomatal conductance (mol/m^2/s)",'FontSize', 16,'FontWeight','bold');
title(t,"G_s-C_a curve",'FontSize', 18,'FontWeight','bold')

filen="figures/tempFigures/SF1_measGsCa22and23.png";
set(gcf, 'PaperPosition', [0 0 4.5 8]); %
saveas(gcf,filen)
%%
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(2, 1, "TileSpacing", "tight");
nexttile
plot(AQ22.x(ind2,1),AQ22.meas(ind2,:),'-o')
xlim([0,1850])
ylim([-3,45])
ax = gca;
ax.FontSize = 14; 

nexttile
plot(AQ23.x(ind2,1),AQ23.meas(ind2,:),'-o')
xlim([0,1850])
ylim([-3,45])
ax = gca;
ax.FontSize = 14; 

xlabel(t,"PAR (µmol/m^2/s)",'FontSize', 16,'FontWeight','bold');
ylabel(t,"Photosynthetic rate (µmol/m^2/s)",'FontSize', 16,'FontWeight','bold');
title(t,"A-PAR curve",'FontSize', 18,'FontWeight','bold')

filen="figures/tempFigures/SF1_measAQ22and23.png";
set(gcf, 'PaperPosition', [0 0 4.5 8]); %
saveas(gcf,filen)
