%% Set path
if exist('/work/xu2/KineticGP/', 'dir')
    folderdir='/work/xu2/KineticGP/';
else
    folderdir='C:/Users/Rudan/Documents/MATLAB_Drive/KineticGP/';
end
% 
addpath(strcat(folderdir,'C4_dynamic_model/'))
addpath(strcat(folderdir,'parameterization/'))
addpath(strcat(folderdir,'generateMainFig/'))
addpath(strcat(folderdir,'generateSupFig/'))



% MF2b_SF2_simulate_initial_result("original","original_parameters_1round")
% 
% 
% %%
% filen="results/parameterization/default/optim_simulation_with_initial_original.mat";
% load(filen)

profiles=readtable("results\parameterization\default\default_single_simulation.csv");
curves=table2array(profiles(:,2:end));
co2=[400, 600, 800, 1000, 1250,300,250,200, 100,75, 25];
par=[1800,1100,500,300,150,50];
%%

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(1, 1, "TileSpacing", "tight");

nexttile
plot(co2,curves(1,1:11),'square','MarkerSize',10,'LineWidth',2)
xlim([0,1300])
ax = gca;
ax.FontSize = 14; 


xlabel(t,"C_a (µmol/mol)",'FontSize', 16,'FontWeight','bold');
ylabel(t,"Photosynthetic rate (µmol/m^2/s)",'FontSize', 16,'FontWeight','bold');
title(t,"A-C_a curve",'FontSize', 18,'FontWeight','bold')

filen="figures/tempFigures/SF2_initial_simACa22and23_original.png";
set(gcf, 'PaperPosition', [0 0 4.5 5]); %
saveas(gcf,filen)
%%
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(1, 1, "TileSpacing", "tight");
nexttile
plot(co2,curves(1,12:22),'square','MarkerSize',10,'LineWidth',2)
xlim([0,1300])
ax = gca;
ax.FontSize = 14; 

xlabel(t,"C_a (µmol/mol)",'FontSize', 16,'FontWeight','bold');
ylabel(t,"Stomatal conductance (mol/m^2/s)",'FontSize', 16,'FontWeight','bold');
title(t,"G_s-C_a curve",'FontSize', 18,'FontWeight','bold')

filen="figures/tempFigures/SF2_initial_simGsCa22and23_original.png";
set(gcf, 'PaperPosition', [0 0 4.5 5]); %
saveas(gcf,filen)
%%

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(1, 1, "TileSpacing", "tight");
nexttile
plot(par,curves(1,23:end),'square','MarkerSize',10,'LineWidth',2)
xlim([0,1850])
ax = gca;
ax.FontSize = 14; 

xlabel(t,"PAR (µmol/m^2/s)",'FontSize', 16,'FontWeight','bold');
ylabel(t,"Photosynthetic rate (µmol/m^2/s)",'FontSize', 16,'FontWeight','bold');
title(t,"A-PAR curve",'FontSize', 18,'FontWeight','bold')

filen="figures/tempFigures/SF2_initial_simAQ22and23_original.png";
set(gcf, 'PaperPosition', [0 0 4.5 5]); %
saveas(gcf,filen)

%% Write into excel file

data=load("data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;

mat=[ ACa22.sim0',ACa23.sim0',GsCa22.sim0',GsCa23.sim0',AQ23.sim0',AQ22.sim0'];


sim_var=["Accession"; 
    strcat("ACa2022_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    strcat("ACa2023_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    strcat("GsCa2022_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    strcat("GsCa2023_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    strcat("AQ2022_",string([1800,1100,500,300,150,50]))';
    strcat("AQ2023_",string([1800,1100,500,300,150,50]))'];

T = array2table([final_acc,mat],"VariableNames",sim_var);

writetable(T,"results/parameterization/default/default_simulation.csv",'WriteVariableNames',true,'WriteRowNames',true);

