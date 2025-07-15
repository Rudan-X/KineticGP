%% Set path
userpath='C:/Users/Rudan/Documents/GitHub/KineticGP/';

addpath(genpath(strcat(userpath,'code/matlab/')))


filen=strcat("results/parameterization/param_top10/optim_simulation.mat");
load(filen);

%%
plotset=struct();
green_c=[0.4660 0.6740 0.1880];
blue_c=[0 0.4470 0.7410];
green_c=[0.4660 0.6740 0.1880]; % Green
orange_c=[0.8500, 0.3250, 0.0980];

plotset.sim_color=green_c; %purple
plotset.sim_shape='*';
plotset.sim_size=9;

plotset.line_width=1.5;
plotset.annot_size=8;
plotset.title_size=15;
plotset.title_color="black";
plotset.plot_measured=1;
plotset.plot_simulated=0;

label_size=15;
test_i=17;
plotset.title_text="";

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(1, 2, "TileSpacing", "tight");

nexttile(1)
plotset.meas_color=blue_c;
MF1_plot_curves_individual_acc(ACa22,"ACI",test_i,plotset)
hold on
plotset.meas_color=orange_c;
MF1_plot_curves_individual_acc(ACa23,"ACI",test_i,plotset)
% leg=legend({"Mean +/- std (2022)","Mean measurement (2022)","Mean +/- std (2023)","Mean measurement (2023)"},'Location','northeast');
ylabel("A (µmol/m^2/s)",'FontSize', label_size,'FontWeight','bold');
xlabel("C_a (µmol/mol)",'FontSize', label_size,'FontWeight','bold','Color',[0 0.4470 0.7410]);
ylim([0 50])

nexttile(2)
plotset.meas_color=blue_c;
MF1_plot_curves_individual_acc(GsCa22,"ACI",test_i,plotset)
hold on
plotset.meas_color=orange_c;
MF1_plot_curves_individual_acc(GsCa23,"ACI",test_i,plotset)
% leg=legend({"Mean +/- std (2022)","Mean measurement (2022)","Mean +/- std (2023)","Mean measurement (2023)"},'Location','northeast');
ylabel("g_s (mol/m^2/s)",'FontSize', label_size,'FontWeight','bold');
xlabel("C_a (µmol/mol)",'FontSize', label_size,'FontWeight','bold','Color',[0 0.4470 0.7410]);
ylim([0 0.3])

%%
filen="figures/tempFigures/MF1_A_and_gs_Ca.png";
set(gcf, 'PaperPosition', [0 0 8 3.5]); %
saveas(gcf,filen)


%%
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(1, 1, "TileSpacing", "tight");

nexttile(1)
plotset.meas_color=blue_c;
MF1_plot_curves_individual_acc(AQ22,"ACI",test_i,plotset)
hold on
plotset.meas_color=orange_c;
MF1_plot_curves_individual_acc(AQ23,"ACI",test_i,plotset)
% leg=legend({"Mean +/- std (2022)","Mean measurement (2022)","Mean +/- std (2023)","Mean measurement (2023)"},'Location','northeast');
ylabel("A (µmol/m^2/s)",'FontSize', label_size,'FontWeight','bold');
xlabel("PAR (µmol/m^2/s)",'FontSize', label_size,'FontWeight','bold','Color',[0 0.4470 0.7410]);
ylim([-5 45])
xlim([0 1850])
%%
filen="figures/tempFigures/MF1_A_PAR.png";
set(gcf, 'PaperPosition', [0 0 4 3.5]); %
saveas(gcf,filen)


%%
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t = tiledlayout(1, 1, "TileSpacing", "tight");

nexttile(1)
plotset.meas_color=blue_c;
plot_curves_individual_acc(AQ22,"ACI",test_i,plotset)
hold on
plotset.meas_color=orange_c;
plot_curves_individual_acc(AQ23,"ACI",test_i,plotset)
ylabel("A (µmol/m2/s)",'FontSize', label_size,'FontWeight','bold');
xlabel("PAR (µmol/m2/s)",'FontSize', label_size,'FontWeight','bold','Color',[0 0.4470 0.7410]);
ylim([-5 40])
xlim([0 1850])
leg=legend({"Mean +/- std (2022)","Mean measurement (2022)","Mean +/- std (2023)","Mean measurement (2023)"},'Location','northeast');

filen="../presentation0530/workflow_legend.png";
set(gcf, 'PaperPosition', [0 0 8 8]); %
saveas(gcf,filen)

