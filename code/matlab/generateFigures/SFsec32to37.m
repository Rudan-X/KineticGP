%% Simulate optimized profiles
userpath='C:/Users/Rudan/Documents/GitHub/';
addpath(genpath(strcat(userpath,'KineticGP/code')))
addpath(genpath(strcat(userpath,'KineticGP/data')))

cd(strcat(userpath,'KineticGP/'))

%% Plot simulated and measured profiles

years=[2022,2023];
y=2;
year=years(y);
files=["SF_01ACi","SF_01gs","SF_01AQ"];
curves=["ACI","GsCa","AQ"];
c=3;
curve=curves(c);

SFsec32to37_plot_curves(year,curves(c))


%%
filen=strcat("figures/SFsec3.",string(1+y+(c-1)*2),"_", curve, string(year),".png");

set(gcf, 'PaperPosition', [0 0 10 15]); %
saveas(gcf,filen)

%%
data=load("data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;
load(strcat("results/parameterization/param_top10/optim_simulation.mat"));


mat=[ ACa22.sim',ACa23.sim',GsCa22.sim',GsCa23.sim',AQ22.sim',AQ23.sim'];

sim_var=["Accession"; 
    strcat("ACa2022_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    strcat("ACa2023_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';

    strcat("GsCa2022_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    strcat("GsCa2023_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';

    strcat("2022_PAR_",string([1800,1100,500,300,150,50]))';
    strcat("2023_PAR_",string([1800,1100,500,300,150,50]))'];

T = array2table([final_acc,mat],"VariableNames",sim_var);

writetable(T,strcat("results/parameterization/param_top10/Simulated_profiles.csv"),'WriteVariableNames',true,'WriteRowNames',true);
%%

mat=[ ACa22.meas',ACa23.meas',GsCa22.meas',GsCa23.meas',AQ22.meas',AQ23.meas'];

sim_var=["Accession"; 
    strcat("ACa2022_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    strcat("ACa2023_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    
    strcat("GsCa2022_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    strcat("GsCa2023_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';

    strcat("AQ2022_",string([1800,1100,500,300,150,50]))';
    strcat("AQ2023_",string([1800,1100,500,300,150,50]))'];

T = array2table([final_acc,mat],"VariableNames",sim_var);

writetable(T,strcat("results/parameterization/param_top10/Measured_profiles.csv"),'WriteVariableNames',true,'WriteRowNames',true);