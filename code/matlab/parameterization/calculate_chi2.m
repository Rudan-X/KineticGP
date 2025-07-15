userpath='C:/Users/Rudan/Documents/GitHub/';
addpath(genpath(strcat(userpath,'KineticGP/code')))
addpath(genpath(strcat(userpath,'KineticGP/data')))

top=10:10:40;


for x=1:4
    topX=top(x);
    estim_ind=optimized_var_ind(topX);
    load(strcat("results/parameterization/param_top",string(topX),"/optim_simulation.mat"))

    error=[real(ACa22.z)',real(GsCa22.z)',real(AQ22.z)',real(ACa23.z)',real(GsCa23.z)',real(AQ23.z')];
    T = array2table(error,'VariableNames',["ACa22","GsCa22","AQ22","ACa23","GsCa23","AQ23"]);
    writetable(T,strcat("results/parameterization/param_top",string(topX),"/original_chi_square.csv"),'WriteVariableNames',true);

    % significant_e=chi2inv(0.05,56-length(estim_ind));
    % error=[real(ACa22.z)',real(GsCa22.z)',real(AQ22.z)',real(ACa23.z)',real(GsCa23.z)',real(AQ23.z')]/significant_e;
    % T = array2table(error,'VariableNames',["ACa22","GsCa22","AQ22","ACa23","GsCa23","AQ23"]);
    % writetable(T,strcat("results/parameterization/param_top",string(topX),"/scaled_chi_square.csv"),'WriteVariableNames',true);
    % 

    % deg=56-length(estim_ind);
    % error=[real(ACa22.z)',real(GsCa22.z)',real(AQ22.z)',real(ACa23.z)',real(GsCa23.z)',real(AQ23.z')]/deg;
    % T = array2table(error,'VariableNames',["ACa22","GsCa22","AQ22","ACa23","GsCa23","AQ23"]);
    % writetable(T,strcat("results/parameterization/param_top",string(topX),"/reduced_chi_square.csv"),'WriteVariableNames',true);
end


thre=zeros(4,1);
degs=zeros(4,1);
nvars=zeros(4,1);
for x=1:4
    topX=top(x);
    estim_ind=optimized_var_ind(topX);

    % thre(x)=chi2inv(0.05,56-length(estim_ind));
    degs(x)=56-length(estim_ind);
    nvars(x)=length(estim_ind);
end

