%% Add path of PESTO optimizer and C4 model
% userpath='C:/Users/Rudan/Documents/GitHub/';

userpath='/work/xu2/';
addpath(genpath(strcat(userpath,'PESTO-master/')))

addpath(genpath(strcat(userpath,'KineticGP/code')))
addpath(genpath(strcat(userpath,'KineticGP/data')))


cd(strcat(userpath,'KineticGP/'))
delete(gcp('nocreate'))

parpool(4);
parfor x=1:4
    SFX2_prepare_figures(x)
end
