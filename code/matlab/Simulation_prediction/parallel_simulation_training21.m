function parallel_simulation_training21()

% userpath='C:/Users/Rudan/Documents/GitHub/';
userpath='/work/xu2/';
addpath(genpath(strcat(userpath,'KineticGP/code')))
addpath(genpath(strcat(userpath,'KineticGP/data')))

cd(strcat(userpath,'KineticGP/'))



method="BGLR";
scale="original";
topx=10:10:40;

delete(gcp('nocreate'))
parpool(4);
parfor x=1:4
    simulation_training21(topx(x),method,scale,"control")
    simulation_training21(topx(x),method,scale,"field")
end

