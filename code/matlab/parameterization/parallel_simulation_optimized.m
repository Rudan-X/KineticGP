function parallel_simulation_optimized()

% userpath='C:/Users/Rudan/Documents/GitHub/';
userpath='/work/xu2/';
addpath(genpath(strcat(userpath,'KineticGP/code')))
addpath(genpath(strcat(userpath,'KineticGP/data')))

cd(strcat(userpath,'KineticGP/'))

topx=10:10:40;


delete(gcp('nocreate'))
parpool(4);
parfor x=1:4
    simulate_optim_result(topx(x))
end