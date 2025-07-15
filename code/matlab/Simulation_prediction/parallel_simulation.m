function parallel_simulation()

% userpath='C:/Users/Rudan/Documents/GitHub/';
userpath='/work/xu2/';
addpath(genpath(strcat(userpath,'KineticGP/code')))
addpath(genpath(strcat(userpath,'KineticGP/data')))

cd(strcat(userpath,'KineticGP/'))

method="rrBLUP";
scale="original";


topx=10:10:40;

delete(gcp('nocreate'))
parpool(4);
parfor x=1:4
    simulation_predicted_parameters_control_mean_and_accuracy(topx(x),method,scale)
end


topx=10:10:40;
delete(gcp('nocreate'))
parpool(4);
parfor x=1:4
    simulation_predicted_parameters_burnin_accuracy(topx(x),method,scale)
end


for x=1:2
    simulation_BLUPpredicted_parameters(topx(x),"rrBLUP",scale)
end


for x=1:2
    ACI_simulation_predicted_parameters(topx(x),"rrBLUP",scale)
end
