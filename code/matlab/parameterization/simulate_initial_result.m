function simulate_initial_result()

KE_type="original";

init_sol=load_initial_solution();

[simA22,simgs22]=simulate_ACI(init_sol,27.7,KE_type);

[simAQ22]=simulate_AQ(init_sol,27.7,KE_type);

mat=[simA22',simA22',simgs22',simgs22',simAQ22',simAQ22'];


sim_var=["Accession";
    strcat("ACa2002_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    strcat("ACa2003_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    strcat("GsCa2022_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';
    strcat("GsCa2023_",string([400, 600, 800, 1000, 1250,300,250,200, 100,75, 25]))';

    strcat("AQ2022_",string([1800,1100,500,300,150,50]))';
    strcat("AQ2023_",string([1800,1100,500,300,150,50]))'];

data=load("data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;

mat=repmat(mat,length(final_acc),1);
T = array2table([final_acc,mat],"VariableNames",sim_var);

writetable(T,"results/parameterization/default/default_single_simulation2seasons.csv",'WriteVariableNames',true,'WriteRowNames',true);

 


end

