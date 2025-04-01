%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART2: Store the remaining testing set from 2022 and 2023 with only AQ curves

acc22=readtable("../data/curves/2022_maize_reference.csv");
acc23=readtable("../data/curves/2023_maize_reference.csv");

testing22=setdiff(acc22{:,"Accession"},final_acc);
testing23=setdiff(acc23{:,"Accession"},final_acc);

% length(intersect(testing22,testing23))
% testing_lines=readtable("../data/processed_data/Testing_237_genotypes.csv",'Delimiter',',','VariableNamingRule','preserve');
% testing_lines=testing_lines{:,1};

%%

measA22=zeros(length(testing22),6);
measA23=zeros(length(testing23),6);

measA22_all=[];
plots22=[];
acc22=[];
measA23_all=[];
plots23=[];
acc23=[];

for argind=1: length(testing22)
    display(argind)
    genotype=testing22{argind};
    
    [A,plots]=load_AQdata22(genotype,"testing");
    if ~isempty(A)
        measA22(argind,:)=mean(A,2,"omitnan");
        measA22_all=[measA22_all;A'];
        plots22=[plots22;plots'];
        acc22=[acc22;repmat(string(genotype),length(plots),1)];
    end
end


for argind=1: length(testing23)
    display(argind)
    genotype=testing23{argind};
    [A,plots]=load_AQdata23(genotype,"testing");
    if ~isempty(A)
        measA23(argind,:)=mean(A,2,"omitnan");
        measA23_all=[measA23_all;A'];
        plots23=[plots23;plots'];
        acc23=[acc23;repmat(string(genotype),length(plots),1)];
    end
end

%%
measA_all=[measA22_all;measA23_all];
plotrep=[plots22;plots23];
accs=[acc22;acc23];
years=[2022*ones(length(plots22),1);2023*ones(length(plots23),1)];

sim_var=["Year";"Accession";"Plot_rep";strcat("PAR_",string([1800,1100,500,300,150,50]))'];

writetable(array2table([years,accs,plotrep,measA_all],"VariableNames",sim_var),"../data/processed_data/Testing_AQcurves_years22&23_plot.csv",'WriteVariableNames',true,'WriteRowNames',true);

%%

measA=[measA22;measA23];
accs=[testing22;testing23];
years=[2022*ones(length(testing22),1);2023*ones(length(testing23),1)];

sim_var=["Year";"Accession";strcat("PAR_",string([1800,1100,500,300,150,50]))'];

writetable(array2table([years,accs,measA],"VariableNames",sim_var),"../data/processed_data/Testing_AQcurves_years22&23_accession.csv",'WriteVariableNames',true,'WriteRowNames',true);
