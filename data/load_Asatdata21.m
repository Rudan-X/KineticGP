function [Asat,plot_acc]=load_Asatdata21(genotype)

%%

acc_info=readtable("../data/curves/2021_maize_reference.csv",'VariableNamingRule','preserve');
acc=string(acc_info{:,"Accession"});

acc_long=cell(length(acc),1);
for i=1:length(acc)
    temp=split(acc(i),"_");
    acc_long{i}=strcat(temp(1),"_00",temp(2));
end
acc_long=string(acc_long);

plots_acc=acc_info{contains(acc_long,genotype),"Plot"};

curve_data=readtable("../data/curves/2021_lightSaturatedGasExchange_maize.csv",'VariableNamingRule','preserve');
plot_all=curve_data{:,"Plot"};
table_ind=[];
for p=1:length(plots_acc)
    table_ind=[table_ind;find(plot_all==plots_acc(p))];
end

Asat=curve_data{table_ind,"A_sat"};
plot_acc=strcat(string(curve_data{table_ind,"Plot"}),"_",string(curve_data{table_ind,"Repeat"}));
% Asat=mean(Asat,"omitnan");
