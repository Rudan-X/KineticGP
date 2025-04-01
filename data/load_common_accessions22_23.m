function [common_acc,uniq_acc22,acc_plot22,uniq_acc23,acc_plot23]=load_common_accessions22_23()
%% ACI and AQ curves from 2022
acc_info=readtable("../data/curves/2022_maize_reference.csv",'VariableNamingRule','preserve');
acc=string(acc_info{:,"Accession"});
uniq_acc=unique(acc);

ACI22=readtable("../data/curves/2022_ACi_rawData_maize.csv",'VariableNamingRule','preserve');
AQ22=readtable("../data/curves/2022_AQcurves_maize.csv",'VariableNamingRule','preserve');
AQ22=AQ22(AQ22{:,"Flag_removal"}~="x",:); % should be removed 

% check for which accessions we have plots
check_ACI=zeros(length(uniq_acc),2);
check_AQ=zeros(length(uniq_acc),2);
for i=1:length(uniq_acc)
    plots_i=acc_info{strcmp(uniq_acc(i),acc_info{:,"Accession"}),"Plot"};
    for j=1:2
        % check for each unique accession if their ACI/AQ curve is available
        % in curve_data
        if sum(ACI22{:,"Plot"}==plots_i(j))>0
            check_ACI(i,j)=plots_i(j);
        end

        if sum(AQ22{:,"Plot"}==plots_i(j))>0
            check_AQ(i,j)=plots_i(j);
        end
    end
end


uniq_acc22=uniq_acc(sum(check_ACI,2)~=0 & sum(check_AQ,2)~=0);
acc_plot22=check_ACI(sum(check_ACI,2)~=0 & sum(check_AQ,2)~=0,:);


%% ACI and AQ curves from 2023
acc_info=readtable("../data/curves/2023_maize_reference.csv",'VariableNamingRule','preserve');
acc=string(acc_info{:,"Accession"});
uniq_acc=unique(acc);

ACI23=readtable("../data/curves/2023_ACi_rawData_maize.csv",'VariableNamingRule','preserve');
AQ23=readtable("../data/curves/2023_AQcurves&lightSaturatedGasExchange_maize.csv",'VariableNamingRule','preserve');

% check for which accessions we have plots
check_ACI=zeros(length(uniq_acc),2);
check_AQ=zeros(length(uniq_acc),2);
for i=1:length(uniq_acc)
    plots_i=acc_info{strcmp(uniq_acc(i),acc_info{:,"Accession"}),"Plot"};
    for j=1:2
        % check for each unique accession if their ACI/AQ curve is available
        % in curve_data
        if sum(ACI23{:,"Plot"}==plots_i(j))>0
            check_ACI(i,j)=plots_i(j);
        end
        if sum(contains(AQ23{:,'plot_id'},string(plots_i(j))))>0
            check_AQ(i,j)=plots_i(j);
        end
    end
end


uniq_acc23=uniq_acc(sum(check_ACI,2)~=0 & sum(check_AQ,2)~=0);
acc_plot23=check_ACI(sum(check_ACI,2)~=0 & sum(check_AQ,2)~=0,:);

%%
common_acc=intersect(uniq_acc22,uniq_acc23);

remove_acc=["SSA_00002","SSA_00016","SSA_00055","SSA_00246","SSA_00263","SSA_00281","SSA_00285","SSA_00301","SSA_00308","SSA_00343",...
    "SSA_00031","SSA_00046","SSA_00078"];

[~,ind]=ismember(remove_acc,common_acc);
common_acc(ind)=[];
%%
% Curves 2022:
% SSA_00031 no measurements for Cref=1250 for 2 cases
% SSA_00046 with 2 measurements, moderate std
% SSA_00055 with 2 measurements, moderate std
% SSA_00078 no measurements for Cref=1250 for 2 cases
%%
% Curves 2023:
% SSA_00002 with 2 measurements, large std
% SSA_00016 with single measurements
% SSA_00055 with 2 measurements, one weird point
% SSA_00246 weird curves
% SSA_00263 single measurement
% SSA_00281 with 2 measurements, large std
% SSA_00285  with 4 measurements, large std
% SSA_00301 with 2 measurements, large std, not stable
% SSA_00308  with 3 measurements, large std
% SSA_00343 with 4 measurements, large std
