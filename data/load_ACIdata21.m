function [Atot,plot_acc]=load_ACIdata21(genotype)
%%
acc_info=readtable("../data/curves/2021_maize_reference.csv",'VariableNamingRule','preserve');
acc_info{:,"Accession"}=replace(acc_info{:,"Accession"},'_','_00');

acc=string(acc_info{:,"Accession"});
uniq_acc=unique(acc);

curve_data=readtable("../data/curves/2021_ACi_rawData_maize.csv",'VariableNamingRule','preserve');

% check for which accessions we have plots
check_plot=zeros(length(uniq_acc),2);
for i=1:length(uniq_acc)
    plots_i=acc_info{strcmp(uniq_acc(i),acc_info{:,"Accession"}),"Plot"};
    for j=1:2
        % check for each unique accession if their ACI curve is available
        % in curve_data
        if sum(curve_data{:,"Plot"}==plots_i(j))>0
            check_plot(i,j)=plots_i(j);
        end
    end
end

uniq_acc=uniq_acc(sum(check_plot,2)~=0);
acc_plot=check_plot(sum(check_plot,2)~=0,:);


% Change the position of points:
% Some plots are measured with decreasing CO2 first, restabilization and then with increasing CO2,
% which is reverse to the majority of cases

changeplots=[2246,2246,1010,1188,1085,1228];
changereps=[1,2,2,1,2,2];

for k=1:length(changeplots)
    inter=intersect(find(curve_data{:,"Plot"}==changeplots(k)),find(curve_data{:,"Repeat"}==changereps(k)));
    temp=curve_data(inter,:);
    curve_data(inter(1:6),:)=temp(8:13,:);
    curve_data(inter(7:13),:)=temp(1:7,:);
end

inter=intersect(find(curve_data{:,"Plot"}==1064),find(curve_data{:,"Repeat"}==2));
temp=curve_data(inter,:);
curve_data(inter(1:5),:)=temp(8:12,:);
curve_data(inter(6:12),:)=temp(1:7,:);


% Remove the points with unreasonable values
curve_data=curve_data(curve_data{:,"Flag_removal"}~="x",:); % no time data is available

%
% REMOVE ACCESSION 9 "SSA_009", only 2 accessions with large variance
% 18: "SSA_043" single sample available, zero standard deviation
% 45: "SSA_171": bad measurement
% 52:  "SSA_208" bad measurement
% 64: "SSA_272" single sample 
% 71: "SSA_304" single sample
% 78: "SSA_369" single sample

[~,ind_remove]=ismember(["SSA_00009","SSA_00043","SSA_00171","SSA_00208","SSA_00272","SSA_00304","SSA_00369"],uniq_acc);
uniq_acc(ind_remove)=[];
acc_plot(ind_remove,:)=[];

% [~,final_acc,param_name,vmaxind]=optim_initialization_parameters();
% final_acc(9)=[];

%
% cases
% SSA_001: repeated measurements at 400 (usually the flag point)
% SSA_002 (p2,r3), SSA_003(p1,r3): repeated measurements in more CO2 levels
% SSA_002 (p2,r1), SSA_008: missing measurements  at some points 

% SSA_024 (rep3): higher COlevels measured (1400-1600)

% remove specific replicates:
% SSA_139, weird curve

% 1169, r=3, missing at 6nd point
% 2198, r=2 missing the first point
% SSA_156, weird curve

% SSA_235, duplicated replicates

% SSA_235, weird curve
% ind=intersect(find(curve_data{:,"Plot"}==2198),find(curve_data{:,"Repeat"}==2));
% curve_data=curve_data(setdiff(1:size(curve_data,1),ind),:); % no time data is available

removeplots=[2226,1085,1228,2174,2246,2246,1250,1133,2057,2198,1169,2227,1010,1149,1082];
removereps=[1,3,2,2,2,3,3,3,2,2,3,1,2,1,1];

for k=1:length(removeplots)
    ind=intersect(find(curve_data{:,"Plot"}==removeplots(k)),find(curve_data{:,"Repeat"}==removereps(k)));
    curve_data=curve_data(setdiff(1:size(curve_data,1),ind),:); % no time data is available
end

%%

[~,arg_ind]=ismember(genotype,uniq_acc);

plot_acc=[];

if arg_ind~=0
    % figure
    Atot=[];
    Ctot=[];
    gstot=[];
    for p=1:2
        plot_i=acc_plot(arg_ind,p);
        if plot_i~=0
            for r=1:3
                ind=intersect(find(curve_data{:,"Plot"}==plot_i),find(curve_data{:,"Repeat"}==r));
                A=curve_data{ind,"A"};
                C=curve_data{ind,"Ca"};
                gs=curve_data{ind,"gsw"}/1.6;
                if length(A)==12
                    Atot=[Atot,A];
                    Ctot=[Ctot,C];
                    gstot=[gstot,gs];
                    plot_acc=[plot_acc;strcat(string(plot_i),"_",string(r))];
                    % plot(C,A,'o')
                    % hold on
                end
            end
        end
    end
    
    
    Atot(6,:)=[];
    gstot(6,:)=[];
    
    % final_A=mean(Atot,2,"omitnan");
    % final_Asd=std(Atot,0,2,"omitnan");
    % % final_C=mean(C_t,2,"omitnan");
    % % final_Ci=mean(Ci_t,2,"omitnan");
    % final_gs=mean(gstot,2,"omitnan");
    % final_gssd=std(gstot,0,2,"omitnan");
    % 
    
else
    Atot=NaN(11,1);
end

