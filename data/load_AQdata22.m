function [A,Plot_Rep]=load_AQdata22(genotype,data_type)
AQ22=readtable("../data/curves/2022_AQcurves_maize.csv",'VariableNamingRule','preserve');
AQ22=AQ22(AQ22{:,"Flag_removal"}~="x",:); % should be removed 

%%
if data_type=="training"
    [~,uniq_acc22,acc_plot22,~,~]=load_common_accessions22_23();
    plots=acc_plot22(contains(uniq_acc22,genotype),:);
elseif data_type=="testing"
    acc_info=readtable("../data/curves/2022_maize_reference.csv",'VariableNamingRule','preserve'); 
    plots=acc_info{strcmp(acc_info{:,"Accession"},genotype),"Plot"};
    plots=plots';
end
if ~isempty(plots)
    Plot_Rep=strcat(string(repelem(plots,3)),'_',string(repmat(1:3,1,2)));
    % Plot_Rep with 0_1, 0_2 or 0_3 means that plot is not found in the
    % dataset
    
    A=zeros(7,length(Plot_Rep));
    for i=1:length(Plot_Rep)
        ind_i=contains(AQ22{:,"Plot_Repeat"},Plot_Rep(i));
        if sum(ind_i)==7
            A(:,i)=AQ22{ind_i,"Photo"};
        end
    end
    keep=sum(A,1)~=0;
    A=A(1:(end-1),keep);
    
    Plot_Rep=Plot_Rep(keep);
    %%
    % final_A=mean(A,2,"omitnan");
    % final_Asd=std(A,0,2,"omitnan");

else
    final_A=[];
    final_Asd=[];
end
%%
% 31,40