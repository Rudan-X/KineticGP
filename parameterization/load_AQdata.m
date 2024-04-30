function [final_A,final_Asd,final_Q, final_C]=load_AQdata(genotype)

%%
problem=[];
% for a=1:72
    % genotype=final_acc(a);
    uniq_acc=split(genotype,'_'); 
    genotype=strcat(uniq_acc(1),"_00",uniq_acc(2));
    
    acc_info=readtable("../data/2022_maize_reference.csv",'VariableNamingRule','preserve');
    acc=string(acc_info{:,"Accession"});
    
    plots_acc=acc_info{contains(acc,genotype),"Plot"};
    
    
    curve_data=readtable("../data/2022_AQcurves_maize.csv",'VariableNamingRule','preserve');
    curve_data=curve_data(curve_data{:,"Flag_removal"}~="x",:); % no time data is available
    plot_all=curve_data{:,"Plot"};
    table_ind=[];
    for p=1:length(plots_acc)
        table_ind=[table_ind;find(plot_all==plots_acc(p))];
    end
    Plot_Rep=unique(curve_data{table_ind,"Plot_Repeat"});
    % temp=split(Plot_Rep,'_');
    % ploti=temp(:,1);
    % repi=temp(:,2);
    
    A=zeros(7,length(Plot_Rep));
    C=A;
    Q=A;
    for i=1:length(Plot_Rep)
        ind_i=contains(curve_data{:,"Plot_Repeat"},Plot_Rep(i));
        if sum(ind_i)==7
            A(:,i)=curve_data{ind_i,"Photo"};
            Q(:,i)=curve_data{ind_i,"PAR"};
            C(:,i)=curve_data{ind_i,"Ci"};
        else
            problem=[problem;Plot_Rep(i)];
            
        end
    end
    keep=sum(A,1)~=0;
    A=A(1:(end-1),keep);
    C=C(1:(end-1),keep);
    Q=Q(1:(end-1),keep);

    % subplot(6,12,a)
    % 
    % for i=1:size(A,2)
    %     plot(Q(:,i),A(:,i),'o')
    % 
    %     hold on
    % end
    % title(final_acc(a))
% end
    final_A=mean(A,2,"omitnan");
    final_Asd=std(A,0,2,"omitnan");
    final_Q=mean(Q,2,"omitnan");
    final_C=mean(C,2,"omitnan");
