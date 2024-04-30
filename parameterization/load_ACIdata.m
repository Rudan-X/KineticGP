function [final_C,final_Ci,final_A,final_Asd,final_gs,final_gssd]=load_ACIdata(genotype)
%%
acc_info=readtable("../data/2021_maize_reference.csv",'VariableNamingRule','preserve');
acc=string(acc_info{:,"Accession"});
uniq_acc=unique(acc);

curve_data=readtable("../data/2021_ACi_rawData_maize.csv",'VariableNamingRule','preserve');

% curve_data=curve_data(curve_data{:,"Flag_removal"}~="x",:); % no time data is available
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

% REMOVE ACCESSION 9 "SSA_009", only 2 accessions with large variance?
% keep it
% 18: "SSA_043" single sample available, zero standard deviation
% 45: "SSA_171": bad measurement
% 52:  "SSA_208" bad measurement
% 64: "SSA_272" single sample 
% 71: "SSA_304" single sample
% 78: "SSA_369" single sample

[~,ind_remove]=ismember(["SSA_043","SSA_171","SSA_208","SSA_272","SSA_304","SSA_369"],uniq_acc);
uniq_acc(ind_remove)=[];
acc_plot(ind_remove,:)=[];

% same_size=zeros(length(uniq_acc),1);
toremove_new=["1180_3","1250_3","1149_1","2287_3","2313_1","2246_2","2246_3","2226_1"];

% remove outlier points
inter=intersect(find(curve_data{:,"Plot"}==1010),find(curve_data{:,"Repeat"}==1));
inter=intersect(inter,find(curve_data{:,"Time"}==duration(10,03,22)));
curve_data(inter,:)=[];

inter=intersect(find(curve_data{:,"Plot"}==1059),find(curve_data{:,"Repeat"}==1));
inter=intersect(inter,find(curve_data{:,"Time"}==duration(13,09,33)));
curve_data(inter,:)=[];

inter=intersect(find(curve_data{:,"Plot"}==2238),find(curve_data{:,"Repeat"}==2));
inter=intersect(inter,find(curve_data{:,"Time"}==duration(13,44,49)));
curve_data(inter,:)=[];

inter=intersect(find(curve_data{:,"Plot"}==2254),find(curve_data{:,"Repeat"}==1));
inter1=intersect(inter,find(curve_data{:,"Time"}==duration(13,07,31)));
inter2=intersect(inter,find(curve_data{:,"Time"}==duration(13,13,04)));
curve_data([inter1,inter2],:)=[];

inter=intersect(find(curve_data{:,"Plot"}==1228),find(curve_data{:,"Repeat"}==2));
inter1=intersect(inter,find(curve_data{:,"Time"}==duration(14,43,14)));
inter2=intersect(inter,find(curve_data{:,"Time"}==duration(14,44,42)));
curve_data([inter1,inter2],:)=[];

inter=intersect(find(curve_data{:,"Plot"}==2054),find(curve_data{:,"Repeat"}==1));
inter=intersect(inter,find(curve_data{:,"Time"}==duration(14,55,19)));
curve_data(inter,:)=[];


inter=intersect(find(curve_data{:,"Plot"}==1082),find(curve_data{:,"Repeat"}==1));
inter=intersect(inter,find(curve_data{:,"Time"}==duration(13,22,41)));
curve_data(inter,:)=[];


inter=intersect(find(curve_data{:,"Plot"}==2226),find(curve_data{:,"Repeat"}==3));
inter=inter(13:19);
curve_data(inter,:)=[];
 

inter=intersect(find(curve_data{:,"Plot"}==1085),find(curve_data{:,"Repeat"}==3));
inter=inter(14:19);
curve_data(inter,:)=[];

inter=intersect(find(curve_data{:,"Plot"}==2174),find(curve_data{:,"Repeat"}==2));
inter=inter(13:26);
curve_data(inter,:)=[];

% % change the position of points:
inter=intersect(find(curve_data{:,"Plot"}==2246),find(curve_data{:,"Repeat"}==1));
temp=curve_data(inter,:);
curve_data(inter(1:6),:)=temp(8:13,:);
curve_data(inter(7:13),:)=temp(1:7,:);

inter=intersect(find(curve_data{:,"Plot"}==1064),find(curve_data{:,"Repeat"}==2));
temp=curve_data(inter,:);
curve_data(inter(1:5),:)=temp(8:12,:);
curve_data(inter(6:12),:)=temp(1:7,:);

inter=intersect(find(curve_data{:,"Plot"}==1010),find(curve_data{:,"Repeat"}==2));
temp=curve_data(inter,:);
curve_data(inter(1:6),:)=temp(8:13,:);
curve_data(inter(7:13),:)=temp(1:7,:);

inter=intersect(find(curve_data{:,"Plot"}==1188),find(curve_data{:,"Repeat"}==1));
temp=curve_data(inter,:);
curve_data(inter(1:6),:)=temp(8:13,:);
curve_data(inter(7:13),:)=temp(1:7,:);

inter=intersect(find(curve_data{:,"Plot"}==1085),find(curve_data{:,"Repeat"}==2));
temp=curve_data(inter,:);
curve_data(inter(1:6),:)=temp(8:13,:);
curve_data(inter(7:13),:)=temp(1:7,:);


inter=intersect(find(curve_data{:,"Plot"}==1228),find(curve_data{:,"Repeat"}==2));
temp=curve_data(inter,:);
curve_data(inter(1:6),:)=temp([8:13],:);
curve_data(inter(7:13),:)=temp(1:7,:);

% add NA points
toadd=["1064_2","2174_2","1169_3"];
positiontoadd={6,1,6};
% remove extra points
todelete=["1037_2"];
positiontodelete={8};


t_inter_ca=[10:1.5: 17.5, 20:1.5: 29]; % for Ca input, 13 timepoints, including the drop to ambient CO2
t_inter_A=[10:1.5: 16, 20:1.5: 29]; % for A, only 12 timepoints are compared

tt=NaN(1,1);
A_t=NaN(1,1);
C_t=NaN(1,1);
Ci_t=NaN(1,1);
gs_t=NaN(1,1);
col_i=0;
%%
[~,arg_ind]=ismember(genotype,uniq_acc);

for p=1:2
    plot_i=acc_plot(arg_ind,p);
    if plot_i~=0
        for r=1:3
            ind=intersect(find(curve_data{:,"Plot"}==plot_i),find(curve_data{:,"Repeat"}==r));
            plot_repeat=strcat(string(plot_i),"_",string(r));
            if ~isempty(ind) && ~contains(plot_repeat,toremove_new)
                col_i=col_i+1;
                time=curve_data{ind,"Time"};
                ts=seconds(time-time(1));
                tt(length(ts),col_i)=NaN;
                tt(1:length(ts),col_i)=ts;
                A=curve_data{ind,"A"};
                C=curve_data{ind,"Ca"};
                Ci=curve_data{ind,"Ci"};
                gs=curve_data{ind,"gsw"}/1.6;
                if contains(plot_repeat,toadd)
                    [~,ind_add]=ismember(plot_repeat,toadd);
                    ind_add=positiontoadd{ind_add};
                    temp=A;
                    A=zeros(length(A)+length(ind_add),1);
                    A(ind_add)=NaN;
                    A(setdiff(1:length(A),ind_add))=temp;
                    temp=C;
                    C=zeros(length(C)+length(ind_add),1);
                    C(ind_add)=NaN;
                    C(setdiff(1:length(C),ind_add))=temp;
                    temp=Ci;
                    Ci=zeros(length(Ci)+length(ind_add),1);
                    Ci(ind_add)=NaN;
                    Ci(setdiff(1:length(Ci),ind_add))=temp;
                    temp=gs;
                    gs=zeros(length(gs)+length(ind_add),1);
                    gs(ind_add)=NaN;
                    gs(setdiff(1:length(gs),ind_add))=temp;
                end

                if contains(plot_repeat,todelete)
                    [~,ind_delete]=ismember(plot_repeat,todelete);
                    A(positiontodelete{ind_delete})=[];
                    C(positiontodelete{ind_delete})=[];
                    Ci(positiontodelete{ind_delete})=[];
                    gs(positiontodelete{ind_delete})=[];
                end
                
                A_t(length(ts),col_i)=NaN;
                A_t(1:length(A),col_i)=A;
 
                C_t(length(ts),col_i)=NaN;
                C_t(1:length(C),col_i)=C;

                Ci_t(length(ts),col_i)=NaN;
                Ci_t(1:length(Ci),col_i)=Ci;

                gs_t(length(ts),col_i)=NaN;
                gs_t(1:length(gs),col_i)=gs;
            end
        end
    end
end

A_t(isnan(A_t))=0;
C_t(isnan(C_t))=0;
Ci_t(isnan(Ci_t))=0;
gs_t(isnan(gs_t))=0;

A_t(sum(A_t,2)==0,:)=[];
C_t(sum(C_t,2)==0,:)=[];
Ci_t(sum(Ci_t,2)==0,:)=[];
gs_t(sum(gs_t,2)==0,:)=[];

A_t(A_t==0)=NaN;
gs_t(A_t==0)=NaN;

C_t(C_t==0)=NaN;
Ci_t(Ci_t==0)=NaN;

A_t(6,:)=[];
gs_t(6,:)=[];


final_A=mean(A_t,2,"omitnan");
final_Asd=std(A_t,0,2,"omitnan");
final_C=mean(C_t,2,"omitnan");
final_Ci=mean(Ci_t,2,"omitnan");
final_gs=mean(gs_t,2,"omitnan");
final_gssd=std(gs_t,0,2,"omitnan");

% plot(final_C([1:5,7:size(final_C,1)]),final_A,'o')
