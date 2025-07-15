function SFsec32to37_plot_curves(year,y_axis)
%% Load simulated curves


top=10:10:40;

tops1=[11, 20, 34, 44];
degs=56-tops1;

data=load("data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;
orange_c=[0.8500, 0.3250, 0.0980]; % orange
blue_c=[0 0.4470 0.7410];
green_c=[0.4660 0.6740 0.1880]; % Green
purple_c=[0.5647 0.4039 0.6549]; %purple
colors = [
    245/255, 212/255, 145/255;  % #F5D491
     41/255, 114/255,  59/255;  % #29723B
    166/255,  46/255,  56/255;  % #A62E38
     59/255,  91/255, 165/255   % #3B5BA5
];

%% Plot

boxs=[.085, .02];
posx_inter=0:6;

if y_axis=='GsCa'
    xstart=0.2;
else
    xstart=0.12;
end
posx1=posx_inter*0.112+xstart;

posxy=zeros(70,4,4);
for t=1:4
    ystart=[0.91 0.90 0.89 0.88];

    posy_inter=0:9;
    posy1=ystart(t)-posy_inter*0.0835;

    for i=1:10
        if i==1
            posxy(((i-1)*7+1):(i*7),:,t)=[posx1', repmat(ystart(t),7,1), repmat(boxs,7,1)];
        else
            posxy(((i-1)*7+1):(i*7),:,t)=[posx1', repmat(posy1(i),7,1), repmat(boxs,7,1)];
        end
    end
end

sim_shape='*';
shape_size=7;
line_width=1.5;
annot_size=8;
title_size=9;


%%
j=0;
x_count=0:9;
x_count=x_count*7+1;
%%
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
myplot = tiledlayout(10, 7, "TileSpacing", "tight");
for arg_ind=1: 68

    nexttile(arg_ind)
    j=j+1;
    
    acc_lab=replace(final_acc(arg_ind),"_","-");
    title_text=strcat(string(arg_ind),".",acc_lab);
    for t=1:4
        topX=top(t);
        sim_data=load(strcat("results/parameterization/param_top",string(topX),"/optim_simulation.mat"));

        switch y_axis
            case 'ACI'
                data1=sim_data.ACa22;
                data2=sim_data.ACa23;
            case 'GsCa'
                data1=sim_data.GsCa22;
                data2=sim_data.GsCa23;
            case 'AQ'
                data1=sim_data.AQ22;
                data2=sim_data.AQ23;
        end

        if year==2022
            data=data1;
        else
            data=data2;
        end
        
        sim=data.sim(:,arg_ind);
        meas=data.meas(:,arg_ind);
        sd=data.sd(:,arg_ind);
        z=data.z(arg_ind)/degs(t);
        %%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%
        
        yu=meas+sd;
        yl=meas-sd;
        
        x=data.x(:,arg_ind);% x must be a column
        
        yu=yu';
        yl=yl'; % yl AND yu must be rows.
        
        [~,idx] = sort(x);
        x=x(idx);
        yu=yu(idx);
        yl=yl(idx);
        measured=meas(idx);
        simulated=sim(idx);
        % simulated0=simA_ini(idx);

        xfill=[x',flipud(x)']; % xfill must be a single row: the second curve is reversed. It should be like [1 2 3 3 2 1] It fills between two curves. It must be in the same row, and
        yfill=[yl, fliplr(yu)];
        

        meas_c=purple_c;

        if t==1
            plot(x,measured,'o','color',meas_c,'MarkerSize',shape_size-2, 'linewidth',line_width)
            hold on;
            fill(xfill,yfill,meas_c,'LineStyle','none','FaceAlpha',.2);
        end
        % hold on 
        % plot(x,simulated0,'*','color',green_c,'MarkerSize',shape_size, 'linewidth',line_width)
        
        hold on 
        plot(x,simulated,sim_shape,'color',colors(t,:),'MarkerSize',shape_size, 'linewidth',line_width)
        
        str = string(round(z,2));
        a=annotation('textbox',posxy(j,:,t),'String',str, 'EdgeColor','none','FitBoxToText','on','Color',colors(t,:)); %
        a.FontSize = annot_size;
        if strcmp(y_axis,"ACI")
            ylim([0,65])
            xlim([0,1300])
        elseif strcmp(y_axis,"GsCa")
            ylim([0, 0.4])
            xlim([0,1300])
        elseif strcmp(y_axis,"AQ")
            ylim([-5,50])
            xlim([0,1850])
        end
        hold on
    end
    
    title(title_text,'FontSize', title_size,'FontWeight','bold') 

    if j<64
        xticks([])
    end
    if sum(x_count==j)<1
        yticks([])
    end
end


switch y_axis
    case 'ACI'
        ylab="A (µmol/m2/s)";
        xlab="C_a (µmol/mol)";

    case 'GsCa'
        ylab="g_s (mol/m2/s)";
        xlab="C_a (µmol/mol)";
    case 'AQ'
        ylab="A (µmol/m2/s)";
        xlab="PAR (µmol/mol)";
end

%%
xlabel(myplot,xlab,'FontSize', 10,'FontWeight','bold');
ylabel(myplot,ylab,'FontSize', 10,'FontWeight','bold');
leg=legend({"Measurement","Measurement +/- std","Simulation (N=11)","Simulation (N=20)","Simulation (N=34)","Simulation (N=44)"},'Location','northeastoutside');
