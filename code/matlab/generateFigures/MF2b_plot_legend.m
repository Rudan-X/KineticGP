function  [] = MF3b_plot_legend(year,curve,y_axis,arg_ind,sim_color,sim_shape,shape_size,line_width,plot_sim,dim)

sim=curve.sim(:,arg_ind);
meas=curve.meas(:,arg_ind);
sd=curve.sd(:,arg_ind);
simA_ini=curve.sim0;

z=curve.z(arg_ind);
%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%

yu=meas+sd;
yl=meas-sd;

x=curve.x(:,arg_ind);% x must be a column

yu=yu';
yl=yl'; % yl AND yu must be rows.

[~,idx] = sort(x);
x=x(idx);
yu=yu(idx);
yl=yl(idx);
measured=meas(idx);
simulated=sim(idx);
simulated0=simA_ini(idx);

blue_c=[0 0.4470 0.7410];
green_c=[0.4660 0.6740 0.1880]; % Green
orange_c=[0.8500, 0.3250, 0.0980];

xfill=[x',flipud(x)']; % xfill must be a single row: the second curve is reversed. It should be like [1 2 3 3 2 1] It fills between two curves. It must be in the same row, and
yfill=[yl, fliplr(yu)];

if year==2022
    meas_c=blue_c;
else
    meas_c=orange_c;
end

if year==2023 && y_axis=="AQ"
    plot(x,measured,'o','color',blue_c,'MarkerSize',shape_size-2, 'linewidth',line_width)
    hold on;
    fill(xfill,yfill,blue_c,'LineStyle','none','FaceAlpha',0.2);
end

plot(x,measured,'o','color',meas_c,'MarkerSize',shape_size-2, 'linewidth',line_width)
hold on;
fill(xfill,yfill,meas_c,'LineStyle','none','FaceAlpha',.2);


hold on 
plot(x,simulated0,'*','color',green_c,'MarkerSize',shape_size, 'linewidth',line_width)

if plot_sim
    hold on 
    plot(x,simulated,sim_shape,'color',sim_color,'MarkerSize',shape_size, 'linewidth',line_width)
end

data=load("../data/processed_data/final_acc22_23.mat");
final_acc=data.final_acc;
acc_lab=replace(final_acc(arg_ind),"_","-");
% title({[strcat(acc_lab,", GS:",string(round(z,2)))],[strcat("MCMC:",string(round(z_S,2)))]},FontSize=title_size) 


str = strcat("X^2 = ",string(round(z,2)));
a=annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.FontSize = 16;

if strcmp(y_axis,"ACI")
    if max(simulated)>1
        ylim([0,60])
    else
        ylim([0,0.4])
    end
    
    xlim([0,1250])
elseif strcmp(y_axis,"AQ")
    ylim([0,40])
    xlim([0,1850])
end