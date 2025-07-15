function  [] = MF1_plot_curves_individual_acc(curve,y_axis,arg_ind,plotset)


sim_color=plotset.sim_color;
meas_color=plotset.meas_color;
sim_shape=plotset.sim_shape;
shape_size=plotset.sim_size;
line_width=plotset.line_width;
% dim=plotset.annot_dim;
% annot_size=plotset.annot_size;
title_size=plotset.title_size;
title_text=plotset.title_text;

plot_measured=plotset.plot_measured;
plot_simulated=plotset.plot_simulated;

sim=curve.sim(:,arg_ind);
meas=curve.meas(:,arg_ind);
sd=curve.sd(:,arg_ind);
% simA_ini=curve.sim0(:,arg_ind);

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
% simulated0=simA_ini(idx);




xfill=[x',flipud(x)']; % xfill must be a single row: the second curve is reversed. It should be like [1 2 3 3 2 1] It fills between two curves. It must be in the same row, and
yfill=[yl, fliplr(yu)];

if plot_measured==1
    fill(xfill,yfill,meas_color,'LineStyle','none','FaceAlpha',.2);
    hold on;
    plot(x,measured,'o','color',meas_color,'MarkerSize',shape_size-2, 'linewidth',line_width)
    hold on 
end
% hold on 
% plot(x,simulated0,'*','color',green_c,'MarkerSize',shape_size, 'linewidth',line_width)
if plot_simulated==1

    plot(x,simulated,sim_shape,'color',sim_color,'MarkerSize',shape_size, 'linewidth',line_width)
end

title(title_text,'FontSize', title_size,'FontWeight','bold','Color',plotset.title_color) 


% str = strcat("X^2 = ",string(round(z,3)));
% a=annotation('textbox',dim,'String',str, 'EdgeColor','none','FitBoxToText','on'); %
% a.FontSize = annot_size;

if strcmp(y_axis,"ACI")
    if max(simulated)>1
        ylim([0,60])
    else
        ylim([0, 1])
    end
    xlim([0,1300])
elseif strcmp(y_axis,"AQ")
    ylim([0,45])
    xlim([0,1850])
end