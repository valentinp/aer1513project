clear all; close all; clc;
load('msckf_NST_on.mat');
load('msckf_NST_off.mat');

msckf_trans_err = p_C_G_err;
msckf_rot_err = theta_CG_err;
swf_trans_err = p_C_G_err_NSToff;
swf_rot_err = theta_CG_err_NSToff;

transLim = 0.4;
rotLim = 0.4;
fontSize = 14;

figure
subplot(3,1,1)
plot(tPlot, msckf_trans_err(1,:), '-b', 'LineWidth', 1.2)
hold on
plot(tPlot, swf_trans_err(1,:), '-r', 'LineWidth', 1.2)
xlim([tPlot(1) tPlot(end)]);
ylim([-transLim transLim])
legend('With Nullspace Projection', 'Without Nullspace Projection');
title(sprintf('Translational Error'))
ylabel('\delta r_x [m]')
set(gca,'FontSize',fontSize)
grid on; grid minor;

subplot(3,1,2)
plot(tPlot, msckf_trans_err(2,:), '-b', 'LineWidth', 1.2)
hold on
plot(tPlot, swf_trans_err(2,:), '-r', 'LineWidth', 1.2)
xlim([tPlot(1) tPlot(end)]);
ylim([-transLim transLim])
ylabel('\delta r_y [m]')
set(gca,'FontSize',fontSize)
grid on; grid minor;

subplot(3,1,3)
plot(tPlot, msckf_trans_err(3,:), '-b', 'LineWidth', 1.2)
hold on
plot(tPlot, swf_trans_err(3,:), '-r', 'LineWidth', 1.2)
xlim([tPlot(1) tPlot(end)]);
ylim([-transLim transLim])
ylabel('\delta r_z [m]')
xlabel('t_k [s]')
%set(gca,'FontSize',12)
%set(findall(gcf,'type','text'),'FontSize',12)
set(gca,'FontSize',fontSize)
grid on; grid minor;

figure
subplot(3,1,1)
plot(tPlot, msckf_rot_err(1,:), '-b', 'LineWidth', 1.2)
hold on
plot(tPlot, swf_rot_err(1,:), '-r', 'LineWidth', 1.2)
xlim([tPlot(1) tPlot(end)]);
ylim([-rotLim rotLim])
title(sprintf('Rotational Error'))
ylabel('\delta\theta_x')
legend('With Nullspace Projection', 'Without Nullspace Projection');
set(gca,'FontSize',fontSize)
grid on; grid minor;

 
subplot(3,1,2)
plot(tPlot, msckf_rot_err(2,:), '-b', 'LineWidth', 1.2)
hold on
plot(tPlot, swf_rot_err(3,:), '-r', 'LineWidth', 1.2)
xlim([tPlot(1) tPlot(end)]);
ylim([-rotLim rotLim])
ylabel('\delta\theta_y')
set(gca,'FontSize',fontSize)
grid on; grid minor;

subplot(3,1,3)
plot(tPlot, msckf_rot_err(3,:), '-b', 'LineWidth', 1.2)
hold on
plot(tPlot, swf_rot_err(3,:), '-r', 'LineWidth', 1.2)
xlim([tPlot(1) tPlot(end)]);
ylim([-rotLim rotLim])
ylabel('\delta\theta_z')
xlabel('t_k [s]')
set(gca,'FontSize',fontSize)
grid on; grid minor;