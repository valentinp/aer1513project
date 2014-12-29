load('swf_est25.mat');
load('msckf_est.mat');
load('datasets/dataset3.mat');
addpath('/Users/valentinp/Research/MATLAB/export_fig'); %Use Oliver Woodford's awesome export_fig package to get trimmed PDFs
kStart = 500;
kEnd = 1000;
swf_trans_err = transErrVec;
swf_rot_err = rotErrVec;

transLim = 0.4;
rotLim = 0.4;
fontSize = 14;

figure
subplot(3,1,1)
plot(t(kStart:kEnd), msckf_trans_err(1,:), '-b', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), swf_trans_err(1,:), '-g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd)]);
ylim([-transLim transLim])
legend('MSCKF', 'SWF');
title(sprintf('Translational Error (MSCKF vs. SWF)'))
ylabel('\delta r_x [m]')
set(gca,'FontSize',fontSize)
grid on

subplot(3,1,2)
plot(t(kStart:kEnd), msckf_trans_err(2,:), '-b', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), swf_trans_err(2,:), '-g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd)]);
ylim([-transLim transLim])
ylabel('\delta r_y [m]')
set(gca,'FontSize',fontSize)
grid on


subplot(3,1,3)
plot(t(kStart:kEnd), msckf_trans_err(3,:), '-b', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), swf_trans_err(3,:), '-g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd)]);
ylim([-transLim transLim])
ylabel('\delta r_z [m]')
xlabel('t_k [s]')
set(gca,'FontSize',fontSize)
set(findall(gcf,'type','text'),'FontSize',fontSize)
grid on

filename = sprintf('transError.pdf');
export_fig(gcf, filename, '-transparent');

figure
subplot(3,1,1)
plot(t(kStart:kEnd), msckf_rot_err(1,:), '-b', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), swf_rot_err(1,:), '-g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd)]);
ylim([-rotLim rotLim])
title(sprintf('Rotational Error (MSCKF vs. SWF)'))
ylabel('\delta\theta_x')
legend('MSCKF', 'SWF');
set(gca,'FontSize',fontSize)
grid on


 
subplot(3,1,2)
plot(t(kStart:kEnd), msckf_rot_err(2,:), '-b', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), swf_rot_err(3,:), '-g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd)]);
ylim([-rotLim rotLim])
ylabel('\delta\theta_y')
set(gca,'FontSize',fontSize)
grid on

subplot(3,1,3)
plot(t(kStart:kEnd), msckf_rot_err(3,:), '-b', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), swf_rot_err(3,:), '-g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd)]);
ylim([-rotLim rotLim])
ylabel('\delta\theta_z')
xlabel('t_k [s]')
grid on

set(gca,'FontSize',fontSize)
set(findall(gcf,'type','text'),'FontSize',fontSize)

filename = sprintf('rotError.pdf');
export_fig(gcf, filename, '-transparent');