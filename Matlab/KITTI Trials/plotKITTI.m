%% Load data
clear; close all; clc;

fileName = '2011_09_26_drive_0001';
% fileName = '2011_09_26_drive_0051';
% fileName = '2011_09_26_drive_0095';
% fileName = '2011_09_26_drive_0036';

fileSuffix = '_sync_KLT';

swfData = load(sprintf('SWF_%s%s',fileName,fileSuffix));
msckfData = load(sprintf('msckf_%s%s',fileName,fileSuffix));

%% Compute RMSE
swf_trans_rmse = sqrt(mean(swfData.transErrVec.^2,1));
swf_rot_rmse = sqrt(mean(swfData.rotErrVec.^2,1));

msckf_trans_rmse = sqrt(mean(msckfData.msckf_trans_err.^2,1));
msckf_rot_rmse = sqrt(mean(msckfData.msckf_rot_err.^2,1));

imu_trans_rmse = sqrt(mean(msckfData.imu_trans_err.^2,1));
imu_rot_rmse = sqrt(mean(msckfData.imu_rot_err.^2,1));

%% Plot stuff
figure(1); clf;
fontSize = 16;
lineWidth = 2;
pos = [200,200,640,480];
xLim = [1, min([size(msckf_rot_rmse,2),size(swf_rot_rmse,2)])];

% Translational RMSE
subplot(2,1,1);
plot(imu_trans_rmse,'-k','LineWidth',lineWidth); hold on;
plot(swf_trans_rmse,'-r','LineWidth',lineWidth); 
plot(msckf_trans_rmse,'-b','LineWidth',lineWidth);
xlim(xLim);
title(fileName,'Interpreter','none');
ylabel('Trans. RMSE (m)');
legend('IMU Only','SWF','MSCKF','Location','NorthWest');
grid minor; box on;
set(gca,'FontSize',fontSize);
set(gcf,'Position',pos);

% Rotational RMSE
subplot(2,1,2);
plot(imu_rot_rmse,'-k','LineWidth',lineWidth); hold on;
plot(swf_rot_rmse,'-r','LineWidth',lineWidth); 
plot(msckf_rot_rmse,'-b','LineWidth',lineWidth);
xlim(xLim);
xlabel('Timestep'); ylabel('Rot. RMSE (Axis-Angle)');
legend('IMU Only','SWF','MSCKF','Location','NorthWest');
grid minor; box on;
set(gca,'FontSize',fontSize);
set(gcf,'Position',pos);

%% Export figure
figFileName = ['RMSE_',fileName,'.pdf'];
export_fig(gcf, figFileName, '-transparent');