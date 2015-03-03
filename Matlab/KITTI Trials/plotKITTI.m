%% Load data
clear; close all; clc;

% fileName = '2011_09_26_drive_0001';
% fileName = '2011_09_26_drive_0051';
% fileName = '2011_09_26_drive_0095';
fileName = '2011_09_26_drive_0036';

fileSuffix = '_sync_KLT';

swfData = load(sprintf('SWF_%s%s',fileName,fileSuffix));
msckfData = load(sprintf('msckf_%s%s',fileName,fileSuffix));

dists = [0,cumsum(sqrt(sum(diff(msckfData.p_I_G_GT,1,2).^2,1)))];
swfData.transErrVec = swfData.transErrVec(:,1:size(dists,2));
swfData.rotErrVec = swfData.rotErrVec(:,1:size(dists,2));

%% Compute RMSE
swf_trans_rmse = sqrt(mean(swfData.transErrVec.^2,1));
swf_rot_rmse = sqrt(mean(swfData.rotErrVec.^2,1));

msckf_trans_rmse = sqrt(mean(msckfData.msckf_trans_err.^2,1));
msckf_rot_rmse = sqrt(mean(msckfData.msckf_rot_err.^2,1));

imu_trans_rmse = sqrt(mean(msckfData.imu_trans_err.^2,1));
imu_rot_rmse = sqrt(mean(msckfData.imu_rot_err.^2,1));

%% Compute ANEES
% swf_anees = ANEES(swfData.transErrVec,swfData.rotErrVec,
msckf_anees = ANEES(msckfData.msckf_trans_err, msckfData.msckf_rot_err,...
    msckfData.err_sigma);

%% Plot stuff
figure(1); clf;
fontSize = 14;
lineWidth = 2;
pos = [200,200,640,320];
xLim = [dists(1), dists(end)];

% Translational RMSE
subplot(2,1,1);
plot(dists,imu_trans_rmse,'-k','LineWidth',lineWidth); hold on;
plot(dists,swf_trans_rmse,'-r','LineWidth',lineWidth); 
plot(dists,msckf_trans_rmse,'-b','LineWidth',lineWidth);
xlim(xLim);
title(fileName,'Interpreter','none');
ylabel('Trans. RMSE (m)');
legend('IMU Only','SWF','MSCKF','Location','NorthWest');
grid minor; box on;
set(gca,'FontSize',fontSize);
set(gcf,'Position',pos);

% Rotational RMSE
subplot(2,1,2);
plot(dists,imu_rot_rmse,'-k','LineWidth',lineWidth); hold on;
plot(dists,swf_rot_rmse,'-r','LineWidth',lineWidth); 
plot(dists,msckf_rot_rmse,'-b','LineWidth',lineWidth);
xlim(xLim);
xlabel('Distance Travelled (m)'); ylabel('Rot. RMSE (Axis-Angle)');
legend('IMU Only','SWF','MSCKF','Location','NorthWest');
grid minor; box on;
set(gca,'FontSize',fontSize);
set(gcf,'Position',pos);

%% Export figure
figFileName = ['RMSE_',fileName,'.pdf'];
export_fig(gcf, figFileName, '-nocrop','-transparent');

%% Stats
fprintf('Total distance: %f \n', dists(end));
fprintf('MSCKF:\n\t trans armse %f \n\t rot armse %f \n\t anees %f',...
    mean(msckf_trans_rmse),mean(msckf_rot_rmse),msckf_anees);
fprintf('SF:\n\t trans armse %f \n\t rot armse %f \n\t anees %f',...
    mean(swf_trans_rmse),mean(swf_rot_rmse),swf_anees);