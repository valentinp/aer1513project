swf20 = load('swf_1215_1715_50_20lessnoisy.mat');
swf40 = load('swf_1215_1715_50_40lessnoisy.mat');
swf60 = load('swf_1215_1715_50_60lessnoisy.mat');
msckf20 = load('msckf_1215_1715_min20_max100_20lessnoisy.mat');
msckf40 = load('msckf_1215_1715_min20_max100_40lessnoisy.mat');
msckf60 = load('msckf_1215_1715_min20_max100_60lessnoisy.mat');
imu = load('imu_1215_1715.mat');


%Calculate Average RMSE (Root-Mean-Squared Error)
clc
disp('IMU Only RMSE')
mean(sqrt(sum(imu.msckf_trans_err.^2, 1)/3))
mean(sqrt(sum(imu.msckf_rot_err.^2, 1)/3))

disp('MSCKF RMSE')
mean(sqrt(sum(msckf20.msckf_trans_err.^2, 1)/3))
mean(sqrt(sum(msckf40.msckf_trans_err.^2, 1)/3))
mean(sqrt(sum(msckf60.msckf_trans_err.^2, 1)/3))
disp('===============')
mean(sqrt(sum(msckf20.msckf_rot_err.^2, 1)/3))
mean(sqrt(sum(msckf40.msckf_rot_err.^2, 1)/3))
mean(sqrt(sum(msckf60.msckf_rot_err.^2, 1)/3))

disp('SWF RMSE')
mean(sqrt(sum(swf20.swf_trans_err.^2, 1)/3))
mean(sqrt(sum(swf40.swf_trans_err.^2, 1)/3))
mean(sqrt(sum(swf60.swf_trans_err.^2, 1)/3))
disp('===============')
mean(sqrt(sum(swf20.swf_rot_err.^2, 1)/3))
mean(sqrt(sum(swf40.swf_rot_err.^2, 1)/3))
mean(sqrt(sum(swf60.swf_rot_err.^2, 1)/3))

%Calculate Average NEES