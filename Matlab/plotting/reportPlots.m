plotAbsErr = true;

%% 500-1000

swf10 = load('swf_500_1000_10.mat');
swf50 = load('swf_500_1000_50.mat');
swf100 = load('swf_500_1000_100.mat');
msckf5 = load('msckf_500_1000_min5_maxInf');
msckf10 = load('msckf_500_1000_min10_max50');
msckf20 = load('msckf_500_1000_min20_max100');
imu = load('imu_500_1000.mat');

load('../datasets/dataset3.mat');
kStart = 500;
kEnd = 1000;


transLim = [-0.5 0.5];
rotLim = [-0.25 0.25];
fontSize = 14;

if plotAbsErr
    imu.msckf_trans_err = abs(imu.msckf_trans_err);
    msckf5.msckf_trans_err = abs(msckf5.msckf_trans_err);
    msckf10.msckf_trans_err = abs(msckf10.msckf_trans_err);
    msckf20.msckf_trans_err = abs(msckf20.msckf_trans_err);
    swf10.swf_trans_err = abs(swf10.swf_trans_err);
    swf50.swf_trans_err = abs(swf50.swf_trans_err);
    swf100.swf_trans_err = abs(swf100.swf_trans_err);

    imu.msckf_rot_err = abs(imu.msckf_rot_err);
    msckf5.msckf_rot_err = abs(msckf5.msckf_rot_err);
    msckf10.msckf_rot_err = abs(msckf10.msckf_rot_err);
    msckf20.msckf_rot_err = abs(msckf20.msckf_rot_err);
    swf10.swf_rot_err = abs(swf10.swf_rot_err);
    swf50.swf_rot_err = abs(swf50.swf_rot_err);
    swf100.swf_rot_err = abs(swf100.swf_rot_err);
    
    transLim = transLim - transLim(1);
    rotLim = rotLim - rotLim(1);
end

figure
subplot(3,1,1)
plot(t(kStart:kEnd), imu.msckf_trans_err(1,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_trans_err(1,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_trans_err(1,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_trans_err(1,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_trans_err(1,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_trans_err(1,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_trans_err(1,:), '--g', 'LineWidth', 1.2)

xlim([t(kStart) t(kEnd) ]);
ylim(transLim)
h_legend = legend('IMU Only','MSCKF 5-Inf', 'MSCKF 10-50', 'MSCKF 20-100', 'SWF 10', 'SWF 50', 'SWF 100','Location', 'northwest');
set(h_legend,'FontSize',10);

title(sprintf('Absolute Translational Error (Interval 1)'))
ylabel('\delta r_x [m]')
set(gca,'FontSize',fontSize)
grid on
grid minor
box on

subplot(3,1,2)
plot(t(kStart:kEnd), imu.msckf_trans_err(2,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_trans_err(2,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_trans_err(2,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_trans_err(2,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_trans_err(2,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_trans_err(2,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_trans_err(2,:), '--g', 'LineWidth', 1.2)

xlim([t(kStart) t(kEnd) ]);
ylim(transLim)
ylabel('\delta r_y [m]')
set(gca,'FontSize',fontSize)
grid on
grid minor
box on

subplot(3,1,3)
plot(t(kStart:kEnd), imu.msckf_trans_err(3,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_trans_err(3,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_trans_err(3,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_trans_err(3,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_trans_err(3,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_trans_err(3,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_trans_err(3,:), '--g', 'LineWidth', 1.2)

xlim([t(kStart) t(kEnd) ]);
ylim(transLim)
ylabel('\delta r_z [m]')
xlabel('t_k [s]')
set(gca,'FontSize',fontSize)
set(findall(gcf,'type','text'),'FontSize',fontSize)
grid on
grid minor
box on

filename = sprintf('6-Way-Comparison-500-1000-Trans.pdf');
export_fig(gcf, filename, '-transparent');

figure
subplot(3,1,1)
plot(t(kStart:kEnd), imu.msckf_rot_err(1,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_rot_err(1,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_rot_err(1,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_rot_err(1,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_rot_err(1,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_rot_err(1,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_rot_err(1,:), '--g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd) ]);
ylim(rotLim)
title(sprintf('Absolute Rotational Error (Interval 1)'))
ylabel('\delta\theta_x')
h_legend = legend('IMU Only','MSCKF 5-Inf', 'MSCKF 10-50', 'MSCKF 20-100', 'SWF 10', 'SWF 50', 'SWF 100','Location', 'northwest');
set(h_legend,'FontSize',10);
set(gca,'FontSize',fontSize)
grid on
grid minor
box on

 
subplot(3,1,2)
plot(t(kStart:kEnd), imu.msckf_rot_err(2,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_rot_err(2,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_rot_err(2,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_rot_err(2,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_rot_err(2,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_rot_err(2,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_rot_err(2,:), '--g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd) ]);
ylim(rotLim)
ylabel('\delta\theta_y')
set(gca,'FontSize',fontSize)
grid on
grid minor
box on

subplot(3,1,3)
plot(t(kStart:kEnd), imu.msckf_rot_err(3,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_rot_err(3,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_rot_err(3,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_rot_err(3,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_rot_err(3,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_rot_err(3,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_rot_err(3,:), '--g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd) ]);
ylim(rotLim)
ylabel('\delta\theta_z')
xlabel('t_k [s]')
grid on
grid minor
box on

set(gca,'FontSize',fontSize)
set(findall(gcf,'type','text'),'FontSize',fontSize)

filename = sprintf('6-Way-Comparison-500-1000-Rot.pdf');
export_fig(gcf, filename, '-transparent');

%% 1215-1715

swf10 = load('swf_1215_1715_10.mat');
swf50 = load('swf_1215_1715_50.mat');
swf100 = load('swf_1215_1715_100.mat');
msckf5 = load('msckf_1215_1715_min5_maxInf');
msckf10 = load('msckf_1215_1715_min10_max50');
msckf20 = load('msckf_1215_1715_min20_max100');
imu = load('imu_1215_1715.mat');

load('../datasets/dataset3.mat');
kStart = 1215;
kEnd = 1715;

transLim = [-1 1];
rotLim = [-0.5 0.5];
fontSize = 14;

if plotAbsErr
    imu.msckf_trans_err = abs(imu.msckf_trans_err);
    msckf5.msckf_trans_err = abs(msckf5.msckf_trans_err);
    msckf10.msckf_trans_err = abs(msckf10.msckf_trans_err);
    msckf20.msckf_trans_err = abs(msckf20.msckf_trans_err);
    swf10.swf_trans_err = abs(swf10.swf_trans_err);
    swf50.swf_trans_err = abs(swf50.swf_trans_err);
    swf100.swf_trans_err = abs(swf100.swf_trans_err);

    imu.msckf_rot_err = abs(imu.msckf_rot_err);
    msckf5.msckf_rot_err = abs(msckf5.msckf_rot_err);
    msckf10.msckf_rot_err = abs(msckf10.msckf_rot_err);
    msckf20.msckf_rot_err = abs(msckf20.msckf_rot_err);
    swf10.swf_rot_err = abs(swf10.swf_rot_err);
    swf50.swf_rot_err = abs(swf50.swf_rot_err);
    swf100.swf_rot_err = abs(swf100.swf_rot_err);
    
    transLim = transLim - transLim(1);
    rotLim = rotLim - rotLim(1);
end

figure
subplot(3,1,1)
plot(t(kStart:kEnd), imu.msckf_trans_err(1,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_trans_err(1,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_trans_err(1,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_trans_err(1,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_trans_err(1,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_trans_err(1,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_trans_err(1,:), '--g', 'LineWidth', 1.2)

xlim([t(kStart) t(kEnd) ]);
ylim(transLim)
h_legend = legend('IMU Only','MSCKF 5-Inf', 'MSCKF 10-50', 'MSCKF 20-100', 'SWF 10', 'SWF 50', 'SWF 100','Location', 'northwest');
set(h_legend,'FontSize',10);

title(sprintf('Absolute Translational Error (Interval 2)'))
ylabel('\delta r_x [m]')
set(gca,'FontSize',fontSize)
grid on
grid minor
box on

subplot(3,1,2)
plot(t(kStart:kEnd), imu.msckf_trans_err(2,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_trans_err(2,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_trans_err(2,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_trans_err(2,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_trans_err(2,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_trans_err(2,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_trans_err(2,:), '--g', 'LineWidth', 1.2)

xlim([t(kStart) t(kEnd) ]);
ylim(transLim)
ylabel('\delta r_y [m]')
set(gca,'FontSize',fontSize)
grid on
grid minor
box on


subplot(3,1,3)
plot(t(kStart:kEnd), imu.msckf_trans_err(3,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_trans_err(3,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_trans_err(3,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_trans_err(3,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_trans_err(3,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_trans_err(3,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_trans_err(3,:), '--g', 'LineWidth', 1.2)

xlim([t(kStart) t(kEnd) ]);
ylim(transLim)
ylabel('\delta r_z [m]')
xlabel('t_k [s]')
set(gca,'FontSize',fontSize)
set(findall(gcf,'type','text'),'FontSize',fontSize)
grid on
grid minor
box on

filename = sprintf('6-Way-Comparison-1215-1715-Trans.pdf');
export_fig(gcf, filename, '-transparent');

figure
subplot(3,1,1)
plot(t(kStart:kEnd), imu.msckf_rot_err(1,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_rot_err(1,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_rot_err(1,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_rot_err(1,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_rot_err(1,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_rot_err(1,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_rot_err(1,:), '--g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd) ]);
ylim(rotLim)
title(sprintf('Absolute Rotational Error (Interval 2)'))
ylabel('\delta\theta_x')
h_legend = legend('IMU Only','MSCKF 5-Inf', 'MSCKF 10-50', 'MSCKF 20-100', 'SWF 10', 'SWF 50', 'SWF 100','Location', 'northwest');
set(h_legend,'FontSize',10);
set(gca,'FontSize',fontSize)
grid on
grid minor
box on

 
subplot(3,1,2)
plot(t(kStart:kEnd), imu.msckf_rot_err(2,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_rot_err(2,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_rot_err(2,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_rot_err(2,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_rot_err(2,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_rot_err(2,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_rot_err(2,:), '--g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd) ]);
ylim(rotLim)
ylabel('\delta\theta_y')
set(gca,'FontSize',fontSize)
grid on
grid minor
box on

subplot(3,1,3)
plot(t(kStart:kEnd), imu.msckf_rot_err(3,:), '-k', 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), msckf5.msckf_rot_err(3,:), '-b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf10.msckf_rot_err(3,:), '-.b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), msckf20.msckf_rot_err(3,:), '--b', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf10.swf_rot_err(3,:), '-g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf50.swf_rot_err(3,:), '-.g', 'LineWidth', 1.2)
plot(t(kStart:kEnd), swf100.swf_rot_err(3,:), '--g', 'LineWidth', 1.2)
xlim([t(kStart) t(kEnd) ]);
ylim(rotLim)
ylabel('\delta\theta_z')
xlabel('t_k [s]')
grid on
grid minor
box on

set(gca,'FontSize',fontSize)
set(findall(gcf,'type','text'),'FontSize',fontSize)

filename = sprintf('6-Way-Comparison-1215-1715-Rot.pdf');
export_fig(gcf, filename, '-transparent');

%% Fresh2 dataset

msckf_1215_1715_min20_max100_fresh2
swf_1215_1715_50_fresh2
