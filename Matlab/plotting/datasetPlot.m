% load('../datasets/dataset3.mat');
load('../datasets/dataset3_fresh_100lessnoisy.mat');

fontSize = 14;
markerSize = 50;
lineWidth = 1.2;

figure(1); clf; hold on;

scatter3(rho_i_pj_i(1,:), rho_i_pj_i(2,:), rho_i_pj_i(3,:), markerSize, 'or', 'fill');
plot3(r_i_vk_i(1,:), r_i_vk_i(2,:), r_i_vk_i(3,:), '-b', 'LineWidth', lineWidth);
legend('Features', 'Sensor trajectory');

set(gca, 'FontSize', fontSize);
campos([-18.7401, -5.8668, 10.6855]);
camtarget([1.5, 1.5, 1]);
grid on; grid minor;
axis equal;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');