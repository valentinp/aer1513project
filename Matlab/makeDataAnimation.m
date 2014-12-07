load('dataset3');

fontSize = 16;

writerObj = VideoWriter('dataset3movie.mp4', 'MPEG-4');
writerObj.FrameRate=60;

open(writerObj);

h = figure(1); clf; hold on;
set(h, 'Position', [0 0 640 480], 'Color', 'w');

scatter3(rho_i_pj_i(1,:), rho_i_pj_i(2,:), rho_i_pj_i(3,:), 'or', 'fill');
xlim([0,3]); ylim([0,3]); zlim([0,2]);
set(gca, 'FontSize', fontSize);
campos([-18.7401, -5.8668, 10.6855]);
camtarget([1.5, 1.5, 1]);
grid on;

hplot = plot3(r_i_vk_i(1,1), r_i_vk_i(2,1), r_i_vk_i(3,1), '-b');
writeVideo(writerObj, getframe(h));
for i = 2:numel(t)
    set(hplot, 'XData', r_i_vk_i(1,1:i), 'YData', r_i_vk_i(2,1:i), 'ZData', r_i_vk_i(3,1:i));
    drawnow;
    writeVideo(writerObj, getframe(h));
end

close(writerObj);