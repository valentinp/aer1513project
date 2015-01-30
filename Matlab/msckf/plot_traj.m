kNum = length(prunedStates);
p_C_G_est = zeros(3,kNum);
q_CG_est = zeros(4,kNum);

for k=1:kNum
    p_C_G_est(:,k) = prunedStates{k}.p_C_G;
    q_CG_est(:,k)  = prunedStates{k}.q_CG;
end

figure(1); clf; hold on;
plot3(p_C_G_est(1,:),p_C_G_est(2,:),p_C_G_est(3,:),'-b');
if ~isempty(map)
    scatter3(map(1,:),map(2,:),map(3,:),'or');
end
xlabel('x');ylabel('y');zlabel('z');
pause(0.001); drawnow;