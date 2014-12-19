clc;
clear all;
addpath('../msckf/utils');
load('../dataset3.mat')

rng(42);

%Set up appropriate structs
calibParams.c_u = cu;
calibParams.c_v = cv;
calibParams.f_u = fu;
calibParams.f_v = fv;
calibParams.b = b;

%Generate new landmarks
newLmNum = 200;
% xRange = linspace(min(rho_i_pj_i(1,:)), max(rho_i_pj_i(1,:)), newLmNum);
% yRange = linspace(min(rho_i_pj_i(2,:)), max(rho_i_pj_i(2,:)), newLmNum);
% zRange = linspace(min(rho_i_pj_i(3,:)), max(rho_i_pj_i(3,:)), newLmNum);
% [newLmPosX, newLmPosY, newLmPosZ ] = meshgrid(xRange, yRange, zRange);

xRange = [min(rho_i_pj_i(1,:)) - 5, max(rho_i_pj_i(1,:)) + 5];
yRange = [min(rho_i_pj_i(2,:)) - 5, max(rho_i_pj_i(2,:)) + 5];
zRange = [min(rho_i_pj_i(3,:)) - 5, max(rho_i_pj_i(3,:)) - 5];

newLmPosX = range(xRange)*rand(1,newLmNum) + xRange(1);
newLmPosY = range(yRange)*rand(1,newLmNum) + yRange(1);
newLmPosZ = range(zRange)*rand(1,newLmNum) + zRange(1);

newLMPos = [newLmPosX(:)'; newLmPosY(:)'; newLmPosZ(:)'];


rho_i_pj_i = newLMPos;
T_cv = [C_c_v -C_c_v*rho_v_c_v; 0 0 0 1];
y_k_j = [];

for k = 1:length(t)
    C_vi = axisAngleToRotMat(theta_vk_i(:,k));
    T_vi = [C_vi -C_vi*r_i_vk_i(:,k); 0 0 0 1];
    T_ci = T_cv*T_vi;
    
    %Add observations
    for lm_i = 1:size(newLMPos, 2)
        p_li_i = newLMPos(:,lm_i);
        p_lc_c = homo2cart(T_ci*cart2homo(p_li_i));
        [yMeas] = stereoCamProject(p_lc_c, calibParams);
        if all(yMeas > 0) && all(yMeas([1,3]) <= 640) && all(yMeas([2,4]) <= 480)
                y_k_j(:,k, lm_i) = yMeas;
        else
                y_k_j(:,k, lm_i) = -1*ones(4,1);
        end
    end
end


save('../dataset3_fresh2.mat');