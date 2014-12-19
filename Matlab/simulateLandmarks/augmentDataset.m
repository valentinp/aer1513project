clc;
clear all;
addpath('../msckf/utils');
load('../dataset3.mat')

%Set up appropriate structs
calibParams.c_u = cu;
calibParams.c_v = cv;
calibParams.f_u = fu;
calibParams.f_v = fv;
calibParams.b = b;

%Generate new landmarks
newLmNum = 40;
xRange = [min(rho_i_pj_i(1,:)) max(rho_i_pj_i(1,:))];
yRange = [min(rho_i_pj_i(2,:)) max(rho_i_pj_i(2,:))];
zRange = [min(rho_i_pj_i(3,:)) max(rho_i_pj_i(3,:))];

newLmPosX = range(xRange)*rand(1,newLmNum) + xRange(1);
newLmPosY = range(yRange)*rand(1,newLmNum) + yRange(1);
newLmPosZ = range(zRange)*rand(1,newLmNum) + zRange(1);

newLMPos = [newLmPosX; newLmPosY; newLmPosZ];

rho_i_pj_i = [rho_i_pj_i newLMPos];
T_cv = [C_c_v -C_c_v*rho_v_c_v; 0 0 0 1];

for k = 1:length(t)
    C_vi = axisAngleToRotMat(theta_vk_i(:,k));
    T_vi = [C_vi -C_vi*r_i_vk_i(:,k); 0 0 0 1];
    T_ci = T_cv*T_vi;
    
    %Add observations
    for lm_i = 1:size(newLMPos, 2)
        p_li_i = newLMPos(:,lm_i);
        p_lc_c = homo2cart(T_ci*cart2homo(p_li_i));
        [yMeas] = stereoCamProject(p_lc_c, calibParams);
        if all(yMeas > 0) && all(yMeas([1,3]) < 700) && all(yMeas([2,4]) < 500)
                y_k_j(:,k, 20+lm_i) = yMeas + randn(4,1);
        else
                y_k_j(:,k, 20+lm_i) = -1*ones(4,1);
        end
    end
end


save('../dataset3_augmented.mat');