addpath('utils');
load('dataset3.mat')

%Dataset window bounds
k1 = 1215;
k2 = 1714;

%Set up the camera parameters
camera.c_u      = cu;
camera.c_v      = cv;
camera.f_u      = fu;
camera.f_v      = fv;
camera.b        = b;
camera.q_CI     = rotMatToQuat(C_c_v);
camera.p_C_I    = rho_v_c_v;

% Notation: X_sub_super, q_FromTo, p_ofWhat_expressedInWhatFrame

%Set up the noise parameters
noiseParams.z_1 = 1;
noiseParams.z_2 = 1;


imuState = cell{1,numel(t)};
imuState{1}.q_IG  = [zeros(3,1); 1];    %Global to IMU rotation quaternion
imuState{1}.p_I_G = zeros(3,1);         %IMU Position in the Global frame
imuState{1}.b_g   = zeros(3,1);         %Gyro bias
imuState{1}.b_v   = zeros(3,1);         %Velocity bias
imuState{1}.covar = zeros(12,12);       %IMU state covariance

camState = cell{1,numel(t)};
camStates{1}.q_CG  = [zeros(3,1); 1];
camStates{1}.p_C_G = zeros(3,1);

%msckfState.imuState
%msckfState.imuCovar
%msckfState.camCovar
%msckfState.imuCamCovar
%msckfState.camStates

% Measurements as cells
dT = [0, diff(t)];
measurements = cell{1,numel(t)};
for k = k1:k2 
    measurements{k}.dT    = dT(k);           % sampling times
    measurements{k}.y     = y_k_j(1:2,k,:);  % left camera only
    measurements{k}.omega = w_vk_vk_i(:,k);  % ang vel
    measurements{k}.v     = v_vk_vk_i(:,k);  % lin vel
end

Nmax = 50;


%Propagate state and covariance
msckfState = propagateMsckfCovar(msckfState, measurements_k, noiseParams);

%Add camera pose to msckfState
msckfState = augmentState(msckfState, camera);

%Continue until a feature is no longer seen or there are Nmax camera frames

%Estimate feature 3D location through Gauss Newton inverse depth
%IMPORTANT: Use 'ideal' measurements:
%u' = (u - c_u)/f_u;
%v' = (v - c_v)/f_v;

[p_f_G] = calcGNPosEst(camStates, observations)

%Compute rh, Th
[H_o_j, A] = calcHoj(p_f_G, msckfState, camStateIndex)
[T_H, Q_1] = calcTH(H_o)
[r_n_j] = calcResidual(r__j, A_j, Q_1)

%Calculate Kalman gain
K = P

%Correct covariance

%Remove any camera states with no tracked features