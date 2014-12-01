addpath('utils');
load('dataset3.mat')

%Set up the camera parameters
camera.c_u = cu;
camera.c_v = cv;
camera.f_u = fu;
camera.f_v = fv;
camera.b = b;

% Notation: X_sub_super

imuState.q_iG = [zeros(3,1); 1];  %Global to IMU rotation quaternion
imuState.p_G = zeros(3,1); %IMU Position in the Global frame
imuState.b_g = zeros(3,1); %Gyro bias
imuState.b_v = zeros(3,1); %Velocity bias

camState.q_CG = [zeros(3,1); 1];
camState.p_C_G = zeros(3,1);

%msckfState = {imuState, camState1, camState2,....}
