%% Test feature Gauss Newton triangulation
%Generate a forward trajectory
i = 1;
T_wCam_GT = [];
for t = 0:0.2:2
    T_wCam_GT(:,:,i) = [eye(3) [t 0 0]'; 0 0 0 1];
    i=i+1;
end
%Generate the true landmarks
landmarks_w = [1 0 5]';
%for x_i = 0:0.1:2
%     for y_i = -1:0.1:1
%         landmarks_w = [landmarks_w [x_i y_i 5]'];
%     end
% end
simSetup.pixelNoiseStd = 1; %pixels
simSetup.cameraResolution = [1280, 960]; %pixels

%Set the camera intrinsics
pixelWidth = 4.8*1e-6;
focalLength = 4*1e-3/pixelWidth; 
c_u = 640;
c_v = 480;

K  = [focalLength 0 c_u;
    0 focalLength c_v;
    0 0 1];

invK = inv(K);

imageMeasurements = genFeatureMeasurements(T_wCam_GT, landmarks_w, K, simSetup)



%[p_f_G] = calcGNPosEst(camStates, observations, camera)
