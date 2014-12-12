function [msckfState, featureTracks, trackedFeatureIds] = initializeMSCKF(firstImuState, firstMeasurements, camera, state_k)
%INITIALIZEMSCKF Initialize the MSCKF with tracked features and ground
%truth


%Compute the first state
firstImuState.b_g = zeros(3,1);
firstImuState.b_v = zeros(3,1);
firstImuState.covar = eye(12);
msckfState.imuState = firstImuState;
msckfState.imuCovar = 0.01*eye(12);


% Compute camera pose from current IMU pose
camState = {};
C_IG = quatToRotMat(msckfState.imuState.q_IG);
camState.q_CG = quatLeftComp(camera.q_CI) * msckfState.imuState.q_IG;
camState.p_C_G = msckfState.imuState.p_I_G + C_IG' * camera.p_C_I;
camState.trackedFeatureIds = [];
camState.state_k = state_k;

msckfState.camStates = {};
msckfState.camStates{1} = camState;

msckfState.camCovar = 0.1*eye(6);
msckfState.imuCamCovar = zeros(12, 6);

%Compute all of the relevant feature tracks
featureTracks = {};
trackedFeatureIds = [];

 for featureId = 1:20
        meas_k = firstMeasurements.y(:, featureId);
        if meas_k(1,1) ~= -1
                %Track new feature
                track.featureId = featureId;
                track.observations = meas_k;
                featureTracks{end+1} = track;
                trackedFeatureIds(end+1) = featureId;
                %Add observation to current camera
                msckfState.camStates{end}.trackedFeatureIds(end+1) = featureId;
        end
 end
 
end

