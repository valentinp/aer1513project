%% ============================Notation============================ %%
% X_sub_super
% q_FromTo
% p_ofWhat_expressedInWhatFrame


%% =============================Setup============================== %%
clear;
close all;
clc
addpath('utils');
load('dataset3.mat')

%Dataset window bounds
kStart = 600;
kEnd = 1714;

%Set up the camera parameters
camera.c_u      = cu;                   % Principal point [pixels]
camera.c_v      = cv;                   % |
camera.f_u      = fu;                   % Focal length [u pixels]
camera.f_v      = fv;                   % Focal length [v pixels]
camera.b        = b;                    % Stereo baseline [m]
camera.q_CI     = rotMatToQuat(C_c_v);  % 4x1 IMU-to-Camera rotation quaternion
camera.p_C_I    = rho_v_c_v;            % 3x1 Camera position in IMU frame

%Set up the noise parameters
noiseParams.z_1 = 1;
noiseParams.z_2 = 1;
noiseParams.Q_imu = eye(12);
noiseParams.imageVariance = mean([noiseParams.z_1, noiseParams.z_2]);  % Slightly hacky. Used to compute the Kalman gain and corrected covariance in the EKF step

%MSCKF parameters
MSCKFParams.minTrackLength = 5;


% IMU state for plotting etc. Structures indexed in a cell array
imuStates = cell(1,numel(t));
% imuStates{k}.q_IG         4x1 Global to IMU rotation quaternion
% imuStates{k}.p_I_G        3x1 IMU Position in the Global frame
% imuStates{k}.b_g          3x1 Gyro bias
% imuStates{k}.b_v          3x1 Velocity bias
% imuStates{k}.covar        12x12 IMU state covariance

% We don't really need these outside of msckfState, do we?
% camState = cell(1,numel(t));
% camStates{k}.q_CG        4x1 Global to camera rotation quaternion
% camStates{k}.p_C_G       3x1 Camera Position in the Global frame
% camStates{k}.trackedFeatureIds  1xM List of feature ids that are currently being tracked from that camera state

%msckfState.imuState
%msckfState.imuCovar
%msckfState.camCovar
%msckfState.imuCamCovar
%msckfState.camStates


% Measurements as structures all indexed in a cell array
dT = [0, diff(t)];
measurements = cell(1,numel(t));
groundTruthStates = cell(1,numel(t));
groundTruthMap = rho_i_pj_i;

for state_k = kStart:kEnd 
    measurements{state_k}.dT    = dT(state_k);                      % sampling times
    measurements{state_k}.y     = squeeze(y_k_j(1:2,state_k,:));    % left camera only
    measurements{state_k}.omega = w_vk_vk_i(:,state_k);             % ang vel
    measurements{state_k}.v     = v_vk_vk_i(:,state_k);             % lin vel
    
    %Idealize measurements
    validMeas = (measurements{state_k}.y(1,:) ~= -1);
    measurements{state_k}.y(1,validMeas) = (measurements{state_k}.y(1,validMeas) - camera.c_u)/camera.f_u;
    measurements{state_k}.y(2,validMeas) = (measurements{state_k}.y(2,validMeas) - camera.c_v)/camera.f_v;
    
    %Ground Truth
    q_IG = rotMatToQuat(axisAngleToRotMat(theta_vk_i(:,state_k)));
    p_I_G = r_i_vk_i(:,state_k);
    
    groundTruthStates{state_k}.imuState.q_IG = q_IG;
    groundTruthStates{state_k}.imuState.p_I_G =p_I_G ;
    
    % Compute camera pose from current IMU pose
    C_IG = quatToRotMat(q_IG);
    q_CG = quatLeftComp(camera.q_CI) * q_IG;
    p_C_G = p_I_G + C_IG' * camera.p_C_I;
    
    groundTruthStates{state_k}.camState.q_CG = q_CG;
    groundTruthStates{state_k}.camState.p_C_G = p_C_G;
    
end

% Other constants
% Currently not used
Nmax = 50;                      % max number of poses before triggering an update

%Struct used to keep track of features
featureTracks = {};
trackedFeatureIds = [];

% featureTracks = {track1, track2, ...}
% track.featureId 
% track.observations



%% ==========================Initial State======================== %%
%Use ground truth for first state and initialize feature tracks with
%feature observations
%Use ground truth for the first state

firstImuState.q_IG = rotMatToQuat(axisAngleToRotMat(theta_vk_i(:,kStart)));
firstImuState.p_I_G = r_i_vk_i(:,kStart);

[msckfState, featureTracks, trackedFeatureIds] = initializeMSCKF(firstImuState, measurements{kStart}, camera, kStart);

%% ============================MAIN LOOP========================== %%

for state_k = kStart:kEnd
    %% ==========================STATE PROPAGATION======================== %%

    %Propagate state and covariance
    msckfState = propagateMsckfStateAndCovar(msckfState, measurements{state_k}, noiseParams);

    %Add camera pose to msckfState
    msckfState = augmentState(msckfState, camera, state_k+1);
    
    
    %% ==========================FEATURE TRACKING======================== %%
    % Add observations to the feature tracks, or initialize a new one
    % If an observation is -1, add the track to featureTracksToResidualize
    featureTracksToResidualize = {};
    for featureId = 1:20
        %IMPORTANT: state_k + 1 not state_k
        meas_k = measurements{state_k+1}.y(:, featureId);
        if ismember(featureId, trackedFeatureIds)
            if meas_k(1,1) == -1
                %Feature is not in view, remove from the tracked features
                [msckfState, camStates, camStateIndices] = removeTrackedFeature(msckfState, featureId);
                
                %Add the track, with all of its camStates, to the
                %residualized list
                track = featureTracks{trackedFeatureIds == featureId};
                if length(camStates) > MSCKFParams.minTrackLength
                    track.camStates = camStates;
                    track.camStateIndices = camStateIndices;
                    featureTracksToResidualize{end+1} = track;
                end
               
                %Remove the track
                featureTracks = featureTracks(trackedFeatureIds ~= featureId);
                trackedFeatureIds(trackedFeatureIds == featureId) = [];
                
            else
                %Append observation and append id to cam states
                featureTracks{trackedFeatureIds == featureId}.observations(:, end+1) = meas_k;
                
                %Add observation to current camera
                msckfState.camStates{end}.trackedFeatureIds(end+1) = featureId;
            end
        else
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
    %% ==========================FEATURE RESIDUAL CORRECTIONS======================== %%
    if ~isempty(featureTracksToResidualize)
        %H_o has more than 1 row, but it will be grown in our for loop like
        %a pet
        H_o = zeros(0, 12 + 6*length(msckfState.camStates));
        
        r_stacked = [];
        A = [];
        A_index = [1,1]; %Helper for building the A matrix

        for f_i = 1:length(featureTracksToResidualize)

            track = featureTracksToResidualize{f_i};

            %Estimate feature 3D location through Gauss Newton inverse depth
            %optimization
            
            camStatesGT = {};
            for c_i_temp = 1:length(track.camStates)
                camStatesGT{end+1} = groundTruthStates{track.camStates{c_i_temp}.state_k}.camState;
            end
            
            %[p_f_G] = calcGNPosEst(track.camStates, track.observations, noiseParams);
            p_f_G = groundTruthMap(:, track.featureId);
            
            %Calculate residual and Hoj 
            [r_j] = calcResidual(p_f_G, track.camStates, track.observations);
            [H_o_j, A_j] = calcHoj(p_f_G, msckfState, track.camStateIndices);

            % Stacked residuals and friends
            iStart = size(H_o,1)+1;
            iEnd = iStart+size(H_o_j, 1) - 1;

            H_o(iStart:iEnd, :) = H_o_j;
            r_stacked((end+1):(end+size(r_j,1)),1) = r_j;
            
            A(A_index(1):(A_index(1)+size(A_j,1)-1),A_index(2):(A_index(2)+size(A_j,2)-1)) = A_j;
            A_index = A_index + size(A_j);
        end

        [T_H, Q_1] = calcTH(H_o);
        r_n = Q_1'*A'*r_stacked;

        % Build MSCKF covariance matrix
        P = [msckfState.imuCovar, msckfState.imuCamCovar;
               msckfState.imuCamCovar', msckfState.camCovar];

        R_n = noiseParams.imageVariance*eye( size(Q_1, 2) );

        % Calculate Kalman gain
        K = P*T_H' / ( T_H*P*T_H' + R_n );

        % State correction
        deltaX = K * r_n;
        msckfState = updateState(msckfState, deltaX);

        % Covariance correction
        tempMat = (eye(12 + 6*size(msckfState.camStates,2)) - K*T_H);
        P_corrected = tempMat * P * tempMat' + K * R_n * K';

        msckfState.imuCovar = P_corrected(1:12,1:12);
        msckfState.camCovar = P_corrected(13:end,13:end);
        msckfState.imuCamCovar = P_corrected(1:12, 13:end);

        %% ==========================STATE HISTORY======================== %% 
        imuStates = updateStateHistory(imuStates, msckfState, camera, state_k);
        
        %% ==========================STATE PRUNING======================== %%
        %Remove any camera states with no tracked features
        msckfState = pruneStates(msckfState);
    end
end %for state_K = ...