%% ============================Notation============================ %%
% X_sub_super
% q_FromTo
% p_ofWhat_expressedInWhatFrame


%% =============================Setup============================== %%
addpath('utils');
load('dataset3.mat')

%Dataset window bounds
kStart = 1215;
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
noiseParams.imageVariance = mean([noiseParams.z_1, noiseParams.z_2]);  % Slightly hacky. Used to compute the Kalman gain and corrected covariance in the EKF step

% IMU state for plotting etc. Structures indexed in a cell array
imuState = cell(1,numel(t));
% imuState{k}.q_IG         4x1 Global to IMU rotation quaternion
% imuState{k}.p_I_G        3x1 IMU Position in the Global frame
% imuState{k}.b_g          3x1 Gyro bias
% imuState{k}.b_v          3x1 Velocity bias
% imuState{k}.covar        12x12 IMU state covariance

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

for state_k = kStart:kEnd 
    measurements{state_k}.dT    = dT(state_k);                      % sampling times
    measurements{state_k}.y     = squeeze(y_k_j(1:2,state_k,:));    % left camera only
    measurements{state_k}.omega = w_vk_vk_i(:,state_k);             % ang vel
    measurements{state_k}.v     = v_vk_vk_i(:,state_k);             % lin vel
    
    %Idealize measurements
    validMeas = (measurements{state_k}.y(1,:) ~= -1);
    measurements{state_k}.y(1,validMeas) = (measurements{state_k}.y(1,validMeas) - camera.c_u)/camera.f_u;
    measurements{state_k}.y(2,validMeas) = (measurements{state_k}.y(2,validMeas) - camera.c_v)/camera.f_v;
end

% Other constants
Nmax = 50;                      % max number of poses before triggering an update

%Struct used to keep track of features
featureTracks = {};
trackedFeatureIds = [];

% featureTrack = {track1, track2, ...}
% track.featureId 
% track.observations



%% ==========================Initial State======================== %%
%Use ground truth for first state and initialize feature tracks with
%feature observations



%% ============================MAIN LOOP========================== %%

for state_k = kStart:kEnd
    %% ==========================STATE PROPAGATION======================== %%

    %Propagate state and covariance
    msckfState = propagateMsckfCovar(msckfState, measurements{state_k}, noiseParams);

    %Add camera pose to msckfState
    msckfState = augmentState(msckfState, camera);
    
    
    %% ==========================FEATURE TRACKING======================== %%
    % Add observations to the feature tracks, or initialize a new one
    % If an observation is -1, add the track to featureTracksToResidualize
    featureTracksToResidualize = {};
    for featureId = 1:20
        meas_k = measurements{state_k}.y(:, featureId);
        if ismember(featureId, trackedFeatureIds)
            if meas_k(1,1) == -1
                %Add to residualize queue and remove from the original
                %struct
                featureTracksToResidualize{end+1} = featureTracks{trackedFeatureIds == featureId};
                featureTracks = featureTracks(trackedFeatureIds ~= featureId);
                trackedFeatureIds(trackedFeatureIds == featureId) = [];
            else
                %Append observation and increase k2
                featureTracks{trackedFeatureIds == featureId}.observations(:, end+1) = meas_k;
                featureTracks{trackedFeatureIds == featureId}.k2 = featureTracks{trackedFeatureIds == featureId}.k2 + 1;
            end
        else
            %Track new feature
            track.featureId = featureId;
            track.observations = meas_k;
            track.cameraIndices = [length(msckfState.camStates)];
            featureTracks{end+1} = track;
            trackedFeatureIds(end+1) = featureId;
        end
    end
    
    %% ==========================FEATURE RESIDUAL CORRECTIONS======================== %%
    H_o = zeros( 2*length(featuresToResidualize) , 12 + size(msckfState.camStates,2) );
    r_stacked = [];
    A = [];

    for f_i = 1:length(featureTracksToResidualize)
        
        track = featureTracksToResidualize{f_i};

        %Estimate feature 3D location through Gauss Newton inverse depth
        %optimization
        [p_f_G] = calcGNPosEst(camStates, track.observations);
        
        %Calculate residual and Hoj 
        [r_j] = calcResidual(p_f_G, camStates, observations);
        [H_o_j, A_j] = calcHoj(p_f_G, msckfState, camStateIndices);

        % Stacked residuals and friends
        iStart = 2*(f_i-1)+1;
        iEnd = iStart+2;

        H_o(iStart:iEnd, :) = H_o_j;
        r_stacked(end+1,:) = r_j;
        A(:, end+1) = A_j;
    end

    [T_H, Q_1] = calcTH(H_o);
    r_n_j = Q_1'*A'*r_stacked;

    % Build MSCKF covariance matrix
    P = [msckfState.imuCovar, msckfState.imuCamCovar;
           msckfState.imuCamCovar', msckfState.camCovar];

    R_n = imageVariance*eye( size(Q_1, 2) );

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

    %% ==========================STATE PRUNING======================== %%
    %Remove any camera states with no tracked features
    
end %for state_K = ...