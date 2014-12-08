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
camera.c_u      = cu;
camera.c_v      = cv;
camera.f_u      = fu;
camera.f_v      = fv;
camera.b        = b;
camera.q_CI     = rotMatToQuat(C_c_v);
camera.p_C_I    = rho_v_c_v;

%Set up the noise parameters
noiseParams.z_1 = 1;
noiseParams.z_2 = 1;

% IMU state for plotting etc. Structures indexed in a cell array
imuState = cell{1,numel(t)};
imuState{1}.q_IG  = [zeros(3,1); 1];    %Global to IMU rotation quaternion
imuState{1}.p_I_G = zeros(3,1);         %IMU Position in the Global frame
imuState{1}.b_g   = zeros(3,1);         %Gyro bias
imuState{1}.b_v   = zeros(3,1);         %Velocity bias
imuState{1}.covar = zeros(12,12);       %IMU state covariance

% We don't really need these outside of msckfState, do we?
% camState = cell{1,numel(t)};
% camStates{1}.q_CG  = [zeros(3,1); 1];
% camStates{1}.p_C_G = zeros(3,1);

%msckfState.imuState
%msckfState.imuCovar
%msckfState.camCovar
%msckfState.imuCamCovar
%msckfState.camStates

% Measurements as structures all indexed in a cell array
dT = [0, diff(t)];
measurements = cell{1,numel(t)};
for k = kStart:kEnd 
    measurements{k}.dT    = dT(k);           % sampling times
    measurements{k}.y     = y_k_j(1:2,k,:);  % left camera only
    validMeas = (measurements{k}.y(1,:) ~= -1);
    measurements{k}.y(1,validMeas) = (measurements{k}.y(1,validMeas) - camera.c_u)/f_u; %Idealize measurements
    measurements{k}.y(2,validMeas) = (measurements{k}.y(2,validMeas) - camera.c_v)/f_v;
    measurements{k}.omega = w_vk_vk_i(:,k);  % ang vel
    measurements{k}.v     = v_vk_vk_i(:,k);  % lin vel
end

% Other constants
Nmax = 50;                      % max number of poses before triggering an update
imageVariance = mean(y_var);    % Slightly hacky. Used to compute the Kalman gain and corrected covariance in the EKF step

%Struct used to keep track of features
featureTracks = {};
trackedFeatureIds = [];
% featureTrack = {track1, track2, ...}
% track.featureId 
% track.observations
% track.k1
% track.k2


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
            track.k1 = state_k;
            track.k2 = state_k;
            featureTracks{end+1} = track;
            trackedFeatureIds(end+1) = featureId;
        end
    end
    
    %% ==========================FEATURE RESIDUAL CORRECTIONS======================== %%
    
    featuresToResidualize = []; %1xN matrix of feature ids (this is just the column of y_k_j) 
    featureTracksStartEnd = []; %2xN matrix of k1_j and k2_j for the jth feature track

    
    H_o = zeros( 2*length(featuresToResidualize) , 12 + size(msckfState.camStates,2) );
    r_stacked = [];
    A = [];

    for f_i = 1:length(featuresToResidualize)
        k1 = featureTracksStartEnd(1, f_i);
        k2 = featureTracksStartEnd(2, f_i);
        featureId = featuresToResidualize(f_i);

        %TODO: REPLACE THIS WITH featureTracksToResidualize which contains
        %observations
        observations = NaN(2, k2 - k1 + 1);

        for k = k1:k2
            observations(:, k-k1+1) =  measurements{k}.y(:, featureId);
        end

        %Estimate feature 3D location through Gauss Newton inverse depth
        %optimization
        [p_f_G] = calcGNPosEst(camStates, observations);
        
        %Calculate residual and Hoj 
        [r_j] = calcResidual(p_f_G, camStates, observations);
        [H_o_j, A_j] = calcHoj(p_f_G, msckfState, camStateIndex);

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