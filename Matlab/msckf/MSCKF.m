%% ============================Notation============================ %%
% X_sub_super
% q_ToFrom
% p_ofWhat_expressedInWhatFrame


%% =============================Setup============================== %%
clear;
close all;
clc;
addpath('utils');

% load('../datasets/dataset3.mat')
% load('../datasets/dataset3_fresh_10noisy.mat')
% load('../datasets/dataset3_fresh_10lessnoisy.mat')
% load('../datasets/dataset3_fresh_20lessnoisy.mat')
% load('../datasets/dataset3_fresh_40lessnoisy.mat')
load('../datasets/dataset3_fresh_60lessnoisy.mat')
% load('../datasets/dataset3_fresh_80lessnoisy.mat')
% load('../datasets/dataset3_fresh_100lessnoisy.mat')

%Dataset window bounds
% kStart = 500; kEnd = 1000;
kStart = 1215; kEnd = 1715;

%Set constant
numLandmarks = size(y_k_j,3);

%Set up the camera parameters
camera.c_u      = cu;                   % Principal point [u pixels]
camera.c_v      = cv;                   % Principal point [v pixels]
camera.f_u      = fu;                   % Focal length [u pixels]
camera.f_v      = fv;                   % Focal length [v pixels]
camera.b        = b;                    % Stereo baseline [m]
camera.q_CI     = rotMatToQuat(C_c_v);  % 4x1 IMU-to-Camera rotation quaternion
camera.p_C_I    = rho_v_c_v;            % 3x1 Camera position in IMU frame

%Set up the noise parameters
y_var = [1,1,1,1];
noiseParams.u_var_prime = 0.01;%y_var(1)/camera.f_u^2;
noiseParams.v_var_prime = 0.01;%y_var(2)/camera.f_v^2;

% [w, w bias change, v, v bias change]
noiseParams.Q_imu = 0.1*diag([w_var', 1e-4*ones(1,3), v_var', 1e-4*ones(1,3)]);
% [q, w bias, v bias, p]
noiseParams.initialIMUCovar = 1e-4 * eye(12);
    
%MSCKF parameters
msckfParams.minTrackLength = 10;     % Set to inf to dead-reckon only
msckfParams.maxTrackLength = 50;     % Set to inf to wait for features to go out of view
msckfParams.maxGNCost      = 10;     % Set to inf to allow any triangulation, no matter how bad
msckfParams.minRCOND       = 0;
msckfParams.doNullSpaceTrick = true;
msckfParams.doQRdecomp = true;


% IMU state for plotting etc. Structures indexed in a cell array
imuStates = cell(1,numel(t));
prunedStates = {};

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

% Important: Because we're idealizing our pixel measurements and the
% idealized measurements could legitimately be -1, replace our invalid
% measurement flag with NaN
y_k_j(y_k_j == -1) = NaN;

for state_k = kStart:kEnd 
    measurements{state_k}.dT    = dT(state_k);                      % sampling times
    measurements{state_k}.y     = squeeze(y_k_j(1:2,state_k,:));    % left camera only
    measurements{state_k}.omega = w_vk_vk_i(:,state_k);             % ang vel
    measurements{state_k}.v     = v_vk_vk_i(:,state_k);             % lin vel
    
    %Idealize measurements
    validMeas = ~isnan(measurements{state_k}.y(1,:));
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

[msckfState, featureTracks, trackedFeatureIds] = initializeMSCKF(firstImuState, measurements{kStart}, camera, kStart, noiseParams);
imuStates = updateStateHistory(imuStates, msckfState, camera, kStart);


%% ============================MAIN LOOP========================== %%

numFeatureTracksResidualized = 0;

for state_k = kStart:(kEnd-1)
    fprintf('state_k = %4d\n', state_k);
    
    
    %% ==========================STATE PROPAGATION======================== %%
    
    %Propagate state and covariance
    msckfState = propagateMsckfStateAndCovar(msckfState, measurements{state_k}, noiseParams);

    %Add camera pose to msckfState
    msckfState = augmentState(msckfState, camera, state_k+1);
    
    
    %% ==========================FEATURE TRACKING======================== %%
    % Add observations to the feature tracks, or initialize a new one
    % If an observation is -1, add the track to featureTracksToResidualize
    featureTracksToResidualize = {};
    
    for featureId = 1:numLandmarks
        %IMPORTANT: state_k + 1 not state_k
        meas_k = measurements{state_k+1}.y(:, featureId);
        
        outOfView = isnan(meas_k(1,1));
        
        if ismember(featureId, trackedFeatureIds)

            if ~outOfView
                %Append observation and append id to cam states
                featureTracks{trackedFeatureIds == featureId}.observations(:, end+1) = meas_k;
                
                %Add observation to current camera
                msckfState.camStates{end}.trackedFeatureIds(end+1) = featureId;
            end
            
            track = featureTracks{trackedFeatureIds == featureId};
            
            if outOfView ...
                    || size(track.observations, 2) >= msckfParams.maxTrackLength ...
                    || state_k+1 == kEnd
                                
                %Feature is not in view, remove from the tracked features
                [msckfState, camStates, camStateIndices] = removeTrackedFeature(msckfState, featureId);
                
                %Add the track, with all of its camStates, to the
                %residualized list
                if length(camStates) >= msckfParams.minTrackLength
                    track.camStates = camStates;
                    track.camStateIndices = camStateIndices;
                    featureTracksToResidualize{end+1} = track;
                end
               
                %Remove the track
                featureTracks = featureTracks(trackedFeatureIds ~= featureId);
                trackedFeatureIds(trackedFeatureIds == featureId) = [];                
            end
            
        elseif ~outOfView && state_k+1 < kEnd % && ~ismember(featureId, trackedFeatureIds)
            %Track new feature
            track.featureId = featureId;
            track.observations = meas_k;
            featureTracks{end+1} = track;
            trackedFeatureIds(end+1) = featureId;

            %Add observation to current camera
            msckfState.camStates{end}.trackedFeatureIds(end+1) = featureId;
        end
    end
     
    
    %% ==========================FEATURE RESIDUAL CORRECTIONS======================== %%
    if ~isempty(featureTracksToResidualize)
        H_o = [];
        r_o = [];
        R_o = [];

        for f_i = 1:length(featureTracksToResidualize)

            track = featureTracksToResidualize{f_i};

            %Estimate feature 3D location through Gauss Newton inverse depth
            %optimization
            [p_f_G, Jcost, RCOND] = calcGNPosEst(track.camStates, track.observations, noiseParams);
            % Uncomment to use ground truth map instead
             %p_f_G = groundTruthMap(:, track.featureId); Jcost = 0; RCOND = 1;
            
            if Jcost > msckfParams.maxGNCost || RCOND < msckfParams.minRCOND
                break;
            else
                numFeatureTracksResidualized = numFeatureTracksResidualized + 1;
                fprintf('Using new feature track with %d observations. Total track count = %d.\n',...
                    size(track.observations,2), numFeatureTracksResidualized);
            end
            
            %Calculate residual and Hoj 
            [r_j] = calcResidual(p_f_G, track.camStates, track.observations);
            R_j = diag(repmat([noiseParams.u_var_prime, noiseParams.v_var_prime], [1, numel(r_j)/2]));
            [H_o_j, A_j, H_x_j] = calcHoj(p_f_G, msckfState, track.camStateIndices);

            % Stacked residuals and friends
            if msckfParams.doNullSpaceTrick
                H_o = [H_o; H_o_j];

                if ~isempty(A_j)
                    r_o_j = A_j' * r_j;
                    r_o = [r_o ; r_o_j];

                    R_o_j = A_j' * R_j * A_j;
                    R_o(end+1 : end+size(R_o_j,1), end+1 : end+size(R_o_j,2)) = R_o_j;
                end
                
            else
                H_o = [H_o; H_x_j];
                r_o = [r_o; r_j];
                R_o(end+1 : end+size(R_j,1), end+1 : end+size(R_j,2)) = R_j;
            end
        end
        
        if ~isempty(r_o)
            % Put residuals into their final update-worthy form
            if msckfParams.doQRdecomp
                [T_H, Q_1] = calcTH(H_o);
                r_n = Q_1' * r_o;
                R_n = Q_1' * R_o * Q_1;
            else
                T_H = H_o;
                r_n = r_o;
                R_n = R_o;
            end           
            
            % Build MSCKF covariance matrix
            P = [msckfState.imuCovar, msckfState.imuCamCovar;
                   msckfState.imuCamCovar', msckfState.camCovar];

            % Calculate Kalman gain
            K = (P*T_H') / ( T_H*P*T_H' + R_n ); % == (P*T_H') * inv( T_H*P*T_H' + R_n )
%             K = (P*T_H') / ( T_H*P*T_H' + R_n + LMLambda*eye(size(R_n,1)));

            % State correction
            deltaX = K * r_n;
            msckfState = updateState(msckfState, deltaX);

            % Covariance correction
            tempMat = (eye(12 + 6*size(msckfState.camStates,2)) - K*T_H);
%             tempMat = (eye(12 + 6*size(msckfState.camStates,2)) - K*H_o);

            P_corrected = tempMat * P * tempMat' + K * R_n * K';

            msckfState.imuCovar = P_corrected(1:12,1:12);
            msckfState.camCovar = P_corrected(13:end,13:end);
            msckfState.imuCamCovar = P_corrected(1:12, 13:end);
           
%             figure(1); clf; imagesc(deltaX); axis equal; axis ij; colorbar;
%             drawnow;
            
        end
        
    end
    
        %% ==========================STATE HISTORY======================== %% 
        imuStates = updateStateHistory(imuStates, msckfState, camera, state_k+1);
        
        
        %% ==========================STATE PRUNING======================== %%
        %Remove any camera states with no tracked features
        [msckfState, deletedCamStates] = pruneStates(msckfState);
        
        if ~isempty(deletedCamStates)
            prunedStates(end+1:end+length(deletedCamStates)) = deletedCamStates;
        end    
        
        if max(max(msckfState.imuCovar(1:12,1:12))) > 1
            disp('omgbroken');
        end
        
%     figure(1); imagesc(msckfState.imuCovar(1:12,1:12)); axis equal; axis ij; colorbar;
%     drawnow;
end %for state_K = ...




%% ==========================PLOT ERRORS======================== %%
kNum = length(prunedStates);
p_C_G_est = NaN(3, kNum);
p_C_G_GT = NaN(3, kNum);
theta_CG_err = NaN(3,kNum);
err_sigma = NaN(6,kNum); % cam state is ordered as [rot, trans]

tPlot = NaN(1, kNum);

for k = 1:kNum
    state_k = prunedStates{k}.state_k;
    
    p_C_G_GT(:,k) = groundTruthStates{state_k}.camState.p_C_G;
    p_C_G_est(:,k) = prunedStates{k}.p_C_G;
    q_CG_est  = prunedStates{k}.q_CG;    
    
    theta_CG_err(:,k) = crossMatToVec( eye(3) ...
                    - quatToRotMat(q_CG_est) ...
                        * ( C_c_v * axisAngleToRotMat(theta_vk_i(:,kStart+k-1)) )' );
      
    err_sigma(:,k) = prunedStates{k}.sigma;
                    
    tPlot(k) = t(state_k);
end

rotLim = [-0.5 0.5];
transLim = [-0.5 0.5];

% Save estimates
msckf_trans_err = p_C_G_est - p_C_G_GT;
msckf_rot_err = theta_CG_err;
save('msckf_est.mat', 'msckf_trans_err', 'msckf_rot_err', 'tPlot', 'err_sigma');

mean(sqrt(sum(msckf_trans_err.^2, 1)/3))
mean(sqrt(sum(msckf_rot_err.^2, 1)/3))

% Translation Errors
figure
subplot(3,1,1)
plot(tPlot, p_C_G_est(1,:) - p_C_G_GT(1,:), 'LineWidth', 2)
hold on
plot(tPlot, 3*err_sigma(4,:), '--r')
plot(tPlot, -3*err_sigma(4,:), '--r')
ylim(transLim)
xlim([tPlot(1) tPlot(end)])
title('Translational Error')
ylabel('\delta r_x')


subplot(3,1,2)
plot(tPlot, p_C_G_est(2,:) - p_C_G_GT(2,:), 'LineWidth', 2)
hold on
plot(tPlot, 3*err_sigma(5,:), '--r')
plot(tPlot, -3*err_sigma(5,:), '--r')
ylim(transLim)
xlim([tPlot(1) tPlot(end)])
ylabel('\delta r_y')

subplot(3,1,3)
plot(tPlot, p_C_G_est(3,:) - p_C_G_GT(3,:), 'LineWidth', 2)
hold on
plot(tPlot, 3*err_sigma(6,:), '--r')
plot(tPlot, -3*err_sigma(6,:), '--r')
ylim(transLim)
xlim([tPlot(1) tPlot(end)])
ylabel('\delta r_z')
xlabel('t_k')

% Rotation Errors
figure
subplot(3,1,1)
plot(tPlot, theta_CG_err(1,:), 'LineWidth', 2)
hold on
plot(tPlot, 3*err_sigma(1,:), '--r')
plot(tPlot, -3*err_sigma(1,:), '--r')
ylim(rotLim)
xlim([tPlot(1) tPlot(end)])
title('Rotational Error')
ylabel('\delta \theta_x')


subplot(3,1,2)
plot(tPlot, theta_CG_err(2,:), 'LineWidth', 2)
hold on
plot(tPlot, 3*err_sigma(2,:), '--r')
plot(tPlot, -3*err_sigma(2,:), '--r')
ylim(rotLim)
xlim([tPlot(1) tPlot(end)])
ylabel('\delta \theta_y')

subplot(3,1,3)
plot(tPlot, theta_CG_err(3,:), 'LineWidth', 2)
hold on
plot(tPlot, 3*err_sigma(3,:), '--r')
plot(tPlot, -3*err_sigma(3,:), '--r')
ylim(rotLim)
xlim([tPlot(1) tPlot(end)])
ylabel('\delta \theta_z')
xlabel('t_k')