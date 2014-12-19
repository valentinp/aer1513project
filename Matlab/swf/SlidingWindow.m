% AER1513 A3
% Author: Valentin Peretroukhin
% December 2014
% Sliding Window Gauss Newton Optimization
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
addpath('utils')
load('../dataset3.mat')

%Set up appropriate structs
calibParams.c_u = cu;
calibParams.c_v = cv;
calibParams.f_u = fu;
calibParams.f_v = fv;
calibParams.b = b;

vehicleCamTransform.C_cv = C_c_v;
vehicleCamTransform.rho_cv_v = rho_v_c_v;

%Set up sliding window
kStart = 1215;
kEnd = 1515; 
kappa = 50; %Sliding window size
maxOptIter = 10;

k1 = kStart;
k2 = k1+kappa;
K = k2 - k1;  %There are K + 1 total states, since x0 is the k1th state

initialStateStruct = {};

%% Extract noise values
Q = diag([v_var; w_var]);
R = diag(y_var);

%% First create the initial guess using dead reackoning

%Use ground truth for the first state
firstState.C_vi = Cfrompsi(theta_vk_i(:,k1));
firstState.r_vi_i = r_i_vk_i(:,k1);

initialStateStruct{1} = firstState;
%There are K + 1 states (there is a '0' state)
for kIdx = 1:K
    k = kIdx + k1;
    imuMeasurement.omega = w_vk_vk_i(:, k-1);
    imuMeasurement.v = v_vk_vk_i(:, k-1);
    deltaT = t(k) - t(k-1);
    
    %Propagate the state forward
    initialStateStruct{kIdx+1} = propagateState(initialStateStruct{kIdx}, imuMeasurement, deltaT);
end

%%

%Slide the window along
stateVecHistStruct = {};
stateSigmaHistMat = [];

%Keep track of triangulated landmarks
rho_i_pj_i_est = -1*ones(3, 20);


for k1 = kStart:kEnd    
tic    
k2 = k1+kappa;

%How many exteroceptive measurements do we have?
%NOTE: k1 is the 0th state
totalLandmarkObs = 0;
lmObsVec = zeros(1, K);

for k = (k1+1):k2
    lmObsVec(k-k1) = sum(y_k_j(1, k, :) > -1);
    totalLandmarkObs = totalLandmarkObs + sum(y_k_j(1, k, :) > -1);
end
totalLandmarkObs

%To initialize G-N, we propagate all the states for the first window
%and then only propagate the most recent state, re-using the rest
if k1 == kStart
    currentStateStruct = initialStateStruct;
else
        currentStateStruct = currentStateStruct(2:end);

        %Extract the measurement
        imuMeasurement.omega = w_vk_vk_i(:, k2-1);
        imuMeasurement.v = v_vk_vk_i(:, k2-1);
        deltaT = t(k2) - t(k2-1);

        %Propagate the state forward
        currentStateStruct{end+1} = propagateState(currentStateStruct{end}, imuMeasurement, deltaT);
end

%Define the optimal state
optimalStateStruct = currentStateStruct;
Jbest = Inf;
dx = Inf;

for optIdx = 1:maxOptIter

%Error Vector
errorVector = NaN(6*K+4*totalLandmarkObs, 1);
%This helper index will keep track of where we need to insert our next
%errors
errorVectorHelperIdx = 1;

%H and T
H = sparse(6*K+4*totalLandmarkObs, 6*(K+1));
T = sparse(6*K+4*totalLandmarkObs, 6*K+4*totalLandmarkObs);

%Helper indices that keep track of the row number of the last block entry
%into H and T
HHelperIdx = 1;
THelperIdx = 1;


for kIdx = 1:K
    k = kIdx + k1;
    imuMeasurement.omega = w_vk_vk_i(:, k-1);
    imuMeasurement.v = v_vk_vk_i(:, k-1);
    deltaT = t(k) - t(k-1);
    
    %==== Build the interoceptive error and Jacobians=====%
    %Note that there are K+1 states (the 0th state is the 1st element)
    kState = currentStateStruct{kIdx+1};
    kMinus1State = currentStateStruct{kIdx};
    intErrorVec = imuError(kState, kMinus1State, imuMeasurement, deltaT);
    H_x_k =  H_xfn(kMinus1State, imuMeasurement, deltaT );
    H_w_k = H_wfn(kMinus1State);
    
    
    %==== Build the exteroceptive error and Jacobians=====%
    validLmObsId = find(y_k_j(1, k, :) > -1);
    if ~isempty(validLmObsId)
        
        extErrorVec = NaN(4*length(validLmObsId), 1);
        G_x_k = NaN(4*length(validLmObsId),6);
        
        idx = 1;
        for lmId = validLmObsId'
            
            yMeas = y_k_j(:, k, lmId);
            %Find the ground truth position of the observed landmark
            %rho_pi_i_check = rho_i_pj_i(:, lmId);
            
            %if (rho_i_pj_i_est(1, lmId) == -1)
                 %Use triangulation to find the position of the landmark
                rho_pc_c = triangulate(yMeas, calibParams);
                rho_pi_i = kState.C_vi'*(vehicleCamTransform.C_cv'*rho_pc_c + vehicleCamTransform.rho_cv_v) +  kState.r_vi_i;
             %   rho_i_pj_i_est(:, lmId) = rho_pi_i;
            %else
                %rho_pi_i = rho_i_pj_i_est(:, lmId);
            %end
           
            extErrorVec(idx:idx+3, 1) = stereoCamError(yMeas, kState, vehicleCamTransform, rho_pi_i, calibParams);
            G_x_k(idx:idx+3, :) = G_xfn(kState, vehicleCamTransform, rho_pi_i, calibParams );
            idx = idx + 4;
        end
    else
        extErrorVec = [];
    end
    
    %Update matrices 
    %==== Error vector =====
    combinedErrorVec = [intErrorVec; extErrorVec];
    errorVector(errorVectorHelperIdx:(errorVectorHelperIdx + length(combinedErrorVec) - 1) ,1) = combinedErrorVec;
    errorVectorHelperIdx = errorVectorHelperIdx + length(combinedErrorVec);
    
    %==== H matrix =====
    Hblock = zeros(6+4*length(validLmObsId), 12);
    Hblock(1:6,1:6) = -H_x_k;
    Hblock(1:6,7:12) = eye(6);
    if ~isempty(validLmObsId)
        Hblock(7:(7+4*length(validLmObsId) - 1), 7:12) = -G_x_k;
    end
    Hblockrows = size(Hblock, 1);
    
    H(HHelperIdx:(HHelperIdx + Hblockrows - 1), 1+6*(kIdx-1):12+6*(kIdx-1) ) = Hblock;
    HHelperIdx = HHelperIdx + Hblockrows;
    
    %==== T matrix =====
    T_k = zeros(6+4*length(validLmObsId), 6+4*length(validLmObsId));
    T_k(1:6, 1:6) = H_w_k*Q*deltaT^2*H_w_k';
    
    %Here, G_n_k is identity, so we can just repeat the variances along the
    %diagonal
    obsVar = diag(R);
    T_k(7:end, 7:end) = diag(repmat(obsVar, [length(validLmObsId),1]));
    Tksize = size(T_k, 1);
    T(THelperIdx:(THelperIdx + Tksize - 1), THelperIdx:(THelperIdx + Tksize - 1)) = T_k;
    THelperIdx = THelperIdx + Tksize;
end

    %Calculate scalar objective
    Jnew = 0.5*errorVector'*(T\errorVector);
    
    
    if Jnew < Jbest 
        optimalStateStruct = currentStateStruct;
        Jbest = Jnew;
    else 
        break;
    end
    
    %Check for convergence
    if norm(dx) < 1e-2
        disp('Converged!')
        break;
    end

    % Solve for the optimal step size!
    dx = (H'*(T\H))\(-H'*(T\errorVector));
    currentStateStruct = updateStateStruct(currentStateStruct, dx);

   
end %End optimization iterations

currentStateStruct = optimalStateStruct;


if optIdx == maxOptIter
    fprintf('Warning: Failed to converge! \n');
end
    fprintf('%d done. J = %.5f. %d iterations. \n', k1, Jbest, optIdx)

%Extract variance of states
stateCov = inv(H'*(T\H));
stateVar = diag(stateCov);

%Keep track of the first state in the window
if ~all(stateVar > 0)
    warning('Variances not positive');
end
stateVecHistStruct{k1 - kStart + 1} = currentStateStruct{1};
stateSigmaHistMat(:,k1 - kStart + 1) = sqrt(abs(stateVar(1:6)));
toc
end %End Sliding window

sigma_x = (stateSigmaHistMat(1,:));
sigma_y = (stateSigmaHistMat(2,:));
sigma_z = (stateSigmaHistMat(3,:));
sigma_th1 = (stateSigmaHistMat(4,:));
sigma_th2 = (stateSigmaHistMat(5,:));
sigma_th3 = (stateSigmaHistMat(6,:));


%% Plot error and variances
addpath('/Users/valentinp/Research/MATLAB/export_fig'); %Use Oliver Woodford's awesome export_fig package to get trimmed PDFs

rotErrVec = zeros(3, length(stateVecHistStruct));
transErrVec = zeros(3, length(stateVecHistStruct));


for stIdx = 1:length(stateVecHistStruct)
    k = kStart+stIdx-1;
    transErrVec(:,stIdx) = stateVecHistStruct{stIdx}.r_vi_i - r_i_vk_i(:,k);
    eRotMat = eye(3) - stateVecHistStruct{stIdx}.C_vi*Cfrompsi(theta_vk_i(:,k))';
    rotErrVec(:, stIdx) = [eRotMat(3,2); eRotMat(1,3); eRotMat(2,1)];
end

transLim = 0.2;
rotLim = 0.2;
recycleStates = 'Yes';

figure
subplot(3,1,1)
plot(t(kStart:kEnd), transErrVec(1,:), 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), 3*sigma_x, '--r')
plot(t(kStart:kEnd), -3*sigma_x, '--r')
ylim([-transLim transLim])
title(sprintf('Translational Error (\\kappa = %d, Recycle States? %s)', kappa, recycleStates))
ylabel('\delta r_x')


subplot(3,1,2)
plot(t(kStart:kEnd), transErrVec(2,:), 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), 3*sigma_y, '--r')
plot(t(kStart:kEnd), -3*sigma_y, '--r')
ylim([-transLim transLim])
ylabel('\delta r_y')

subplot(3,1,3)
plot(t(kStart:kEnd), transErrVec(3,:), 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), 3*sigma_z, '--r')
plot(t(kStart:kEnd), -3*sigma_z, '--r')
ylim([-transLim transLim])
ylabel('\delta r_z')
xlabel('t_k')
%set(gca,'FontSize',12)
%set(findall(gcf,'type','text'),'FontSize',12)

figure
subplot(3,1,1)
plot(t(kStart:kEnd), rotErrVec(1,:), 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), 3*sigma_th1, '--r')
plot(t(kStart:kEnd), -3*sigma_th1, '--r')
ylim([-rotLim rotLim])
title(sprintf('Rotational Error (\\kappa = %d, Recycle States? %s)', kappa, recycleStates))
ylabel('\delta\theta_x')

 
subplot(3,1,2)
plot(t(kStart:kEnd), rotErrVec(2,:), 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), 3*sigma_th2, '--r')
plot(t(kStart:kEnd), -3*sigma_th2, '--r')
ylim([-rotLim rotLim])
ylabel('\delta\theta_y')

subplot(3,1,3)
plot(t(kStart:kEnd), rotErrVec(3,:), 'LineWidth', 1.2)
hold on
plot(t(kStart:kEnd), 3*sigma_th3, '--r')
plot(t(kStart:kEnd), -3*sigma_th3, '--r')
ylim([-rotLim rotLim])
ylabel('\delta\theta_z')
xlabel('t_k')
