function imuState_prop = propagateImuState(imuState_k, measurements_k)
% prop == propagated to k+1

    C_IG = quatToRotMat(imuState_k.q_IG);

    % Rotation state
    arg = measurements_k.omega - imuState_k.b_g;
    imuState_prop.q_IG = 0.5 * omegaMat(arg) * imuState_k.q_IG * measurements_k.dT;
    
    % Bias states
    imuState_prop.b_g = imuState_k.b_g;
    imuState_prop.b_v = imuState_k.b_v;
    
    % Translation state
    imuState_prop.p_I_G = C_IG' * (measurements_k.v - imuState_k.b_v) * measurements_k.dT;
    
end