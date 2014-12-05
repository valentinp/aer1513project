function msckfState_aug = augmentState(msckfState, camera)
    
    C_IG = quatToRotMat(msckfState.imuState.q_IG);
    q_CG = quatLeftComp(camera.q_CI) * msckfState.imuState.q_IG;
    p_C_G = msckfState.imuState.p_I_G + C_IG' * camera.p_C_I;

    P = [msckfState.imuCovar, msckfState.imuCamCovar;
        msckfState.imuCamCovar', msckfState.camCovar];
    
    J = calcJ(camera, msckfState.imuState, msckfState.camStates);
    
    N = size(msckfState.cameraStates,2);
    
    tempMat = [eye(15+6*N); J];
    
    P_aug = tempMat * P * tempMat';
    
    msckfState_aug = msckfState;
    msckfState_aug.camStates{N+1}.p_C_G = p_C_G;
    msckfState_aug.camStates{N+1}.q_CG = q_CG;
    msckfState_aug.imuCovar = P_aug(1:15,1:15);
    msckfState_aug.camCovar = P_aug(16:end,16:end);
    msckfState_aug.imuCamCovar = P_aug(1:15, 16:end);
end