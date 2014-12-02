function F = calcF(imuState_k, measurements_k)
    F = zeros(12,12);
    
    omegaHat = measurements_k.omega - imuState_k.b_g;
    C_IG = quatToRotMat(imuState_k.q_IG);
    
    F(1:3,1:3) = -crossMat(omegaHat);
    F(1:3,4:6) = -eye(3);
    F(4:6,10:12) = -C_IG;
  
end