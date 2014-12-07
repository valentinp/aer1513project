function G = calcG()
% Multiplies the noise vector in the linearized continuous-time
% error state model

    G = zeros(12,12);
    
    G(1:3,1:3) = -eye(3);
    G(4:6,4:6) = eye(3);
    G(7:9,10:12) = eye(3);
    
end