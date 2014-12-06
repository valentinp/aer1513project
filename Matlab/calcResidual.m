function [r_j] = calcResidual(p_f_G, camStates, measurements)
%CALCRESIDUAL Calculates the residual for a feature position

% measurements is 2 x M_j
% camStates is a cell array of the camState structs for the states
%   included in measurements

    for i = 1:size(camStates,2)
        C_CG = quatToRotMat(camStates{i}.q_CG);
        zhat_i_j = C_CG * (p_f_G - camStates{i}.p_C_G);
        
        iStart = 2*(i-1)+1;
        iEnd = 2*i;
        r_j = measurements(iStart:iEnd) - zhat_i_j;
    end
        
end