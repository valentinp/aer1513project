function [H_o_j, A_j] = calcHoj(p_f_G, msckfState, camStateIndices)
%CALCHOJ Calculates H_o_j according to Mourikis 2007
% Inputs: p_f_G: feature location in the Global frame
%         msckfState: the current window of states
%         camStateIndex: i, with camState being the ith camera pose in the window       
% Outputs: H_o_j, A

N = length(camStateIndices);
H_x_j = zeros(2, 12 + 6*N);
H_f_j = zeros(2, 3*N);

c_i = 1;
for camStateIndex = camStateIndices
    camState = msckfState.camStates{camStateIndex};

    %The feature position in the camera frame
    p_f_C = quatToRotMat(camState.q_CG)*(p_f_G - camState.p_C_G);

    X = p_f_C(1);
    Y = p_f_C(2);
    Z = p_f_C(3);

    J_i = (1/Z)*[1 0 -X/Z; 0 1 -Y/Z];

    H_f_j(:, (3*c_i - 2):3*c_i) = J_i*quatToRotMat(camState.q_CG);

    H_x_j(:,12+6*(camStateIndex-1) + 1:12+6*(camStateIndex-1) + 3) = J_i*crossMat(p_f_C);
    H_x_j(:,(12+6*(camStateIndex-1) + 4):(12+6*(camStateIndex-1) + 6)) = -J_i*quatToRotMat(camState.q_CG);

    c_i = c_i + 1;
end

A_j = null(H_f_j');
H_o_j = A_j'*H_x_j;

end

