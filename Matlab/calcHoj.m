function [H_o_j, A] = calcHoj(p_f_G, msckfState, camStateIndex)
%CALCHOJ Calculates H_o_j according to Mourikis 2007
% Inputs: p_f_G: feature location in the Global frame
%         camState: camera pose in the Global frame
%         camStateIndex: i, with camState being the ith camera pose in the window       
% Outputs: H_o_j, A


camState = msckfState{camStateIndex + 1};
N = length(msckfState);

%The feature position in the camera frame
p_f_C = quatToRotMat(camState.q_CG)*(p_f_G - camState.p_C_G);

X = p_F_C(1);
Y = p_F_C(2);
Z = p_F_C(3);

J_i = (1/Z)*[1 0 -X/Z; 0 1 -Y/Z];


H_x_j = zeros(2, 12 + 6*N);
H_x_j(:,12+6*(camStateIndex-1) + 1:15+6*(camStateIndex-1) + 3) = J_i*crossMat(p_f_C);
H_x_j(:,(12+6*(camStateIndex-1) + 4):(15+6*(camStateIndex-1) + 6)) = -J_i*quatToRotMat(camState.q_CG);


A = null(H_x_j);
H_o_j = A'*H_x_j;
end

