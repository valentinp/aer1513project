function [p_f_G] = calcGNPosEst(camStates, observations)
%CALCGNPOSEST Calculate the position estimate of the feature using Gauss
%Newton optimization
%   INPUT:
%   observations: 2xM matrix of pixel values
%   camStates: Cell array of M structs of camera poses
%   camera: intrinsic calibration
%   OUTPUT:
%   p_f_G: 3x1 feature vector in the global frame

%K is not needed if we assume observations are not pixels but x' = (u -
%c_u)/f_u
%K = [camera.f_u 0 camera.c_u; 0 camera.f_v camera.c_v; 0 0 1];

%Get initial estimate through intersection
%Use the first 2 camStates
C_12 = quatToRotMat(camStates{1}.q_CG)*quatToRotMat(camStates{2}.q_CG)';
t_21_1 = quatToRotMat(camStates{1}.q_CG)*(camStates{2}.p_C_G - camStates{1}.p_C_G);

p_f1_1_bar = triangulate(observations(:,1),observations(:,2),C_12, t_21_1);

xBar = p_f1_1_bar(1);
yBar = p_f1_1_bar(2);
zBar = p_f1_1_bar(3);

alphaBar = xBar/zBar;
betaBar = yBar/zBar;
rhoBar = 1/zBar;

xEst = [alphaBar; betaBar; rhoBar];

Cnum = length(camStates);

%Optimize
maxIter = 10;
for optI = 1:maxIter

    E = zeros(3, 2*Cnum);
    W = zeros(2*Cnum, 2*Cnum);
    errorVec = zeros(2*Cnum, 1);

    for iState = 1:Cnum
        %Form the weight matrix
        W((3*State - 2):(3*State),(3*State - 2):(3*State)) = diag([camera.n_u^2 camera.n_v^2]);

        C_i1 = quatToRotMat(camStates{iState}.q_CG)*quatToRotMat(camStates{1}.q_CG)';


        %Form the error vector
        zHat = observations(:, iState);
        h = C_i1*[alphaBar; betaBar; 1] + rhoBar*p_f1_1_bar;

        errorVec((2*State - 1):(2*State),1) = zHat - 1/h(3)*[h(1); h(2)];


        %Form the Jacobian
        dEdalpha = [-1/h(3)*C_i1(1,1) + (1/h(3)^2)*C_i1(3,1)*h(1); ...
                    -1/h(3)*C_i1(2,1) + (1/h(3)^2)*C_i1(3,1)*h(2)];

        dEdbeta = [-1/h(3)*C_i1(1,2) + (1/h(3)^2)*C_i1(3,2)*h(1); ...
                    -1/h(3)*C_i1(2,2) + (1/h(3)^2)*C_i1(3,2)*h(2)];

        dEdrho = [-1/h(3)*xBar + (1/h(3)^2)*zBar*h(1); ...
                    -1/h(3)*yBar + (1/h(3)^2)*zBar*h(2)];


        Eblock = [dEdalpha; dEdbeta; dEdrho];
        E((2*State - 1):(2*State), (3*State - 2):(3*State)) = Eblock;
    end
    
    %Calculate the cost function
    J = 0.5*errorVec'*(W\errorVec)
    
    %Solve!
    dx_star =  (E'*(W\E))\(-E'*(W\errorVec)); 
    
    xEst = xEst + dx_star;
    
    if norm(dx_star) < 0.01
        break;
    else
        alphaBar = xEst(1);
        betaBar = xEst(2);
        rhoBar = xEst(3);
        
        zBar = 1/rhoBar;
        yBar = betaBar*zBar;
        xBar = alphaBar*zBar;
    end
end

p_f_G = 1/xEst(3)*quatToRotMat(camStates{1}.q_CG)'*[xEst(1:2); 1] + camStates{1}.p_C_G; 


    function [p_f1_1] = triangulate(obs1, obs2, C_12, t_21_1)
        % triangulate Triangulates 3D points from two sets of feature vectors and a
        % a frame-to-frame transformation

           %Calculate unit vectors 
           v_1 = [obs1;1];
           v_2 = [obs2;1];
           v_1 = v_1/norm(v_1);
           v_2 = v_2/norm(v_2);

           A = [v_1 -C_12*v_2];
           b = t_21_1;

           scalar_consts = abs(A\b);
           p_f1_1 = (scalar_consts(1)*v_1 + scalar_consts(2)*C_12*v_2 + t_21_1)/2;
    end

end

