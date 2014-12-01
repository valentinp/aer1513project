function [r_n_j] = calcResidual(r__j, A_j, Q_1)
%CALCRESIDUAL Calculates the left null space residual

r_n_j = Q_1'*A_j'*r__j;

end

