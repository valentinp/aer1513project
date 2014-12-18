function [G_x] = G_xfn(kState, vehicleCamTransform, rho_pi_i, calibParams )
%G_XFN Calculates 4x6 G_x matrix for the stereo camera model

C_cv = vehicleCamTransform.C_cv;
rho_cv_v = vehicleCamTransform.rho_cv_v;
C_vi = kState.C_vi;
r_vi_i = kState.r_vi_i;

%Calculate the nominal feature vector
p_pc_c_bar = C_cv*(C_vi*(rho_pi_i - r_vi_i) - rho_cv_v); 


%Caculate the Jacobian of the transformation (L) and then pass it through
%the model g

L = [-C_cv*C_vi, C_cv*crossMat(C_vi*(rho_pi_i - r_vi_i))];
G_x = gJacob(p_pc_c_bar, calibParams)*L;

end

