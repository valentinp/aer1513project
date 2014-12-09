function [rotMat] = rotMatFromPsi(psiVec)
%ROTMATFROMPSI Compute rotation matrix from minimal rotation parameters

psiMag = norm(psiVec);
rotMat = cos(psiMag)*eye(3) + (1 - cos(psiMag))*(psiVec/psiMag)*(psiVec/psiMag)' - sin(psiMag)*crossMat(psiVec/psiMag);

end

