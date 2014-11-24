function C = quatToRotMat(quat)
%
% Converts a quaternion into a 3x3 rotation matrix
% using the {i,j,k,1} convention
%
    if( size(quat,1) ~= 4 || size(quat,2) ~= 1 )
        error('Input quaternion must be 4x1');
    end
    
    if( abs(norm(quat) - 1) > eps )
        error('Input quaternion must be unit-length');
    end
    
    R = quatRightComp(quat)' * quatLeftComp(quat);
    C = R(1:3,1:3);
end