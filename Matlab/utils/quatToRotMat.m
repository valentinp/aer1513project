function rotMat = quatToRotMat(quat)
%
% Converts a quaternion into a 3x3 rotation matrix using {1,i,j,k} convention
%
    if( size(quat,1) ~= 4 || size(quat,2) ~= 1 )
        error('Input quaternion must be 4x1');
    end
    
    if( abs(norm(quat) - 1) > eps )
        error('Input quaternion must be unit-length');
    end
    
    a = quat(1); % scalar component
    b = quat(2); % i component
    c = quat(3); % j component
    d = quat(4); % k component
    
    rotMat = [  a*a + b*b - c*c - d*d,  2*(b*c - a*d),          2*(b*d + a*c);
                2*(b*c + a*d),          a*a - b*b + c*c - d*d,  2*(c*d - a*b);
                2*(b*d - a*c),          2*(c*d + a*b),          a*a - b*b - c*c + d*d   ];
end