function quatProd = quatMult(quat1, quat2)
% 
% Multiplies two quaternions using {1,i,j,k} convention
%
    if( size(quat1,1) ~= 4 || size(quat1,2) ~= 1 ...
        || size(quat2,1) ~= 4 || size(quat2,2) ~= 1 )
        error('Input quaternions must be 4x1');
    end
    
    a = quat1(1); % scalar component
    b = quat1(2); % i component
    c = quat1(3); % j component
    d = quat1(4); % k component
    
    Qmat1 = [   a, -b, -c, -d;
                b,  a,  -d, c;
                c,  d,  a,  -b;
                d,  -c, b,  a   ];
                            
    quatProd = Qmat1 * quat2;
end