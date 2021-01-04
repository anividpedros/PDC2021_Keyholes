function M = rotationM(a,na)

% Angle verification
% if abs(a) > 2*pi    
%      error('Rotation Angle must be in RADIANS!!!')
% end

% Direction of rotation setup
if na < 0
    a = a * -1;
    na = na * -1;
end

% Rotation Matrices
if na==1
    M = [1 0 0;
        0 cos(a) sin(a);
        0 -sin(a) cos(a)];
elseif na==2
    M = [cos(a) 0 -sin(a)
        0 1 0;
        sin(a) 0 cos(a)];
elseif na==3
    M = [cos(a) sin(a) 0;
        -sin(a) cos(a) 0;
        0 0 1];
else
    error('Invalid Rotation Matrix Input')
end

end