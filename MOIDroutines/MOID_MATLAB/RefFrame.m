function [incliB, omegaB, argpB] = RefFrame (A,B)
%% Reference frame rotation
% 3-1-3 Tranformation (orbital elements)
theta1 = A.Omega;
theta2 = A.i;
theta3 = A.argp;
RA = zeros (3,3); %from reference to A frame
RA(1,1) = cos(theta3)*cos(theta1)-sin(theta3)*sin(theta1)*cos(theta2);
RA(1,2) = sin(theta1)*cos(theta3)+cos(theta1)*cos(theta2)*sin(theta3);
RA(1,3) = sin(theta2)*sin(theta3);
RA(2,1) = -cos(theta1)*sin(theta3) - sin(theta1)*cos(theta2)*cos(theta3);
RA(2,2) = -sin(theta3)*sin(theta1) + cos(theta3)*cos(theta2)*cos(theta1);
RA(2,3) = sin(theta2)*cos(theta3);
RA(3,1) = sin(theta2)*sin(theta1);
RA(3,2) = -sin(theta2)*cos(theta1);
RA(3,3) = cos(theta2);

%% Preparing the orbit
x = [cos(B.Omega)*cos(B.argp) - sin(B.argp)*cos(B.i)*sin(B.Omega);
     sin(B.Omega)*cos(B.argp) + sin(B.argp)*cos(B.i)*cos(B.Omega);
     sin(B.i)*sin(B.argp)];
 
y = [-cos(B.Omega)*sin(B.argp) - cos(B.argp)*cos(B.i)*sin(B.Omega);
     -sin(B.Omega)*sin(B.argp) + cos(B.argp)*cos(B.i)*cos(B.Omega);
     sin(B.i)*cos(B.argp)];

z = [sin(B.i)*sin(B.Omega);
    -sin(B.i)*cos(B.Omega);
    cos(B.i)];

xn = RA * x;
yn = RA * y;
zn = RA * z;

incliB = atan2(sqrt(zn(1).^2 + zn(2).^2),zn(3));
omegaB = -atan2(zn(1),-zn(2));
argpB = -atan2(xn(3),yn(3));
end
