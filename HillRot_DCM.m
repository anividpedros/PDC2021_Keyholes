function [NO, NO6] = HillRot_DCM(X_earth0)

% Setup of the frame:
% Y as direction of velocity
% Z as angular momentum
% X completes the frame

X = X_earth0(1:3); X = X/norm(X);
Y = X_earth0(4:6); Y = Y/norm(Y);
He = cross( X,Y );

Z  = He/norm(He);
X  = cross( Y, Z );

NO = [X Y Z];

NO6 = [NO zeros(3);
    zeros(3) NO];