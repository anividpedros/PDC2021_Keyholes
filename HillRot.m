function [Xo] = HillRot( Xp,Xa )

% Setup of the frame:
% Y as direction of velocity
% Z as angular momentum
% X completes the frame

Rp = Xp(1:3); rp = norm(Rp); 
Vp = Xp(4:6); vp = norm(Vp);
He = cross( Rp,Vp );
hp = norm(He);

Z  = He/hp;
Y  = Vp/vp;
X  = cross( Y, Z );
ON = [X Y Z]';

Xo = ON*(Xp(1:3) - Xa(1:3));

wp = (hp/rp/rp)*Z ;
vr = cross(wp,Xo) ;

Xo(4:6) = ON*(Xa(4:6) - Xp(4:6) - vr);