function [Xo] = HillRotInv( Xp,Xa )

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
NO = [X Y Z]';

dx = NO*Xa(1:3);
Xo = Xp(1:3) + dx;

wp = (hp/rp/rp)*Z ;
vr = cross(wp,dx) ;

Xo(4:6) = NO*Xa(4:6) + Xp(4:6) + vr;