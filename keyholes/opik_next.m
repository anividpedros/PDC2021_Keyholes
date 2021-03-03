function [zeta2,theta1,phi1,xi1,zeta1,b1_XYZ,U1_XYZ] = opik_next(U,theta,phi,xi,zeta,t0,h,m)
% Propagates Öpik state vector from one encounter to the one taking place
% "h" asteroid periods afterwards.

% Earth mass (solar masses)
% m = 3.0034896149157658e-06;

%% 1. Propagate from pre- to post-encounter quantities

% gamma
c = m/U^2;
bSq = xi^2 + zeta^2;
cSq = c^2; b = sqrt(bSq);
gamma = atan2(2*b*c,bSq - cSq);

% Trig functions
cga = cos(gamma); sga = sin(gamma);
cth = cos(theta); sth = sin(theta);
cphi = cos(phi);   sphi = sin(phi);
tphi = sphi/cphi;
cpsi = zeta/b;   spsi = xi/b;

% Pre-encounter velocity vector, impact parameter, and angular momentum
% in Cartesian coordinates.
U_XYZ = U*[sth*sphi; cth; sth*cphi];
R_XYZ2b = [cphi, 0, -sphi; sth*sphi, cth, cphi*sth; cth*sphi, -sth, cth*cphi]; % See AAS 15-355
b_XYZ = R_XYZ2b'*[xi;0;zeta];
hn_XYZ = cross(b_XYZ,U_XYZ);
hn_XYZ = hn_XYZ/norm(hn_XYZ);

% Pre-encounter SMA
a = 1/(1 - U^2 - 2*U*cth);

% Post-encounter velocity
U1_XYZ = ax_angle(U_XYZ,gamma,hn_XYZ);
U1 = U;
cth1 = U1_XYZ(2)/U1; theta1 = acos(cth1);
phi1 = atan2(U1_XYZ(1),U1_XYZ(3));
sphi1 = sin(phi1);
cphi1 = cos(phi1);

% Post-encounter b-plane coordinates
term1 = ((bSq - cSq)*sth - 2*c*zeta*cth)^2;
term2 = 4*cSq*xi^2;
den = sqrt(term1 + term2);
xi1 = ((bSq + cSq)*xi*sth)/den;
zeta1 = ((bSq - cSq)*zeta*sth - 2*bSq*c*cth)/den;

% Post-encounter semi-major axis
term1 = (bSq + cSq)*(1 - U^2);
term2 = -2*U*( (bSq - cSq)*cth + 2*c*zeta*sth );
den = term1 + term2;
a1 = (bSq + cSq)/den;

%% 2. From post-encounter to next pre-encounter

U2 = U1; theta2 = theta1; phi2 = phi1;

% Assume MOID unchanged (Keplerian propagation)
xi2 = xi1;
term1 = mod(2*pi*h*(a1)^1.5 + pi, 2*pi);
term2 = (term1 - pi)*sin(theta1);
zeta2 = zeta1 - term2;

%% Additional Outputs
b1_XYZ = ax_angle(b_XYZ,gamma,hn_XYZ);

end

