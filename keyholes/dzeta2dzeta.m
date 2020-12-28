function dzeta2dzeta = dzeta2dzeta(U,theta,phi,xi,zeta,m,h,theta1,phi1,xi1,zeta1)

bSq = xi^2 + zeta^2;
cSq = m/U^2;
c = sqrt(cSq);
sth = sin(theta); cth = cos(theta);
sph = sin(phi);   cph = cos(phi);

%% Post-first-encounter quantities

U1 = U;
cth1 = cos(theta1); sth1 = sin(theta1);
sph1 = sin(phi1); cph1 = cos(phi1);
X01 = xi1/cph1;
Y01 = (xi1*cth1*(sph1/cph1) - zeta1)/sth1;

UX1 = U*sth1*sph1;
UY1 = U*cth1;
UZ1 = U*sth1*cph1;

Ztb1 = -xi1*sph1 + zeta1*cth1*cph1; % (see p. 1006 of AAS 15-355, with eta = 0)

%% Partial derivatives (post-first-encounter trigonometric functions)

% d(cos(theta'))/d(zeta)
bSqPcSq = bSq + cSq;
denomin = bSqPcSq^2;
numerat = 2*c*( bSqPcSq*sth + 2*zeta*( c*cth - zeta*sth ) );
dcth1dzeta = numerat/denomin;

% d(sin(phi'))/d(zeta)
bSqMcSq = bSq - cSq;
fourcSqxiSq = 4*cSq*xi^2;
aux1 = bSqMcSq*sth - 2*c*zeta*cth;
inside = ( aux1^2 + fourcSqxiSq );
denomin = sqrt(inside);
numerat = 2*sph*( zeta*sth - c*cth );
term1   = numerat/denomin;

aux2 = aux1*( zeta*sth - c*cth );
aux3 = aux1*sph - 2*c*xi*cph;
numerat = 2*aux2*aux3;
denomin = inside^1.5;
term2  = numerat/denomin;

dsph1dzeta = term1 - term2;

% d(sin(theta'))/d(zeta)
denomin = sqrt(inside);
numerat = bSqMcSq*cth + 2*c*zeta*sth;
dsth1dzeta = - numerat/denomin * dcth1dzeta;

% d(cos(phi'))/d(zeta)
denomin = sqrt(inside);
numerat = 2*cph* ( zeta*sth - c*cth );
term1 = numerat/denomin;

aux2 = aux1*( zeta*sth - c*cth );
aux3 = aux1*cph + 2*c*xi*sph;
numerat = 2*aux2*aux3;
denomin = inside^1.5;
term2 = numerat/denomin;

dcph1dzeta = term1 - term2;

%% Partial derivatives (post-first-encounter coordinates)

% d(X'(tb))/d(zeta)
fourcSqzeta = 4*cSq*zeta;
term1 = fourcSqzeta*( c*sth*sph + xi*cph );
term2 = ( bSq^2 - cSq^2 + fourcSqzeta*zeta )*cth*sph;
numerat = term1 + term2;
denomin = ( bSq + cSq )^2;
dX1tbdzeta = numerat/denomin;

% d(Y'(tb))/d(zeta)
numerat = fourcSqzeta*( c*cth - zeta*sth ) + ( cSq^2 - bSq^2 )*sth;
denomin = ( bSq + cSq )^2;
dY1tbdzeta = numerat/denomin;

% d(Z'(tb))/d(zeta)
term1 = fourcSqzeta*( c*sth*cph - xi*sph );
term2 = ( bSq^2 - cSq^2 + fourcSqzeta*zeta )*cth*cph;
numerat = term1 + term2;
denomin = ( bSq + cSq )^2;
dZ1tbdzeta = numerat/denomin;

% d(UX')/d(zeta)
dUX1dzeta = U*( dsth1dzeta*sph1 + sth1*dsph1dzeta );

% d(UY')/d(zeta)
dUY1dzeta = U*dcth1dzeta;

% d(UZ')/d(zeta)
dUZ1dzeta = U * ( dsth1dzeta*cph1 + sth1*dcph1dzeta );

% d(X0')/d(zeta)
dX01dzeta = dX1tbdzeta - (Ztb1/UZ1) * dUX1dzeta - ... 
     UX1 * ( dZ1tbdzeta/UZ1 - (Ztb1/UZ1^2) * ( dUZ1dzeta ) );
 
% d(Y0')/d(zeta)
dY01dzeta = dY1tbdzeta - (Ztb1/UZ1) * dUY1dzeta - ... 
     UY1 * ( dZ1tbdzeta/UZ1 - (Ztb1/UZ1^2) * ( dUZ1dzeta ) );

%% d(zeta')/d(zeta)

dzeta1dzeta = dX01dzeta * cth1 * sph1 + X01 * dcth1dzeta * sph1 + ...
    X01 * cth1 * dsph1dzeta + dY01dzeta * sth1 + Y01 * dsth1dzeta;

%% d(zeta'')/d(zeta)

% s(U', theta')
numerat = 2*pi * ( U1*cth1^2 + cth1*( 1 - U1^2 ) - 3*U1 );
denomin = ( 1 - U1^2 - 2*U1*cth1 )^2.5 * sth1;
sU1th1 = numerat/denomin;

dzeta2dzeta = h*sU1th1*dcth1dzeta + dzeta1dzeta;

end