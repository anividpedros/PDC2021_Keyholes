function Kep = Equi_2_Keplerian( A )

e = sqrt( A(2)*A(2) + A(3)*A(3) );
w = atan2( A(2), A(3) );

tani = sqrt( A(4)*A(4) + A(5)*A(5) );
th= atan2( A(4), A(5) );

% A(2) = e;
% A(3) = tani;
% A(4) = th;
% A(5) = w - th;

Kep = [A(1); e; tani; th; w-th; A(6) ];

% h = Eq(2);
% k = Eq(3);
% p = Eq(4);
% q = Eq(5);
% 
% e = sqrt( h*h + k*k );
% w = atan2( h, k );
% 
% tani = sqrt( p*p + q*q );
% i = tani;
% th= atan2( p, q );
% 
% Eq(2) = e;
% Eq(3) = i;
% Eq(4) = th;
% Eq(5) = w - th;
% 
% Kep = [Eq(1); e; i; th; w-th; Eq(6) ];