addpath(genpath(pwd))
format longG

OE1 = [1 0.1 0.1 1 1 1];
OE2 = [1.5 0.3 0.5 .5 .5 .5];

A.sma   = OE1(1);
A.e     = OE1(2);
A.i     = OE1(3);
A.Omega = OE1(4);
A.argp  = OE1(5);

B.sma   = OE2(1);
B.e     = OE2(2);
B.i     = OE2(3);
B.Omega = OE2(4);
B.argp  = OE2(5);

MOID1 = MOID_SDG_win(OE1([1 2 4 3 5]), OE2([1 2 4 3 5]))

[MOID2, vA, vB] = ComputeMOID(A,B);

% % Visualization fxns from Oscar's library
% MOID_plotter(OE1,OE2,0,5);
% grid on
% plot( wrapToPi(vA)*180/pi,wrapToPi(vB)*180/pi, 'r+' )
% 
% %%
% ri = coordinates_AnAstro(0, [OE1(1:5) 0 vA], 0,0,1,1);
% rj = coordinates_AnAstro(0, [OE2(1:5) 0 vB], 0,0,1,1);
% D = norm( ri(1:3) - rj(1:3) );

[MOID3, vA, vB] = ComputeMOID_mex_win(A,B)

