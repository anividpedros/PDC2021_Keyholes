% Code that integrates Sun - Earth - Asteroid System
clc
clear
close all
format shortG

filename = 'main_sect3_multimethod_dist_circularEarth.m';
filepath = matlab.desktop.editor.getActiveFilename;
addpath(genpath(pwd))

cspice_furnsh( 'SPICEfiles/naif0012.tls.pc' )
cspice_furnsh( 'SPICEfiles/gm_de431.tpc' )
cspice_furnsh( 'SPICEfiles/pck00010.tpc' )
cspice_furnsh( 'SPICEfiles/de431_part-1.bsp' )
cspice_furnsh( 'SPICEfiles/de431_part-2.bsp' )

% 2021 PDC Should we still test the same case
cspice_furnsh( 'SPICEfiles/2021_PDC-s11-merged-DE431.bsp' )
ast_id= '-937014';
epoch = '2021 October 20 TDB';
et = cspice_str2et( epoch ) - 24*0.2*3600;

% Constants
cons.AU  = cspice_convrt(1, 'AU', 'KM');
cons.GMs = cspice_bodvrd( 'SUN', 'GM', 10 );
cons.GMe = cspice_bodvrd( '399', 'GM', 1 ); % gravitational parameter Earth
cons.Re  = 6378.140; % Earth Radius
cons.Rs = 696340;
cons.mu = cons.GMe/cons.GMs; % massa ratio
cons.yr  = 365.25 * 24 * 3600;
cons.Day = 3600*24; 
lstar = cons.AU;
rearth = cons.Re/lstar; % Adimensional Radius for Earth
rsun = cons.Rs/lstar;
tstar = sqrt(lstar^3/(cons.GMe+cons.GMs));
n = 1/tstar;

% initial state for asteroid
state_ast = cspice_spkezr(ast_id, et, 'ECLIPJ2000', 'NONE', 'SUN' );
state_earth = cspice_spkezr('399', et, 'ECLIPJ2000', 'NONE', 'SUN' );
barycenter = cons.GMe * state_earth(1:3) / (cons.GMe+cons.GMs);
w_vec = [0;0;1] * n;

state_earth(1:3) = state_earth(1:3) - barycenter;
theta = atan2(state_earth(2),state_earth(1));
rotmat = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];

state_ast(1:3) = state_ast(1:3) - barycenter;

state_ast(1:3) = rotmat * state_ast(1:3);
state_ast(4:6) = rotmat * state_ast(4:6);
state_ast(4:6) = state_ast(4:6) - cross(w_vec,state_ast(1:3));


state_ast(1:3) = state_ast(1:3)./lstar;
state_ast(4:6) = (state_ast(4:6)./lstar).*tstar;


t0 = 0; %? PENDING TO CHANGE
tf = 2; %? PENDING TO CHANGE, add event function based on distance, when minimum happens, distance starts increasing again maybe?
tv = [t0 tf];

% Propagation using CR3BP assumption
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',ones(1,6)*tol);
[t,X] = ode113(@(t,X) CR3BP(t,X,cons.mu),tv,state_ast,options);


%% Plotting Constants
[Xs,Ys,Zs] = sphere;
X2 = Xs * rearth;
Y2 = Ys * rearth;
Z2 = Zs * rearth;

X3 = Xs * rsun;
Y3 = Ys * rsun;
Z3 = Zs * rsun;

% get state post-encounter
figure
hold on
surf(X2+(1-cons.mu),Y2,Z2,'FaceColor',[0 0.569 1])
% surf(X3+(-cons.mu),Y3,Z3,'FaceColor',[1 0.5 0])
plot3(X(:,1),X(:,2),X(:,3))
axis equal