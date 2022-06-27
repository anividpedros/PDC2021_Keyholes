clc
clear
close all
format shortG

addpath(genpath(pwd))

%% Script steps
% 1. Read heliocentric elements (spice file) - date close to the encounter
% 2. Modify for circular Earth and adjust NEO orbit for encounter still happening
% 3. Switch to geocentric Keplerian elements and find time of perigeetp
% 4. Coordinates of the modified target plane(ξp,ζp,vp)
% 5. Compute the corresponding target plane pre-encounter coordinates (input for Öpik).(ξ,η,ζ)
% 6. Use Öpik theory to compute post-encounter coordinates
% 7. Convert to Heliocentric: Use Öpik theory formulae

%% 1. Read heliocentric elements

cspice_furnsh( 'SPICEfiles/naif0012.tls.pc' )
cspice_furnsh( 'SPICEfiles/gm_de431.tpc' )
cspice_furnsh( 'SPICEfiles/pck00010.tpc' )
cspice_furnsh( 'SPICEfiles/de431_part-1.bsp' )
cspice_furnsh( 'SPICEfiles/de431_part-2.bsp' )

% 2021 PDC
cspice_furnsh( 'SPICEfiles/2021_PDC-s11-merged-DE431.bsp' )
ast_id= '-937014';
epoch = '2021 October 20 TDB';

et = cspice_str2et( epoch ) - 24*0.715*3600;

cons.AU  = cspice_convrt(1, 'AU', 'KM');
cons.GMs = cspice_bodvrd( 'SUN', 'GM', 10 );
cons.GMe = cspice_bodvrd( '399', 'GM', 1 );
cons.Re  = 6378.140;

cons.yr  = 365.25 * 24 * 3600 ;
cons.Day = 3600*24; 

state_eat = cspice_spkezr( '399', et, 'ECLIPJ2000', 'NONE', '10' );
kep_eat = cspice_oscelt( state_eat, et, cons.GMs );

state_ast = cspice_spkezr( ast_id, et, 'ECLIPJ2000', 'NONE', '10' );
kep_ast = cspice_oscelt( state_ast, et, cons.GMs );
sma_ast = kep_ast(1)/(1-kep_ast(2));

%% 2. Modify elements for encounter with simpler Earth model
% Make Earth Circular-Ecliptic -----
kep_eat(2:3) = 0;
cons.yr = 2*pi/sqrt(cons.GMs/kep_eat(1)^3);

% Modify NEO to try to make dCA smaller ----
kep_ast(3) = kep_ast(3)+.2;     % Increase inclination
kep_ast(5) = kep_ast(5)+.06;    % Adjust arg of perihelion to low MOID
kep_ast(6) = kep_ast(6)-0.0111; % Adjust timing for very close flyby

% Generate states moments before the flyby
dt = -24*24 * 3600;
t0 = et + dt;

state_ast = cspice_conics(kep_ast, t0);
state_eat = cspice_conics(kep_eat, t0);

%% 3. Switch to geocentric and find time periapses
% Planetocentric frame definition (Valsecchi 2003):
% Y: Direction of motion of the planet;
% X: Sun on negative -X;
% Z: Completes
state_ast_O = HillRot(state_eat, state_ast);
kep_ast_O   = cspice_oscelt( state_ast_O, t0, cons.GMe );

% Find periapsis time
MA_per = 0;

a_ast_O = kep_ast_O(1)/(1 - kep_ast_O(2));
n_ast_O = sqrt( -cons.GMe/a_ast_O(1)^3 );
dt_per  = (MA_per - kep_ast_O(6))/n_ast_O;


%% 4. Using ephemeris at date of perigee, obtain MTP coordinates
state_ast = cspice_conics(kep_ast, t0 + dt_per);
state_eat = cspice_conics(kep_eat, t0 + dt_per);
DU = norm( state_eat(1:3) );
TU = sqrt(DU^3/cons.GMs);
cons.AU = DU;

ap    = kep_eat(1)/(1-kep_eat(2))/DU ;
longp = mod( kep_eat(4)+kep_eat(5)+kep_eat(6), 2*pi ) ; % In general sense should be longitude

state_ast_per = HillRot(state_eat, state_ast);
r_ast_O = state_ast_per(1:3);
v_ast_O = state_ast_per(4:6);

% Coordinates at MTP - Perigee
b_p   = norm( r_ast_O );
V     = norm( v_ast_O );
theta_p = acos( v_ast_O(2)/V );
phi_p   = atan2( v_ast_O(1), v_ast_O(3) );
 
cp = cos(phi_p);   sp = sin(phi_p);
ct = cos(theta_p); st = sin(theta_p);
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];

xi_p   = r_ast_O(1)*cp - r_ast_O(3)*sp;
eta_p  = r_ast_O(1)*st*sp + r_ast_O(2)*ct + r_ast_O(3)*st*cp;
zeta_p = r_ast_O(1)*ct*sp - r_ast_O(2)*st + r_ast_O(3)*ct*cp; 


%% 5. From MTP to TP
% Asymptotic velocity: Energy conservation
GME = cons.GMe;
U = sqrt( V*V - 2*GME/b_p );

U_nd = U/(DU/TU);
m = cons.GMe/cons.GMs;

% Rescaling b-plane coordinates: Angular momentum conservation
b    = (V/U)*b_p  ;
xi   = (V/U)*xi_p ;
zeta = (V/U)*zeta_p ;
 
% Rotation of velocity vector
hgamma = atan( cons.GMe/(b*U*U) );
hv     = cross( r_ast_O, v_ast_O );
DCM    = PRV_2_DCM( -hgamma, hv/norm(hv) );
UVec   = (U/V)*DCM*v_ast_O;
theta  = acos(UVec(2)/U);
phi    = atan2(UVec(1),UVec(3));

cp = cos(phi);   sp = sin(phi);
ct = cos(theta); st = sin(theta);
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];



%% Grid
N = 200;
xi_val = linspace(-20,20,N)/cons.Re;
zeta_val = linspace(-20,20,N)/cons.Re;
orb_grid = NaN(6,N,N);
val_grid = NaN(2,N*N);
h = 2;
m = 3.003489614915766e-6;
count = 1;

for i = 1:N
    for j = 1:N
        
        xi0 = xi_val(i);
        zeta0 = zeta_val(i);
        [zeta2,theta1,phi1,xi1,zeta1,r_b_post,v_b_post] = opik_next(U_nd,theta,phi,xi0,zeta0,t0,h,m);
        kep_opik_post = opik_bplane_2_oe(theta1,phi1,zeta1,xi1,U_nd,phi,longp,ap)';
        
        % Adjust orbit elements
        kep_opik_post(1) = kep_opik_post(1)*(1-kep_opik_post(2))*DU ; % rp to sma
        kep_opik_post(6) = TA_2_MA(kep_opik_post(6),kep_opik_post(2));
        kep_opik_post(7:8) = [t0 + dt_per;cons.GMs]; % format for SPICE, perigee time added
        
        orb_grid(:,i,j) = kep_opik_post(1:6);
        val_grid(:,count) = [xi0; zeta0];
        count = count+1;
    
    end
end

figure
subplot(1,3,1)
contourf(xi_val,zeta_val,squeeze(orb_grid(3,:,:)))
colormap('jet')
colorbar
subplot(1,3,2)
contourf(xi_val,zeta_val,squeeze(orb_grid(2,:,:)))
colormap('jet')
colorbar
subplot(1,3,3)
contourf(xi_val,zeta_val,squeeze(orb_grid(1,:,:)))
colormap('jet')
colorbar
