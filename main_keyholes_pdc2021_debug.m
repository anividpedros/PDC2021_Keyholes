% GENERATION OF MODIFIED KEYHOLES ON THE B-PLANE
% 7th IAA Planetary Defense Conference 2021
clc
clear
close all
format shortG

filename = 'main_keyholes_pdc2021_debug.m';
filepath = matlab.desktop.editor.getActiveFilename;
currPath = filepath(1:(end-length(filename)));
cd(currPath)
addpath(genpath(pwd))

%% Script steps
% 1. Read heliocentric elements (spice file) - date close to the encounter
% 2. Modify for circular Earth and adjust NEO orbit for encounter still happening
% 3. Switch to geocentric Keplerian elements and find time of perigeetp
% 4. Coordinates of the modified target plane(ξp,ζp,vp)
% 5. Compute the corresponding target plane pre-encounter coordinates (input for Öpik).(ξ,η,ζ)
% 6. Use Öpik theory to compute post-encounter coordinates
% 7. Convert to Heliocentric: Use Öpik theory formulae

%% Code To-Do
% - Substitute cspice_conics calls
% - Delta t's of the encounter
% - !!! If Earth not circular, review definition of Hill Frame

%% 1. Read heliocentric elements
mice_local_path = 'C:\Users\Oscar\odrive\Google Drive (2)\MSc-ASEN\Research\Code2020\EncountersCode\mice';
% mice_local_de431 = '2nd Author write here your directory';

dir_local_de431 = 'C:\Users\Oscar\Documents\Spice-Kernels\';
% dir_local_de431 = '2nd Author write here your directory';

addpath(genpath(mice_local_path))
cspice_furnsh( 'naif0012.tls.pc' )
cspice_furnsh( 'gm_de431.tpc' )
cspice_furnsh( 'pck00010.tpc' )
cspice_furnsh( [dir_local_de431 'de431_part-1.bsp'] )
cspice_furnsh( [dir_local_de431 'de431_part-2.bsp'] )

% 2021 PDC
cspice_furnsh( '2021_PDC-s11-merged-DE431.bsp' )
ast_id= '-937014';
epoch = '2021 October 20 TDB';

et = cspice_str2et( epoch ) + 24*.715*3600;

cons.AU  = cspice_convrt(1, 'AU', 'KM');
cons.GMs = cspice_bodvrd( 'SUN', 'GM', 10 );
% cons.GMe = cspice_bodvrd( '3', 'GM', 1 );
cons.GMe = 398600.43543609593;
cons.Re  = 6378.140;

cons.yr  = 365.25 * 24 * 3600 ;
cons.Day = 3600*24; 

state_eat = cspice_spkezr( '399', et, 'ECLIPJ2000', 'NONE', '10' );
kep_eat = cspice_oscelt( state_eat, et, cons.GMs );

state_ast = cspice_spkezr( ast_id, et, 'ECLIPJ2000', 'NONE', '10' );
kep_ast = cspice_oscelt( state_ast, et, cons.GMs );


%% 2. Modify elements for encounter with simpler Earth model
% Make Earth Circular-Ecliptic -----
kep_eat(2:3) = 0;

% Modify NEO to try to make dCA smaller ----
kep_ast(3) = kep_ast(3)+.2;     % Increase inclination
kep_ast(5) = kep_ast(5)+.06;    % Adjust arg of perihelion to low MOID
kep_ast(6) = kep_ast(6)-0.0111; % Adjust timing for very close flyby

% Generate states moments before the flyby
dt = -4*24 * 3600;
t0 = et + dt;

state_ast = cspice_conics(kep_ast, t0);
state_eat = cspice_conics(kep_eat, t0);
DU = norm( state_eat(1:3) );
TU = sqrt(DU^3/cons.GMs);


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


%% 4. MTP coordinates
% Kind of explained in Farnocchia 2019 (being closest approach).
% Cartesian to Opik coordinates: generic, as in Valsecchi2003/Valsecchi2015
state_ast_per = cspice_conics( kep_ast_O, t0 + dt_per );
r_ast_O = state_ast_per(1:3);
v_ast_O = state_ast_per(4:6);

% Coordinates at MTP - Perigee
b_p   = norm( r_ast_O )
V     = norm( v_ast_O );
theta = acos( v_ast_O(2)/V );
phi   = atan2( v_ast_O(1), v_ast_O(3) );
 
cp = cos(phi);   sp = sin(phi);
ct = cos(theta); st = sin(theta);
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];

xi_p   = r_ast_O(1)*cp - r_ast_O(3)*sp;
eta_p  = r_ast_O(1)*st*sp + r_ast_O(2)*ct + r_ast_O(3)*st*cp;
zeta_p = r_ast_O(1)*ct*sp - r_ast_O(2)*st + r_ast_O(3)*ct*cp; 


%% 5. MTP to TP
% As in code by Amato, not sure of the reference

% Asymptotic velocity
GME = cons.GMe;
U = sqrt( V*V - 2*GME/b_p );
% Rescaling b-plane coordinates
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

v_b_pre = UVec /(DU/TU)

%% 6. Post-encounter coordinates
% Code from Amato, probably Valsecchi 2003
U_nd = U/(DU/TU);
m = cons.GMe/cons.GMs;

% t0   = et + dt_per ;
h    = 0;   % Number of revolutions until next encounter: only used for zeta2

[~,theta1,phi1,xi1,zeta1] = opik_next(U_nd,theta,phi,xi/DU,zeta/DU,t0,h,m);


%% 7. New heliocentric elements
% Equations in Valsecchi2015
ap    = kep_eat(1)/(1-kep_eat(2))/DU ;
longp = mod( kep_eat(4)+kep_eat(5)+kep_eat(6), 2*pi ) ; % In general sense should be longitude

kep_post_sma = opik_bplane_2_oe( theta1,phi1,zeta1,xi1,U_nd,phi,longp,ap );
kep_ast_post = kep_post_sma';
kep_ast_post(1) = kep_ast_post(1)*(1-kep_ast_post(2))*DU 



%% ========= Plotting for verification =========
% - 3D trajectory in orbit frame
% - 3D trajectory in inertial orbits

% Distance(t) - Validate delta t

figure(1); sc = cons.Re;
plot( dt_per/3600/24, norm(r_ast_O)/sc, 'ro' );
grid on
hold on
xlabel('dt (days)')
ylabel('d (R_E)')

tv = (-24*5:1:24*5) *3600 ;
d  = zeros(length(tv),1);
for i=1:length(tv)
    state_ast_t = cspice_conics(kep_ast_O, t0 + tv(i)); % Why negative?
    d(i) = norm(state_ast_t(1:3)) ;
    
    dr_ast_hyp(i,:) = state_ast_t(1:3);    
end

figure(1);
plot( tv/3600/24, d/sc, 'b' )


%% Compare to numerical integration
% t0 = et - 4*86400 ;
state_ast_t0 = cspice_conics(kep_ast_O, t0);
state_eat_t0 = cspice_conics(kep_eat,   t0);
X0 = HillRotInv( state_eat_t0, state_ast_t0 );

% X0 = cspice_conics(kep_ast, t0);

state_jup = cspice_spkezr( '5',   t0, 'ECLIPJ2000', 'NONE', '10' );
kep_jup = cspice_oscelt( state_jup, t0, cons.GMs );

tint = (0:1:24*8)*3600;

cons_ode.GMcentral = cons.GMs;
cons_ode.GMthird   = cons.GMe;
cons_ode.GMfourth  = 0;
cons_ode.OEp   = kep_eat;
cons_ode.OEp2  = kep_jup;
cons_ode.JD0_p = t0;

tol = 1e-13;
options=odeset('RelTol',tol,'AbsTol',ones(1,6)*tol);
[~,X]=ode113(@(t,X) multibodies_3rd4th_oe(t,X,cons_ode),tint,X0,options);
for i=1:length(tint)
        
    kep_4bp = cspice_oscelt( X(i,:)', t0 + tint(i), cons.GMs );
    kep_4bp_t(i,:) = kep_4bp;
    kep_4bp_t(i,1) = kep_4bp_t(i,1)/(1-kep_4bp_t(i,2));
    
    % moid_4bp_t(i) = MOID_SDG_win( kep0_4bp_t(i,[1 2 4 3 5]), kepE_sma([1 2 4 3 5]) );
    
    % Compute distance to Earth
    xe   = cspice_conics(kep_eat, tint(i)+t0 );
    d_int(i) = norm(xe(1:3) - X(i,1:3)'); % From 4bp integration
    
    state_ast_hel_t = cspice_conics(kep_ast, tint(i)+t0);
    state_eat_hel_t = cspice_conics(kep_eat, tint(i)+t0);

    d2(i)= norm(state_ast_hel_t(1:3) - state_eat_hel_t(1:3)); % From heliocentric elements
    
end

figure(1)
tp = tint ;%+ t0 - et ;
plot( tp(1)/86400, norm(state_ast_t0(1:3))/sc, 'go' )
plot( tp/86400, d_int/sc, 'm' )

state_ast_t0_hel = cspice_conics(kep_ast, t0);
plot( tp(1)/86400, norm(state_ast_t0_hel(1:3)-state_eat_t0(1:3))/sc, 'b+')

plot( tp/86400, d2/sc, 'c--' )

state_ast_0 = cspice_conics(kep_ast, et + dt);
state_eat_0 = cspice_conics(kep_eat, et + dt);
plot( 0/86400, norm(state_ast_0(1:3)-state_eat_0(1:3))/sc, 'k+')


%% Compare orbit elements
plot_oe_evolution( 'oe', tint, kep_4bp_t, 15, 'm', tint([1 end])/86400, 86400 )

tp = [0 dt_per/86400];
kep_ast_sma    = kep_ast; 
kep_ast_sma(1) = kep_ast_sma(1)/(1-kep_ast_sma(2))

kep_ast_plot = [kep_ast_sma kep_ast_sma]';

plot_oe_evolution( 'oe', tp, kep_ast_plot, 15, 'b--', tint([1 end])/86400, 1)


%% Do HillRot and HillRotInv go completely back and forth?

state_ast = cspice_conics(kep_ast, et );
state_eat = cspice_conics(kep_eat, et );

state_ast_O = HillRot(state_eat, state_ast);
state_ast_2 = HillRotInv(state_eat, state_ast_O);


%% How about full conversion between Geocentric and Heliocentric?
state_ast = cspice_conics(kep_ast, et );
state_eat = cspice_conics(kep_eat, et );

state_ast_O = HillRot(state_eat, state_ast);

kep_ast_O   = cspice_oscelt( state_ast_O, et, cons.GMe );

state_ast_O2 = cspice_conics( kep_ast_O, et );

state_ast_2 = HillRotInv( state_eat, state_ast_O2 );

kep_ast_2 = cspice_oscelt( state_ast_2, et, cons.GMs );
kep_ast;


%% Check initial heliocentric elements





















