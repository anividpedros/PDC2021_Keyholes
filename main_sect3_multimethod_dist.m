% SOLUTION OF A PLANETARY ENCOUNTER: PDC 2021
% 7th IAA Planetary Defense Conference 2021
clc
clear
close all
format shortG

filename = 'main_sect3_multimethod_dist.m';
filepath = matlab.desktop.editor.getActiveFilename;
currPath = filepath(1:(end-length(filename)));
cd(currPath)
addpath(genpath(pwd))

%% ===== Description ===== 
% Original script: main_1flyby_pdc2021.m
% Extension: Generate distance(t) plots given post-encounter coordinates
%  =======================

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
sma_ast = kep_ast(1)/(1-kep_ast(2));

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


%% 6. Solving the encounter with Opik formulae
h = 0; % Number of revolutions until next encounter: only used for zeta2

[~,theta1,phi1,xi1,zeta1,r_b_post,v_b_post] = opik_next(U_nd,theta,phi,xi/DU,zeta/DU,t0,h,m);


%% 7. Heliocentric orbit elements from post-encounter coordinates (Opik formulae)
kep_opik_post = opik_bplane_2_oe( theta1,phi1,zeta1,xi1,U_nd,phi,longp,ap )';
kep_opik_post(1) = kep_opik_post(1)*DU ;
% This function returns sma in first element


%% Reading data from REBOUND
dades = csvread(['PDCasteroid']);
        t = dades(1,:);
        d = dades(2,:);
        
F = figure(4);
plot(t/2/pi,d)
xlabel('time (yr)')
ylabel('d (au)')
grid on
        

% Compute MOID of integrated trajectory
nt = size(dades,2);
% IMPLEMENT TIME VARYING ELEMENTS OF EARTH
for i=1:nt
    kep_nbp = dades(3:8,i)';
    moid_nbp_t(i) = MOID_SDG_win( kep_nbp([1 2 4 3 5]), kepE_sma([1 2 4 3 5]) );
end
figure
plot(t/2/pi,moid_nbp_t)
grid on

%% Secular propagation
kep0 = kep_opik_post;
kep0_sma = kep_opik_post';
kep0_sma(1) = kep0(1)/(1-kep0(2));

kepE = kep_eat;
kepE_sma = kepE';
kepE_sma(1) = kepE(1)/(1-kepE(2));

MOID0 = MOID_SDG_win( kep0_sma([1 2 4 3 5]), kepE_sma([1 2 4 3 5]) );

state_jup = cspice_spkezr( '5',  et0, 'ECLIPJ2000', 'NONE', '10' );
kepJ = cspice_oscelt( state_jup,   et0, cons.GMs );
kepJ_sma = kepJ;
kepJ_sma(1) = kepJ(1)/(1-kepJ(2));

% Secular Model: Lagrange-Laplace
cons.OEp = kepJ';
cons.GMp = cspice_bodvrd( '599', 'GM', 1 );
secular_model_LL = secular_model_10BP_s2(kep0_sma, cons, 1);
tv = (0.:.001:20) *cons.yr;

kep_LL_t = kep_opik_post;

for i = 1:length(tv)
    
    [~,kep0_LL_t(i,:)] = drifted_oe_s2( secular_model_LL, tv(i), kep0_sma, kepJ' );
    kep_LL_t(2:6) = kep0_LL_t(i,2:6); 
    
    xa   = cspice_conics(kep_LL_t, et0 + tv(i) );
    xe   = cspice_conics(kep_eat,  et0 + tv(i) );
    d2(i) = norm(xe(1:3) - xa(1:3)); 
    
    moid_LL_t(i) = MOID_SDG_win( kep0_LL_t(i,[1 2 4 3 5]), kepE_sma([1 2 4 3 5]) );
    
end

% Plot moid time evolutions
F = figure(5);
clf

xsc = cons.yr;
ysc = DU; %cons.Re;

plot(tv([1 end])/xsc, MOID0*[1 1]/ysc, '--')
hold on
plot(tv/xsc, moid_LL_t/ysc)
% plot(tv/xsc, moid_4bp_t/ysc)
% plot(tv/xsc, d/ysc)
% plot(tv/xsc, d2/ysc)

grid on
xlabel('time (yr)')
ylabel('MOID (DU)');%(R_\oplus)')

legend('post-encounter','secLL')%,'4BP')

%-------------------------
% Is the next encounter happening?
F = figure(6);

% plot(tv/xsc, d/ysc)
hold on
plot(tv/xsc, d2/ysc)

grid on
xlabel('time (yr)')
ylabel('d (R_\oplus)')

legend('4BP','post-enc')
