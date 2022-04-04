% EFFECT OF NON-KEPLERIAN MOID EVOLUTION 
% ON PRELIMINARY KEYHOLE ANALYSES
%
% Continuing...
%
% SOLUTION OF A PLANETARY ENCOUNTER: PDC 2021
% 7th IAA Planetary Defense Conference 2021
clc
clear
% close all
format shortG

filename = 'main_alttarget.m';
filepath = matlab.desktop.editor.getActiveFilename;
currPath = filepath(1:(end-length(filename)));
cd(currPath)
addpath(genpath(pwd))


%% Read trajectory
% 1. Read heliocentric elements
cspice_furnsh( 'SPICEfiles/naif0012.tls.pc' )
cspice_furnsh( 'SPICEfiles/gm_de431.tpc' )
cspice_furnsh( 'SPICEfiles/pck00010.tpc' )
% cspice_furnsh( 'SPICEfiles/de431_part-1.bsp' )
% cspice_furnsh( 'SPICEfiles/de431_part-2.bsp' )
dir_local_de431 = 'C:\Users\Oscar\Documents\Spice-Kernels\';
cspice_furnsh( [dir_local_de431 'de431_part-1.bsp'] )
cspice_furnsh( [dir_local_de431 'de431_part-2.bsp'] )


cons.AU  = cspice_convrt(1, 'AU', 'KM');
cons.GMs = cspice_bodvrd( 'SUN', 'GM', 10 );
cons.GMe = cspice_bodvrd( '399', 'GM', 1 );
cons.Re  = 6378.140;

cons.yr  = 365.25 * 24 * 3600 ;
cons.Day = 3600*24; 

id = 1;
% load('Horizons-2016UD.mat')
load('Horizons-2022CJ5.mat')

% id = 62;
% epoch = '2014 October 01 TDB';
% epoch = '2014 August 01 TDB';

id = 74;
date_cal = char(CalendarDateTDB(id));
epoch = [date_cal(7:17) ' TDB'];
et = cspice_str2et( epoch );

state_eat = cspice_spkezr( '399', et, 'ECLIPJ2000', 'NONE', '10' );
kep_eat = cspice_oscelt( state_eat, et, cons.GMs );
sma_eat = kep_eat(1)/(1-kep_eat(2));

% kep_ast = cspice_oscelt( state_ast, et, cons.GMs );
kep_ast = [QR(id) EC(id) IN(id)*pi/180 OM(id)*pi/180 W(id)*pi/180 MA(id)*pi/180 et cons.GMs]';
sma_ast = kep_ast(1)/(1-kep_ast(2));

t0 = et ;
state_ast = cspice_conics(kep_ast, t0);
state_eat = cspice_conics(kep_eat, t0);

%%
% 3. Switch to geocentric and find time periapses
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
dt_per  = (MA_per - kep_ast_O(6))/n_ast_O; dt_per_day = dt_per/86400


% 4. Using ephemeris at date of perigee, obtain MTP coordinates
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


% 5. From MTP to TP
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


%% 1. Compute Circles

xi_nd   = xi/DU;
zeta_nd = zeta/DU;
% b_nd = sqrt(xi_nd^2 + zeta_nd^2);
b_nd = b/DU;
c_nd = (cons.GMe/cons.GMs)/U_nd^2;

b2 = b_nd*b_nd;
c2 = c_nd*c_nd;
U2 = U_nd*U_nd;
% aux1 = (b2+c2)*(1-U2);
% aux2 = -2*U_nd*( (b2-c2)*ct + 2*c_nd*zeta_nd*st );

a_pre  = 1/(1-U2-2*U_nd*ct);
% a_post = (b2 + c2)/(aux1 + aux2);

kmax = 20;
hmax = 50;

i = 0;
circ = [];
for k=1:kmax
    for h=1:hmax
        
        i = i+1 ;
        a0p(i) = (k/h)^(2/3);
        if sum( find(a0p(i) == a0p(1:i-1)) ); continue; end
        
        ct0p = ( 1-U2-ap/a0p(i) )/2/U_nd ;
        if abs(ct0p) > 1; continue; end
        st0p = sqrt(1 - ct0p^2);
        
        % Cirle center location and radius
        D = (c_nd*st)/(ct0p - ct);
        R = abs( c_nd*st0p/(ct0p - ct) );
        
%         % Skip if doesn't intersect the (0,Rstar) circle
        % Does not make sense becaues we don't design maneuvers
        % However, if we don't include a filter, the keyhole search
        % explodes for the very small circles
%         Rstar = 1*cons.Re/DU;
%         if ( (abs(D)<abs(R-Rstar)) || (abs(D)>(R+Rstar)) ); continue; end
        
        circ = [circ; k h D*DU R*DU];
        
    end
end


% Plot Valsecchi Circles
F = figure(2); clf;
plot_valsecchi_circles(2,circ,cons,U,kmax,xi,zeta); % first input, number of figure
axis([-25 25 -25 25])

%%
% plot( xi_nd, zeta_nd, 'r+')
i12 = find( circ(:,1)==12 );
for i=1:length(i12)
plot_circle( circ, i12(i), 'r-', cons.Re );
end

%% Plot horizons trajectory...
figure(3); clf;
OEplot = [A,EC,IN,OM,W,MA];
rng = 1:2000;
for i=1:6
subplot(2,3,i)
plot( JDTDB(rng)-JDTDB(1), OEplot(rng,i) ); grid on; hold on;
plot( JDTDB(id)-JDTDB(1), OEplot(id,i), 'r+' )
end

%% Candidates by list of flybys
load('sbdb_cad_flybys_2LD.mat')

reps = find( des(1:end-1)==des(2:end) );
repeated_bodies = des(reps);

T = table(des(reps),cd(reps),cd(reps+1),(jd(reps+1)-jd(reps))/365.25,...
    dist(reps)*cons.AU/cons.Re,v_inf(reps),h(reps),...
    'VariableNames',{'Des','Calendar Date 1','Calendar Date 2','diff (yr)',...
    'dCA (R_E)','vinf(km/s)','H'})



%% Convergence analysis of the B-plane coordinates
% Variants to test
% - Retrieve ephemeris in both cases
% - Making Earth circular along the way - Didn't work
% - Keep constant Earth elements starting way before
%   ~ Try ~SOI?

% Extras to add
% - Propagate trajectories
%  -> Nominal d(t)
%  -> d(t) with 2bp for both

iaux = 60:74;

% Defining a "nominal" time
id_nom = 73;
date_cal = char(CalendarDateTDB(id_nom));
epoch = [date_cal(7:17) ' TDB'];
et = cspice_str2et( epoch );

state_eat = cspice_spkezr( '399', et, 'ECLIPJ2000', 'NONE', '10' );
kep_eat = cspice_oscelt( state_eat, et, cons.GMs );
sma_eat = kep_eat(1)/(1-kep_eat(2));


%------------------------------
% OPTIONAL MAKE EARTH CIRCULAR
%------------------------------
% % Case make Earth circular at every time
% kep_eat(3) = 0;
% kep_eat(2) = 0;
% % kep_eat(1) = sma_eat;

for i = 1:length(iaux)
    
    id = iaux(i);
    
    date_cal = char(CalendarDateTDB(id));
    epoch = [date_cal(7:17) ' TDB'];
    et = cspice_str2et( epoch );
    
%     state_eat = cspice_spkezr( '399', et, 'ECLIPJ2000', 'NONE', '10' );
%     kep_eat = cspice_oscelt( state_eat, et, cons.GMs );
%     sma_eat = kep_eat(1)/(1-kep_eat(2));
%     
%     % Case make Earth circular at every time
%     kep_eat(3) = 0;
%     kep_eat(2) = 0;
% %     kep_eat(1) = sma_eat;
    
    % kep_ast = cspice_oscelt( state_ast, et, cons.GMs );
    kep_ast = [QR(id) EC(id) IN(id)*pi/180 OM(id)*pi/180 W(id)*pi/180 MA(id)*pi/180 et cons.GMs]';
    sma_ast = kep_ast(1)/(1-kep_ast(2));
    
    t0 = et ;
    
    F = figure(20); clf;
    [xiv(i), zetav(i), phiv(i), thetav(i), dtv(i), Uv(i), bv(i), dv(i)] = ...
        Bplane_encounter( kep_ast, kep_eat, t0, cons );
    
    grid on
    axis([-50 50 -50 50])

end

%%
F = figure(28); %clf;

subplot(2,2,1)
plot( JDTDB(iaux)-JDTDB(75), xiv/cons.Re )
grid on; hold on;
xlabel('t (days)'); ylabel('\xi (R_E)')

subplot(2,2,2)
plot( JDTDB(iaux)-JDTDB(75), zetav/cons.Re )
grid on; hold on;
xlabel('t (days)'); ylabel('\zeta (R_E)')

subplot(2,2,3)
plot( JDTDB(iaux)-JDTDB(75), phiv )
grid on; hold on;
xlabel('t (days)'); ylabel('\phi (R_E)')

subplot(2,2,4)
plot( JDTDB(iaux)-JDTDB(75), thetav )
grid on; hold on;
xlabel('t (days)'); ylabel('\theta (R_E)')


F = figure(29);

subplot(1,3,1)
plot( JDTDB(iaux)-JDTDB(75), JDTDB(iaux)'-JDTDB(75)+dtv)
grid on; hold on;
xlabel('t ephem (days)'); ylabel('CA Day')

subplot(1,3,2)
plot( JDTDB(iaux)-JDTDB(75), Uv )
grid on; hold on;
xlabel('t ephem (days)'); ylabel('U')

subplot(1,3,3)
plot( JDTDB(iaux)-JDTDB(75), bv/cons.Re )
grid on; hold on;
xlabel('t ephem (days)'); ylabel('b (R_E)')


% %%
% 
% figure(21);
% legend('E(t)','E(t)+circ','E(t0)','E(t0)+circ')
% 
% figure(22);
% legend('E(t)','E(t)+circ','E(t0)','E(t0)+circ')

%%

%% Figures wrt. Distance
F = figure(30); %clf;

subplot(2,2,1)
plot( dv/cons.AU, xiv/cons.Re )
grid on; hold on;
xlabel('d (au)'); ylabel('\xi (R_E)')

subplot(2,2,2)
plot( dv/cons.AU, zetav/cons.Re )
grid on; hold on;
xlabel('d (au)'); ylabel('\zeta (R_E)')

subplot(2,2,3)
plot( dv/cons.AU, phiv )
grid on; hold on;
xlabel('d (au)'); ylabel('\phi (R_E)')

subplot(2,2,4)
plot( dv/cons.AU, thetav )
grid on; hold on;
xlabel('d (au)'); ylabel('\theta (R_E)')


F = figure(31);

subplot(1,3,1)
plot( dv/cons.AU, JDTDB(iaux)'-JDTDB(75)+dtv)
grid on; hold on;
xlabel('d (au)'); ylabel('CA Day')

subplot(1,3,2)
plot( dv/cons.AU, Uv )
grid on; hold on;
xlabel('d (au)'); ylabel('U')

subplot(1,3,3)
plot( dv/cons.AU, bv/cons.Re )
grid on; hold on;
xlabel('d (au)'); ylabel('b (R_E)')



%% Propagate trajectory until the next encounter
% == Nominal ==
for id = 1:2450
    date_cal = char(CalendarDateTDB(id));
    epoch = [date_cal(7:17) ' TDB'];
    et = cspice_str2et( epoch );

    kep_ast = [QR(id) EC(id) IN(id)*pi/180 OM(id)*pi/180 W(id)*pi/180 MA(id)*pi/180 et cons.GMs]';
    state_ast = cspice_conics(kep_ast, et );
    state_eat = cspice_spkezr( '399', et, 'ECLIPJ2000', 'NONE', '10' );

    dist(id) = norm( state_ast(1:3) - state_eat(1:3) );
    
end

F = figure(25); clf;

plot( (JDTDB(1:2450) - JDTDB(1))/365.25, dist/cons.AU )
grid on; hold on;
xlabel('t (yr)')
ylabel('d (au)')


%% == Post encounter 1 set ==

% Defining a "nominal" time
id_nom = 72;
date_cal = char(CalendarDateTDB(id_nom));
epoch = [date_cal(7:17) ' TDB'];
et = cspice_str2et( epoch );

state_eat = cspice_spkezr( '399', et, 'ECLIPJ2000', 'NONE', '10' );
kep_eat = cspice_oscelt( state_eat, et, cons.GMs );
sma_eat = kep_eat(1)/(1-kep_eat(2));

% Case make Earth circular at every time
kep_eat(3) = 0;
kep_eat(2) = 0;
kep_eat(1) = sma_eat;

id = id_nom;
% kep_ast = cspice_oscelt( state_ast, et, cons.GMs );
kep_ast = [QR(id) EC(id) IN(id)*pi/180 OM(id)*pi/180 W(id)*pi/180 MA(id)*pi/180 et cons.GMs]';
sma_ast = kep_ast(1)/(1-kep_ast(2));

[d_pe, tvv] = distance_post( kep_ast, kep_eat, t0, cons );

figure(25); hold on;
plot( (tvv)/365.25/86400 + 0.5, d_pe/cons.AU )
xlim([-0.1 25])


%% Auxiliar functions
function plot_circle( circ, i, format, sc )

thv= 0:(2*pi/99):2*pi ;
    
D = circ(i,3);
R = circ(i,4);

xi_circ   = R*cos(thv);
zeta_circ = D + R*sin(thv);

plot( xi_circ/sc, zeta_circ/sc, format )

end

%=========================================
% Do everything in one...
function [xi, zeta, phi, theta, dt_per_day, U, b, b_p] = Bplane_encounter( kep_ast, kep_eat, t0, cons )

state_ast = cspice_conics(kep_ast, t0);
state_eat = cspice_conics(kep_eat, t0);

% 3. Switch to geocentric and find time periapses
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
dt_per  = (MA_per - kep_ast_O(6))/n_ast_O; dt_per_day = dt_per/86400


% 4. Using ephemeris at date of perigee, obtain MTP coordinates
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


% 5. From MTP to TP
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

% Compute Valsecchi Circles
xi_nd   = xi/DU;
zeta_nd = zeta/DU;
% b_nd = sqrt(xi_nd^2 + zeta_nd^2);
b_nd = b/DU;
c_nd = (cons.GMe/cons.GMs)/U_nd^2;

b2 = b_nd*b_nd;
c2 = c_nd*c_nd;
U2 = U_nd*U_nd;
% aux1 = (b2+c2)*(1-U2);
% aux2 = -2*U_nd*( (b2-c2)*ct + 2*c_nd*zeta_nd*st );

a_pre  = 1/(1-U2-2*U_nd*ct);
% a_post = (b2 + c2)/(aux1 + aux2);

kmax = 20;
hmax = 50;

i = 0;
circ = [];
for k=1:kmax
    for h=1:hmax
        
        i = i+1 ;
        a0p(i) = (k/h)^(2/3);
        if sum( find(a0p(i) == a0p(1:i-1)) ); continue; end
        
        ct0p = ( 1-U2-ap/a0p(i) )/2/U_nd ;
        if abs(ct0p) > 1; continue; end
        st0p = sqrt(1 - ct0p^2);
        
        % Cirle center location and radius
        D = (c_nd*st)/(ct0p - ct);
        R = abs( c_nd*st0p/(ct0p - ct) );
        
%         % Skip if doesn't intersect the (0,Rstar) circle
        % Does not make sense becaues we don't design maneuvers
        % However, if we don't include a filter, the keyhole search
        % explodes for the very small circles
%         Rstar = 1*cons.Re/DU;
%         if ( (abs(D)<abs(R-Rstar)) || (abs(D)>(R+Rstar)) ); continue; end
        
        circ = [circ; k h D*DU R*DU];
        
    end
end


% Plot the desired circles
i12 = find( circ(:,1)==12 );
for i=1:length(i12)
plot_circle( circ, i12(i), 'r-', cons.Re ); hold on;
end

plot( xi/cons.Re, zeta/cons.Re, 'k+' )

end

%===============================================
% Recycle code to propagate the trajectory
% Taken from main_CMDA_modules point 5
function [d_pe, tvv] = distance_post( kep_ast, kep_eat, t0, cons )

state_ast = cspice_conics(kep_ast, t0);
state_eat = cspice_conics(kep_eat, t0);

% 3. Switch to geocentric and find time periapses
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
dt_per  = (MA_per - kep_ast_O(6))/n_ast_O; dt_per_day = dt_per/86400


% 4. Using ephemeris at date of perigee, obtain MTP coordinates
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


% 5. From MTP to TP
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

xi0 = xi/DU;
zeta0 = zeta/DU;

% 6. Solving the encounter with Opik formulae
h = 0; % Number of revolutions until next encounter: only used for zeta2

% Using keyhole point selected: xi0, zeta0
[zeta2,theta1,phi1,xi1,zeta1,r_b_post,v_b_post] = opik_next(U_nd,theta,phi,xi0,zeta0,t0,h,m);

% 7. Heliocentric orbit elements from post-encounter coordinates (Opik formulae)
longp = mod( kep_eat(4)+kep_eat(5)+kep_eat(6), 2*pi ) ; % In general sense should be longitude
longp = atan2( state_eat(2),state_eat(1) );
kep_opik_post = opik_bplane_2_oe( theta1,phi1,zeta1,xi1,U_nd,phi,longp,ap )';

kep_opik_post(1) = kep_opik_post(1)*(1-kep_opik_post(2))*DU ; 
% 'opik_bplane_2_oe.m' function returns sma in first element

kep_opik_post(6) = TA_2_MA(kep_opik_post(6),kep_opik_post(2));
% 'opik_bplane_2_oe.m' function returns true anomaly 6th element

kep_opik_post(7:8) = [t0 + dt_per;
                        cons.GMs];
% Include GM of the central body and epoch


% Post-encounter constant elements distance
tf = 20;
tv1 = ((-0.1:.001:tf) + 0) *cons.yr;
et0 = t0 + dt_per;

d_pe = zeros(length(tv1),1);
for i=1:length(tv1)
    xa   = cspice_conics(kep_opik_post, et0 + tv1(i) );
    xe   = cspice_conics(kep_eat,       et0 + tv1(i) );
    d_pe(i) = norm(xe(1:3) - xa(1:3)); 
end
tvv = tv1;

end


%% 















