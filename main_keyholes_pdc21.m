% GENERATION OF MODIFIED KEYHOLES ON THE B-PLANE
% 7th IAA Planetary Defense Conference 2021
clc
clear
close all
addpath(genpath(pwd))

%% Script steps
% 1. States of Earth and Asteroid
% 2. Find crossing of SOI
% 3. Coordinates in Earth-Centered Orbital Frame
% 4. Opik Variables at CA
% 5. Scan for resonances (Compute Valsecchi circles)
% 6. Compute Keyholes

%% 0. Constants

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

asteroid = 2; 
%{ 1. 2021 PDC, 2. 2017 PDC }

switch asteroid
    case 1
        cspice_furnsh( '2021_PDC-s11-merged-DE431.bsp' )
        ast_id= '-937014';
%         epoch = '2021 April 20 TDB';
        epoch = '2021 October 18 TDB';
        tv = (0:.25:24*2.2) *3600 ;
        
    case 2
        cspice_furnsh( '2017_PDC-merged-DE431.bsp' )
        ast_id= '-937001';
%         epoch = '2027 July 20, 13:09:38.000 TDB';
%         epoch = '2027 July 21, 08:17:15.999 TDB';
        epoch = '2027 July 19 TDB';
        tv = (0:.25:24*2.2) *3600 ;
        
end


%% 1. States of Earth and Asteroid
%-- Setup for 2021 PDC

et = cspice_str2et( epoch );
EQ_2_EC = cspice_sxform( 'J2000', 'ECLIPJ2000', et );

cons.AU  = cspice_convrt(1, 'AU', 'KM');
cons.GMs = cspice_bodvrd( 'SUN', 'GM', 1 );
cons.GMe = cspice_bodvrd( '3', 'GM', 1 );
cons.Re  = mean( cspice_bodvrd( '399', 'RADII', 3 ) );
cons.yr  = 365.25 * 24 * 3600 ;
cons.Day = 3600*24; 

state_eat = cspice_spkezr( '3', et, 'J2000', 'NONE', '0' );
state_eat = EQ_2_EC * state_eat;
kep_eat = cspice_oscelt( state_eat, et, cons.GMs );

state_ast = cspice_spkezr( ast_id, et, 'J2000', 'NONE', '0' );
state_ast = EQ_2_EC * state_ast;
kep_ast = cspice_oscelt( state_ast, et, cons.GMs );

% Orbital periods
a_eat = kep_eat(1)/(1-kep_eat(2)) ;
T_eat = 2*pi / sqrt(  cons.GMs/a_eat^3 ) ;

a_ast = kep_ast(1)/(1-kep_ast(2)) ;
T_ast = 2*pi / sqrt(  cons.GMs/a_ast^3 ) ;

% Reference units (km,km/s)
% DU: Distance of the Earth from the Sun
% TU: Such that mu = 1
DU = norm( state_eat(1:3) );
TU = sqrt(DU^3/cons.GMs);

state_ast_nd = state_ast/DU ;
state_ast_nd(4:6) = state_ast(4:6)*TU ;
kep_ast_nd = [kep_ast(1)/DU; kep_ast(2:end)] ;
T_ast_nd = T_ast/TU ;

state_eat_nd = state_eat/DU ;
state_eat_nd(4:6) = state_eat(4:6)*TU ;
kep_eat_nd = [kep_eat(1)/DU; kep_eat(2:end)] ;
T_eat_nd = T_eat/TU ;


%% 2. Finding crossing of Sphere of Influence
d_eph = zeros(length(tv),1);

SOI = a_eat*(cons.GMe/cons.GMs)^(2/5);

for i=1:length(tv) % Limit of ephemeris data
    
    state_ast_t = cspice_spkezr( ast_id, et+tv(i), 'J2000', 'NONE', '0' );   
    state_eat_t = cspice_spkezr( '3', et+tv(i), 'J2000', 'NONE', '0' );
    
    d_eph(i) = norm(state_ast_t(1:3) - state_eat_t(1:3));
    
end

F = figure(1);
plot( tv/3600/24, d_eph )
grid on
hold on
xlabel('dt (days)')
ylabel('d (km)')

id_SOI = find( d_eph<SOI,1 );
dt_SOI = (SOI - d_eph(id_SOI-1)) * (tv(id_SOI)-tv(id_SOI-1))/(d_eph(id_SOI)-d_eph(id_SOI-1)) ;
et_SOI = et + dt_SOI + tv(id_SOI-1);

% hold on
% plot( (dt_SOI+tv(id_SOI-1))/3600/24, SOI, 'ro' )
% plot( (tv(id_SOI))/3600/24, d_eph(id_SOI), 'r+' )
% plot( (tv(id_SOI-1))/3600/24, d_eph(id_SOI-1), 'r+' )


%% 3. Coordinates in the Earth-Centered orbital frame
state_eat = cspice_spkezr( '3',    et_SOI, 'J2000', 'NONE', '0' );
state_ast = cspice_spkezr( ast_id, et_SOI, 'J2000', 'NONE', '0' );

[NO, NO6] = HillRot_DCM(state_eat);
state_ast_O = NO'*(state_ast(1:3) - state_eat(1:3));
state_ast_O(4:6) = NO'*(state_ast(4:6) - state_eat(4:6));

r_ast_O = state_ast_O(1:3);
v_ast_O = state_ast_O(4:6);

% Unperturbed prop from SOI crossing to pericenter
tv = 0;
tv = (0:.1:36) *3600 ;
d  = zeros(length(tv),1);

kep_ast_O = cspice_oscelt( state_ast_O, et_SOI, cons.GMe );
for i=1:length(tv)
    
    state_ast_t = cspice_conics(kep_ast_O, et_SOI + tv(i));
    d(i) = norm(state_ast_t(1:3)) ;
    
end

plot( (tv+et_SOI-et)/3600/24, d )

MA_per = 0;
a_ast_O = kep_ast_O(1)/(1 - kep_ast_O(2));
n_ast_O = sqrt( -cons.GMe/a_ast_O(1)^3 );
dt_per  = (MA_per - kep_ast_O(6))/n_ast_O;

state_ast_per = cspice_conics( kep_ast_O, et_SOI + dt_per );
plot( (et_SOI+dt_per-et)/3600/24, norm(state_ast_per), 'r+' )


%% 4. Opik variables at b-plane
U     = norm( v_ast_O );
theta = acos( v_ast_O(2)/U );
phi   = atan2( v_ast_O(1), v_ast_O(3) );

h = cross( r_ast_O, v_ast_O ) ;
b = norm(h)/U ;

% Assume linear motion to find b-plane crossing
xi   = r_ast_O(1)*cos(phi);
zeta = r_ast_O(1)*cos(theta)*sin(phi) - r_ast_O(2)*sin(theta);

%--------------------------------
% There is something wrong here!
%--------------------------------

%% 5. Scan for resonances
xi_nd = xi/DU;
zeta_nd = zeta/DU;
b_nd = sqrt(xi_nd^2 + zeta_nd^2);
U_nd = U/(DU/TU);
c_nd = (cons.GMe/cons.GMs)/U_nd^2;
ct = cos(theta);
st = sin(theta);

b2 = b_nd*b_nd;
c2 = c_nd*c_nd;
U2 = U_nd*U_nd;
aux1 = (b2+c2)*(1-U2);
aux2 = -2*U_nd*( (b2-c2)*ct + 2*c_nd*zeta_nd*st );

a_pre  = 1/(1-U2-2*U_nd*ct);
a_post = (b2 + c2)/(aux1 + aux2);

kmax = 20;
hmax = 50;

i = 0;
circ = [];
for k=1:kmax
    for h=1:hmax
        
        i = i+1 ;
        a0p(i) = (k/h)^(2/3);
        if sum( find(a0p(i) == a0p(1:i-1)) ); continue; end
        
        ct0p = ( 1-U2-1/a0p(i) )/2/U_nd ;
        if abs(ct0p) > 1; continue; end
        st0p = sqrt(1 - ct0p^2);
        
        % Cirle center location and radius
        D = (c_nd*st)/(ct0p - ct);
        R = abs( c_nd*st0p/(ct0p - ct) );
        
        
%         % Skip if doesn't intersect the (0,Rstar) circle
%         % Does not make sense becaues we don't design maneuvers
%         if ( (abs(D)<abs(R-Rstar)) || (abs(D)>(R+Rstar)) ); continue; end
        
        circ = [circ; 
            k h D*DU R*DU];
        
    end
end

% Plot Valsecchi circles!
F  = figure(2);

nr = size(circ,1);
co = winter(22);
thv= 0:(2*pi/99):2*pi ;
RE_km = 6378.140;
sc = RE_km;

focus_factor = sqrt(1 + 2*cons.GMe/(RE_km*U^2));
RE_focussed  = focus_factor;

for i=1:nr
    
    k = circ(i,1);
    D = circ(i,3);    
    R = circ(i,4);
    
    xi_circ   = R*cos(thv);
    zeta_circ = D + R*sin(thv);
    
    cc = co(k,:);
    plot( xi_circ/sc, zeta_circ/sc, 'Color', cc )
    hold on
    
end

fill(RE_focussed*cos(thv), RE_focussed*sin(thv),'white');
plot(RE_focussed*cos(thv), RE_focussed*sin(thv),'k');
plot(cos(thv), sin(thv),'k--');

grid on
axis equal
caxis([1 20])
cb = colorbar;
cb.Label.String = 'k';
axis([-1 1 -1 1]*10)
xlabel('\xi (R_\oplus)');
ylabel('\zeta (R_\oplus)');


%% 6. Keyhole computation
% Section dependencies: scripts in 'keyholes'

t0 = 0;
m  = cons.GMe/cons.GMs ;

F = figure(3);
hold on

sc = 6378.140/DU;
for i=1:nr
    k = circ(i,1);
    h = circ(i,2);
    D = circ(i,3)/RE_km;    
    R = circ(i,4)/RE_km;    
    
    [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
        two_keyholes(k, h, D, R, U_nd, theta, phi, m,t0,DU);
    
    cc = co(k,:);    
    plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
        'Color',cc);
    plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
         'Color',cc);    
    
end

fill(RE_focussed*cos(thv), RE_focussed*sin(thv),'white');
plot(RE_focussed*cos(thv), RE_focussed*sin(thv),'k');
plot(cos(thv), sin(thv),'k--');

grid on
axis equal
caxis([1 20])
cb = colorbar;
cb.Label.String = 'k';
% axis([-1 1 -1 1]*10)
xlabel('\xi (R_\oplus)');
ylabel('\zeta (R_\oplus)');












