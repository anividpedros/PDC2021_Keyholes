% GENERATION OF MODIFIED KEYHOLES ON THE B-PLANE
% 7th IAA Planetary Defense Conference 2021
clc
clear
close all

filename = 'main_keyholes_pdc21.m';
filepath = matlab.desktop.editor.getActiveFilename;
currPath = filepath(1:(end-length(filename)));
cd(currPath)
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

cons.AU  = cspice_convrt(1, 'AU', 'KM');
cons.GMs = cspice_bodvrd( 'SUN', 'GM', 10 );
% cons.GMe = cspice_bodvrd( '3', 'GM', 1 );
cons.GMe = 398600.43543609593;

% cons.Re  = mean( cspice_bodvrd( '399', 'RADII', 3 ) );
cons.Re  = 6378.140;

cons.yr  = 365.25 * 24 * 3600 ;
cons.Day = 3600*24; 

state_eat = cspice_spkezr( '399', et, 'ECLIPJ2000', 'NONE', '10' );
kep_eat = cspice_oscelt( state_eat, et, cons.GMs );

state_ast = cspice_spkezr( ast_id, et, 'ECLIPJ2000', 'NONE', '10' );
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
    
    state_ast_t = cspice_spkezr( ast_id, et+tv(i), 'ECLIPJ2000', 'NONE', '10' );   
    state_eat_t = cspice_spkezr( '399',  et+tv(i), 'ECLIPJ2000', 'NONE', '10' );
    
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
% Take same date as in script by Amato
if asteroid == 2
    epoch_inp   = '2027 July 21, 08:17:15.999 TDB';
    et_inp = cspice_str2et( epoch_inp );
    et_SOI = et_inp;
    dt_SOI = et_SOI - et;
end

state_eat = cspice_spkezr( '399',  et_SOI, 'ECLIPJ2000', 'NONE', '10' );
state_ast = cspice_spkezr( ast_id, et_SOI, 'ECLIPJ2000', 'NONE', '10' );
state_ast_eat = cspice_spkezr( ast_id, et_SOI, 'ECLIPJ2000', 'NONE', '399' );

DU = norm( state_eat(1:3) );
TU = sqrt(DU^3/cons.GMs);

state_ast_O = HillRot(state_eat, state_ast);
r_ast_O = state_ast_O(1:3);
v_ast_O = state_ast_O(4:6);
plot( (et_SOI-et)/3600/24, norm(r_ast_O), 'ro' );


% Unperturbed prop from SOI crossing to pericenter
kep_ast_O = cspice_oscelt( state_ast_O, et_SOI, cons.GMe );

% Plot coordinates given by hyperbolic elements
tv = (-12:.1:36) *3600 ;
d  = zeros(length(tv),1);
for i=1:length(tv)
    state_ast_t = cspice_conics(kep_ast_O, et_SOI + tv(i));
    d(i) = norm(state_ast_t(1:3)) ;
end
plot( (tv+et_SOI-et)/3600/24, d )


% Coordinates at perigee
MA_per = 0;
% Coordinates at node crossing
% TA_node = -kep_ast_O(5);
% MA_node = TA_2_MA( TA_node, kep_ast_O(2) );

a_ast_O = kep_ast_O(1)/(1 - kep_ast_O(2));
n_ast_O = sqrt( -cons.GMe/a_ast_O(1)^3 );
dt_per  = (MA_per - kep_ast_O(6))/n_ast_O;

state_ast_per = cspice_conics( kep_ast_O, et_SOI + dt_per );
r_ast_O = state_ast_per(1:3);
v_ast_O = state_ast_per(4:6);

plot( (et_SOI+dt_per-et)/3600/24, norm(state_ast_per), 'r+' )


%% 4. Opik variables at b-plane
% Coordinates at MTP - Perigee
b     = norm( r_ast_O );
V     = norm( v_ast_O );
theta = acos( v_ast_O(2)/V );
phi   = atan2( v_ast_O(1), v_ast_O(3) );

cp = cos(phi);   sp = sin(phi);
ct = cos(theta); st = sin(theta);

Rphi = rotationM(-phi,2);
Rth  = rotationM(-theta,1);
ROT  = Rth*Rphi ;
r_ast_b = ROT*r_ast_O;

xi   = r_ast_b(1);
eta  = r_ast_b(2);
zeta = r_ast_b(3);

% CODE from processing.f90/CART_TO_OPIK does the rotation with +theta, +phi
% verified by changing the sign in the call of rotationM for Rphi,Rth
% Valsecchi[2003] does the rotation by -theta, -phi!
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];
xi   = r_ast_O(1)*cp - r_ast_O(3)*sp;
eta  = r_ast_O(1)*st*sp + r_ast_O(2)*ct + r_ast_O(3)*st*cp;
zeta = r_ast_O(1)*ct*sp - r_ast_O(2)*st + r_ast_O(3)*ct*cp; 
% At the moment, when trying to replicate the results, this is left as is.
% And it makes sense, the result is eta = 0!

% Coordinates at TP

% Asymptotic velocity
GME = cons.GMe;
U = sqrt( V*V - 2*GME/b );
% Rescaling b-plane coordinates
b    = (V/U)*b  ;
xi   = (V/U)*xi ;
zeta = (V/U)*zeta ;
% Rotation of velocity vector
hgamma = atan( cons.GMe/(b*U*U) );
hv     = cross( r_ast_O, v_ast_O );
DCM    = PRV_2_DCM( -hgamma, hv/norm(hv) );
UVec   = (U/V)*DCM*v_ast_O;
theta  = acos(UVec(2)/U);
phi    = atan2(UVec(1),UVec(3));

cp = cos(phi);   sp = sin(phi);
ct = cos(theta); st = sin(theta);


%% 5. Scan for resonances
xi_nd   = xi/DU;
zeta_nd = zeta/DU;
b_nd = sqrt(xi_nd^2 + zeta_nd^2);
U_nd = U/(DU/TU);
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
        
        ct0p = ( 1-U2-1/a0p(i) )/2/U_nd ;
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
        
        circ = [circ; 
            k h D*DU R*DU];
        
    end
end

% Plot Valsecchi circles!
F  = figure(2);

nr = size(circ,1);
co = winter(22);
thv= 0:(2*pi/99):2*pi ;
sc = cons.Re;

focus_factor = sqrt(1 + 2*cons.GMe/(cons.Re*U^2));
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
colormap(co);
fill(RE_focussed*cos(thv), RE_focussed*sin(thv),'white');
plot(RE_focussed*cos(thv), RE_focussed*sin(thv),'k');
plot(cos(thv), sin(thv),'k--');
plot(3*cos(thv), 3*sin(thv),'k--');

grid on
axis equal
caxis([1 20])
cb = colorbar;
cb.Label.String = 'k';
axis([-1 1 -1 1]*5)
xlabel('\xi (R_\oplus)');
ylabel('\zeta (R_\oplus)');

%% Compare results to file
T = readtable('2017PDC_resonances_all.dat');

circ2 = circ;
circ2(:,3:4) = circ2(:,3:4)/sc ;


%% 6. Keyhole computation
% Section dependencies: scripts in 'keyholes'

RE_au = cons.Re/DU;

%---- Choose new circles ----
m  = cons.GMe/cons.GMs ;
t0 = 0;
circles = circ;

%---- Choose precomputed circles ----
% % Nominal Öpik parameters
% U_nd = 0.44705634713083897;
% theta = 6.573081547432976E+01*pi/180;
% phi = 7.141024642502616E+01*pi/180;
% m = 3.0034896149157654e-06;
% t0 = 0;
% DU = 152003856.97586098;
% 
% circles = T{:,:};
% circles(:,3:4) = circles(:,3:4)*cons.Re;
% ------------------------------------

F = figure(3);
hold on

sc = cons.Re/DU;
nr = size(circles,1);
for i=1:nr
    
    % New circles
    k = circles(i,1);
    h = circles(i,2);
    D = circles(i,3)/cons.Re;    
    R = circles(i,4)/cons.Re;    
    
    [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
        two_keyholes(k, h, D, R, U_nd, theta, phi, m,t0,DU);
    
    cc = co(k,:);    
    plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
        'Color',cc);
    plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
         'Color',cc);    
    
end
colormap(co);
fill(RE_focussed*cos(thv), RE_focussed*sin(thv),'white');
plot(RE_focussed*cos(thv), RE_focussed*sin(thv),'k');
plot(cos(thv), sin(thv),'k--');
plot(3*cos(thv), 3*sin(thv),'k--');

grid on
axis equal
caxis([1 20])
cb = colorbar;
cb.Label.String = 'k';
axis([-1 1 -1 1]*30)
xlabel('\xi (R_\oplus)');
ylabel('\zeta (R_\oplus)');



%% Example of MOID variation over time
%% Pick flyby from the keyholes and generate ICs

kref = 15;
ic = find( circles(:,1) == kref, 1 );
ic = 142;

k = circles(ic,1);
h = circles(ic,2);
D = circles(ic,3)/cons.Re;
R = circles(ic,4)/cons.Re;

[kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
    two_keyholes(k, h, D, R, U_nd, theta, phi, m,t0,DU);

% 1. Find a point close to xi=0 (arbitrary choice)
% [~,ik] = min( abs(kh_up_xi) );
% 2. Find the first point of the keyhole arc (arbirary as well)
ik = find( ~isnan(kh_up_xi), 1 ); 

xi0   = kh_up_xi(ik);
zeta0 = kh_up_zeta(ik);

% Find post-encounter Opik coordinates
[zeta2,theta1,phi1,xi1,zeta1] = opik_next(U_nd,theta,phi,xi0(1),zeta0(1),t0,h);

% Post-encounter Opik to Cartesisan
cp = cos(phi1);   sp = sin(phi1);
ct = cos(theta1); st = sin(theta1);
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp]';

Rphi = rotationM(-phi1,2);
Rth  = rotationM(-theta1,1);
auxR  = Rth*Rphi ;

r_ast_O1 = DU*auxR'*[xi1; 0; zeta1]
v_ast_O1 = U_nd*(DU/TU)*[st*sp; ct; st*cp]

% Assuming the encounter is instantaneous, et_1 = et_SOI + dt_per
% Generate heliocentric elements
state_ast_post = [r_ast_O1; v_ast_O1];
state_eat = cspice_spkezr( '399',  et_SOI+dt_per, 'ECLIPJ2000', 'NONE', '10' );

state_ast_post_sun = HillRotInv( state_eat, state_ast_post );
kep_ast_post = cspice_oscelt( state_ast_post_sun, et_SOI+dt_per, cons.GMs )

% Propagation
et0      = et_SOI+dt_per;
kep0     = kep_ast_post;
kep0_sma = kep_ast_post';
kep0_sma(1) = kep0(1)/(1-kep0(2));

kepE = cspice_oscelt( state_eat, et_SOI+dt_per, cons.GMs );
kepE_sma = kepE';
kepE_sma(1) = kepE(1)/(1-kepE(2));


MOID0 = MOID_SDG_win( kep0_sma([1 2 4 3 5]), kepE_sma([1 2 4 3 5]) );



state_jup = cspice_spkezr( '5',  et_SOI+dt_per, 'ECLIPJ2000', 'NONE', '10' );
kepJ = cspice_oscelt( state_jup, et_SOI+dt_per, cons.GMs );
kepJ_sma = kepJ;
kepJ_sma(1) = kepJ(1)/(1-kepJ(2));

% Secular Model: Lagrange-Laplace
cons.OEp = kepJ';
cons.GMp = cspice_bodvrd( '5', 'GM', 1 );
secular_model_LL = secular_model_10BP_s2(kep0_sma, cons, 1);
tv = (0.:.05:50) *cons.yr;
for i = 1:length(tv)
    
    [~,kep0_LL_t(i,:)] = drifted_oe_s2( secular_model_LL, tv(i), kep0_sma, kepJ' );
    moid_LL_t(i) = MOID_SDG_win( kep0_LL_t(i,[1 2 4 3 5]), kepE_sma([1 2 4 3 5]) );
    
end

% Numerical integration: 3rd-body perturbers: Earth and Jupiter
% X0 = state_ast_post;
X0 = cspice_conics(kep0, et0+tv(1) );

cons_ode.GMcentral = cons.GMs;
cons_ode.GMthird  = cons.GMe;
cons_ode.GMfourth = cons.GMp(1);
cons_ode.OEp   = kepE;
cons_ode.OEp2  = kepJ;
cons_ode.JD0_p = et0;

tol = 1e-13;
options=odeset('RelTol',tol,'AbsTol',ones(1,6)*tol);
[t,X]=ode113(@(t,X) multibodies_3rd4th_oe(t,X,cons_ode),tv,X0,options);
for i=1:length(tv)
    kep0_4bp = cspice_oscelt( X(i,:)', et0+tv(i), cons.GMs );
    kep0_4bp_t(i,:) = kep0_4bp;
    kep0_4bp_t(i,1) = kep0_4bp_t(i,1)/(1-kep0_4bp_t(i,2));
    
    moid_4bp_t(i) = MOID_SDG_win( kep0_4bp_t(i,[1 2 4 3 5]), kepE_sma([1 2 4 3 5]) );
    
    % Compute distance to Earth
    xe   = cspice_conics(kepE, et0+tv(i) );
    d(i) = norm(xe(1:3) - X(i,1:3)'); % From 4bp integration
    
    xa   = cspice_conics(kep0, et0+tv(i) );
    d2(i)= norm(xe(1:3) - xa(1:3)); % From post-encounter elements
    
end

% Plot moid time evolutions
F = figure(5);
clf

xsc = cons.yr;
ysc = cons.Re;

plot(tv([1 end])/xsc, MOID0*[1 1]/ysc, '--')
hold on
plot(tv/xsc, moid_LL_t/ysc)
plot(tv/xsc, moid_4bp_t/ysc)
% plot(tv/xsc, d/ysc)
% plot(tv/xsc, d2/ysc)

grid on
xlabel('time (yr)')
ylabel('MOID (R_\oplus)')

legend('post-encounter','secLL','4BP')

% Is the next encounter happening?
F = figure(6);
clf

plot(tv/xsc, d/ysc)
hold on
plot(tv/xsc, d2/ysc)

grid on
xlabel('time (yr)')
ylabel('d (R_\oplus)')

legend('4BP','post-enc')


%% New Keyhole radius with that MOID:
F = figure(10); clf;
hold on

% kref = 15; % takes 141:158
% ic = find( circles(:,1) == kref, 1 );

profile clear
profile on


co = summer(22);
sc = cons.Re/DU;
nr = size(circles,1);
% for i=find( circles(:,1) == kref )
for i= 1:nr %142%141:158 %ic%:nr
    
    % New circles
    k = circles(i,1);
    h = circles(i,2);
    D = circles(i,3)/cons.Re;    
    R = circles(i,4)/cons.Re;    
    
    it = find( tv/cons.yr == k );
    
    [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
        two_keyholes(k, h, D, R, U_nd, theta, phi, m,t0,DU );

    cc = co(k,:);
    plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
        'Color',cc, 'LineWidth',2.5);
    plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
         'Color',cc, 'LineWidth',2.5);    
    
    % We are assuming for this preliminar computation that the same
    % deltaMOID happens for all the points of the circle
    dxi_num = moid_4bp_t(it)/DU - xi1;
    dxi_LL  = moid_LL_t(it)/DU - xi1;
    
    [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
        two_keyholes_dxi(k, h, D, R, U_nd, theta, phi, m,t0,DU, dxi_num);
    
    plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
        'r');
    plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
        'r');   

     [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
        two_keyholes_dxi(k, h, D, R, U_nd, theta, phi, m,t0,DU, dxi_LL);
    
    plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
        'b');
    plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
        'b');   
    
    
end

colormap(co);
fill(RE_focussed*cos(thv), RE_focussed*sin(thv),'white');
plot(RE_focussed*cos(thv), RE_focussed*sin(thv),'k');
plot(cos(thv), sin(thv),'k--');
plot(3*cos(thv), 3*sin(thv),'k--');

grid on
axis equal
caxis([1 20])
cb = colorbar;
cb.Label.String = 'k';
% axis([-1 1 -1 1]*5)
xlabel('\xi (R_\oplus)');
ylabel('\zeta (R_\oplus)');

% xlim([-20 20 ])

profile viewer
profile off