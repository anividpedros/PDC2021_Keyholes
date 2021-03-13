% GENERATION OF MODIFIED KEYHOLES ON THE B-PLANE
% 7th IAA Planetary Defense Conference 2021
clc
clear
close all
format shortG

filename = 'main_keyholes_pdc21_detailsFlyby_testing.m';
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

% 2021 PDC

cspice_furnsh( '2021_PDC-s11-merged-DE431.bsp' )
ast_id= '-937014';
epoch = '2021 October 20 TDB';


%% 1. States of Earth and Asteroid
%-- Setup for 2021 PDC

et = cspice_str2et( epoch ) + 24*.715*3600;

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


% 1. Propagation to find dCA with modified elements
% Modifying the elements

% Make Earth Circular-Ecliptic -----
kep_eat(2:3) = 0;
state_eat = cspice_conics(kep_eat, et);
% Modify NEO to try to make dCA smaller ----
kep_ast(3) = kep_ast(3)+.2;     % Increase inclination
kep_ast(5) = kep_ast(5)+.06;    % Adjust arg of perihelion to low MOID
% kep_ast(6) = kep_ast(6)-0.0111; % Adjust timing for very close flyby

kep_ast(6) = kep_ast(6)-0.01075;

state_ast = cspice_conics(kep_ast, et);
%-------------------------

kep0_sma = kep_ast';
kep0_sma(1) = kep_ast(1)/(1-kep_ast(2));

kepE_sma = kep_eat';
kepE_sma(1) = kep_eat(1)/(1-kep_eat(2));

MOID0 = MOID_SDG_win( kep0_sma([1 2 4 3 5]), kepE_sma([1 2 4 3 5]) );
MOID0_AU = MOID0/cons.AU 

%-------------------------


state_ast_O = HillRot(state_eat, state_ast);
r_ast_O = state_ast_O(1:3);
v_ast_O = state_ast_O(4:6);

figure(1); 
plot( 0/3600/24, norm(r_ast_O), 'ro' );
grid on
hold on
xlabel('dt (days)')
ylabel('d (km)')

% Unperturbed prop from SOI crossing to pericenter
kep_ast_O = cspice_oscelt( state_ast_O, et, cons.GMe );

% Plot coordinates given by hyperbolic elements
tv = (-24*5:1:24*5) *3600 ;
d  = zeros(length(tv),1);
for i=1:length(tv)
    state_ast_t = cspice_conics(kep_ast_O, et - tv(i)); % Why negative?
    d(i) = norm(state_ast_t(1:3)) ;
    
    dr_ast_hyp(i,:) = state_ast_t(1:3);
    
    V   = norm( state_ast_t(4:6) );
    thv(i) = acos( state_ast_t(5)/V );
    phv(i) = atan2( state_ast_t(4), state_ast_t(6) );
    
    

    
end

sc3d = cons.Re;
Dmin_RE = min(d)/sc3d
Dmin_AU = min(d)/cons.AU

figure(1);
plot( tv/3600/24, d )

% Coordinates at perigee
MA_per = 0;
% Coordinates at node crossing
% TA_node = -kep_ast_O(5);
% MA_node = TA_2_MA( TA_node, kep_ast_O(2) );

a_ast_O = kep_ast_O(1)/(1 - kep_ast_O(2));
n_ast_O = sqrt( -cons.GMe/a_ast_O(1)^3 );
dt_per  = -(MA_per - kep_ast_O(6))/n_ast_O;

state_ast_per = cspice_conics( kep_ast_O, et - dt_per );
r_ast_O = state_ast_per(1:3);
v_ast_O = state_ast_per(4:6);

plot( dt_per/3600/24, norm(state_ast_per(1:3)), 'r+' )

%------
F = figure(11); clf;
F.Position = [-1114 6 944 714];

plot3( dr_ast_hyp(:,1)/sc3d, dr_ast_hyp(:,2)/sc3d, dr_ast_hyp(:,3)/sc3d,...
    'LineWidth', 1.5 )
hold on
[EX,EY,EZ] = sphere(20); R = cons.Re;
surf( EX*R/sc3d, EY*R/sc3d, EZ*R/sc3d, ...
    'EdgeColor', 'none', 'FaceAlpha', .2 )
colormap winter
grid on
axis equal
axis([-1 1 -1 1 -1 1]*25*cons.Re/sc3d)
view(-50,20)
xlabel('X_{H} (km) (~ -Sun)')
ylabel('Y_{H} (km) (// V_p)')
zlabel('Z_{H} (km) (// H_p)')
plot3( r_ast_O(1)/sc3d, r_ast_O(2)/sc3d, r_ast_O(3)/sc3d , 'c+')


figure(12);
subplot(2,1,1)
plot( tv/86400, thv ); ylabel('\theta')
grid on; hold on
subplot(2,1,2)
plot( tv/86400, phv ); ylabel('\phi')
grid on; hold on


%% 3.B. Validate call of oscelt with heliocentric elemenets
kep_ast_S = cspice_oscelt( state_ast, et, cons.GMs );
kep_eat_S = cspice_oscelt( state_eat, et, cons.GMs );

d  = zeros(length(tv),1);
for i=1:length(tv)
    state_ast_t = cspice_conics(kep_ast_S, et + tv(i));
    state_eat_t = cspice_conics(kep_eat_S, et + tv(i));
    
    d(i) = norm(state_ast_t(1:3)-state_eat_t(1:3)) ;
    
    state_ast_O = HillRot(state_eat_t, state_ast_t);
    dr_ast_sun(i,:) = state_ast_O(1:3);    
    
end

figure(1);
plot( tv/3600/24, d )
% legend('et','Hyp','dt_per','Hel')

% figure(11);
% plot3( dr_ast_sun(:,1)/sc3d, dr_ast_sun(:,2)/sc3d, dr_ast_sun(:,3)/sc3d ,...
%     'LineWidth', 1.5)
% 
% L = legend('ast-HYP','Earth','ast-HEL');
% L.Position = [0.78031 0.80714 0.19286 0.16548];


%% 4. Opik variables at b-plane
% Coordinates at MTP - Perigee
b     = norm( r_ast_O );
V     = norm( v_ast_O );
theta = acos( v_ast_O(2)/V );
phi   = atan2( v_ast_O(1), v_ast_O(3) );
 
% plotBplane( phi, theta, 5, 11, 'r' );

tp = dt_per;
figure(12); 
subplot(2,1,1); plot( tp/86400, theta, 'r+' );
subplot(2,1,2); plot( tp/86400, phi, 'r+' );

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

v_b_pre = UVec;

cp = cos(phi);   sp = sin(phi);
ct = cos(theta); st = sin(theta);

plotBplane( phi, theta, 20, 11, 'b' );
L.String{5} = 'B-plane';
L.String{6} = 'V_{inf}';

tp = dt_per;

figure(12); 
tp2 = (tv([1 end]));

subplot(2,1,1); 
plot( tp/86400, theta, 'b+' );
plot( tp2/86400,[1 1]*theta,'b--')
legend('\theta','\theta-Per','\theta-MTP')

subplot(2,1,2); 
plot( tp/86400, phi, 'b+' );
plot( tp2/86400,[1 1]*phi,'b--')
legend('\phi','\phi-Per','\phi-MTP')


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
% plot(3*cos(thv), 3*sin(thv),'k--');

grid on
axis equal
caxis([1 kmax])
cb = colorbar;
cb.Label.String = 'k';
axis([-1 1 -1 1]*15)
xlabel('\xi (R_\oplus)');
ylabel('\zeta (R_\oplus)');

plot( xi/sc, zeta/sc, '+r' )

%% 6. Keyhole computation
% Section dependencies: scripts in 'keyholes'

RE_au = cons.Re/DU;

%---- Choose new circles ----
m  = cons.GMe/cons.GMs ;
t0 = 0;
circles = circ;


% Pick flyby from the keyholes and generate ICs
kref = 2;
ic = find( circles(:,1) == kref, 1 );
% ic = ic + 1;
% ic = 142;

% kref = 3;
% ic = 7;

k = circles(ic,1);
h = circles(ic,2);
D = circles(ic,3)/cons.Re;
R = circles(ic,4)/cons.Re;

[kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
    two_keyholes(k, h, D, R, U_nd, theta, phi, m,t0,DU);

% Plot selected keyhole
figure(2)
cc = [1 0 0];
sc = cons.Re/DU;

plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
    'Color',cc);
plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
    'Color',cc);

figure(11)
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];
for i=1:length(kh_up_zeta(:,1))
    kh3d(i,:) = auxR'*[kh_up_xi(i,1); 0; kh_up_zeta(i,1)]/sc ;
end
plot3( kh3d(:,1),kh3d(:,2),kh3d(:,3), 'Color', cc )

auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];
for i=1:length(kh_down_zeta(:,1))
    kh3d(i,:) = auxR'*[kh_down_xi(i,1); 0; kh_down_zeta(i,1)]/sc ;
end
plot3( kh3d(:,1),kh3d(:,2),kh3d(:,3), 'Color', cc )


% 1. Find a point close to xi=0 (arbitrary choice)
% [~,ik] = min( abs(kh_up_xi) );
% 2. Find the first point of the keyhole arc (arbirary as well)
ik = find( ~isnan(kh_up_xi), 1 ); 
xi0   = kh_up_xi(ik);
zeta0 = kh_up_zeta(ik);

% ik = find( ~isnan(kh_down_xi), 1 ); 
% xi0   = kh_down_xi(ik);
% zeta0 = kh_down_zeta(ik);

auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];
r0   = auxR'*[xi0; 0; zeta0];

figure(11);
plot3( r0(1)*DU/sc3d, r0(2)*DU/sc3d, r0(3)*DU/sc3d , 'g+')

% Find post-encounter Opik coordinates
[zeta2,theta1,phi1,xi1,zeta1] = opik_next(U_nd,theta,phi,xi0(1),zeta0(1),t0,h,m);

%--- Validate post-encounter phi and theta
figure(12)
subplot(2,1,1); plot( tp2/86400,[1 1]*theta1,'m--')
subplot(2,1,2); plot( tp2/86400,[1 1]*phi1,'m--')

% Post-encounter Opik to Cartesisan
cp = cos(phi1);   sp = sin(phi1);
ct = cos(theta1); st = sin(theta1);
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];

% Alternative rotation of the relative velocity
hgamma = atan( cons.GMe/(b*U*U) );
hv     = cross( r_ast_O, v_ast_O );
DCM    = PRV_2_DCM( hgamma, hv/norm(hv) ); %%% WRONG hv!?!?!?!?
UVec   = (U/V)*DCM*v_ast_O;
theta1 = acos(UVec(2)/U);
phi1   = atan2(UVec(1),UVec(3));

%--- Validate post-encounter phi and theta
figure(12)
subplot(2,1,1); plot( tp2/86400,[1 1]*theta1,'g--')
subplot(2,1,2); plot( tp2/86400,[1 1]*phi1,'g--')

% plotBplane( phi1, theta1, 20, 11, 'r' );

cp = cos(phi1);   sp = sin(phi1);
ct = cos(theta1); st = sin(theta1);
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];

% Rphi = rotationM(-phi1,2);
% Rth  = rotationM(-theta1,1);
% auxR  = Rth*Rphi ;

r_ast_O1 = DU*auxR'*[xi1; 0; zeta1];
v_ast_O1 = U_nd*(DU/TU)*[st*sp; ct; st*cp];

figure(11)
plot3( r_ast_O1(1)/sc3d, r_ast_O1(2)/sc3d, r_ast_O1(3)/sc3d, 'r+' )

state_ast_post = [r_ast_O1; v_ast_O1];
state_eat = cspice_conics(kep_eat, et + dt_per);

state_ast_post_sun = HillRotInv( state_eat, state_ast_post );
kep_ast_post = cspice_oscelt( state_ast_post_sun, et+dt_per, cons.GMs )

kep_ast_post_O = cspice_oscelt( state_ast_post, et+dt_per, cons.GMe )


% Propagation
et0      = et+dt_per;
kep0     = kep_ast_post;
kep0_sma = kep_ast_post';
kep0_sma(1) = kep0(1)/(1-kep0(2));

kepE = cspice_oscelt( state_eat, et+dt_per, cons.GMs );
kepE_sma = kepE';
kepE_sma(1) = kepE(1)/(1-kepE(2));

MOID0 = MOID_SDG_win( kep0_sma([1 2 4 3 5]), kepE_sma([1 2 4 3 5]) );

state_jup = cspice_spkezr( '5',  et+dt_per, 'ECLIPJ2000', 'NONE', '10' );
kepJ = cspice_oscelt( state_jup, et+dt_per, cons.GMs );
kepJ_sma = kepJ;
kepJ_sma(1) = kepJ(1)/(1-kepJ(2));


% Analytical equations for post-encounter Keplerian elements
U = U_nd*(DU/TU);
sma_post = kepE_sma(1)/( 1 - U_nd*U_nd - 2*U_nd*ct )

% Switch to heliocentric after a few days

dt_SOI = 0*3600*24;

state_ast_post = [r_ast_O1; v_ast_O1];
kep_ast_post_O = cspice_oscelt( state_ast_post, et+dt_per, cons.GMe );

state_ast_SOI  = cspice_conics( kep_ast_post_O, et + dt_per - dt_SOI ); % Note the minus only here
state_eat      = cspice_conics( kep_eat,        et + dt_per + dt_SOI );

% Switch to Heliocentric frame
state_ast_post_sun = HillRotInv( state_eat, state_ast_SOI );
kep_ast_post = cspice_oscelt( state_ast_post_sun, et+dt_per+dt_SOI, cons.GMs );



%% New post-encounter elements generation
% 1. Choose keyhole point in the modified target plane
% 2. Compute post-encounter coordinates for the keyhole
% 3. Transport to perigee, that uniquely defines the hyperbola
% 4. Propagate to outgoing assymptote
% 5. Switch to heliocentric elements
% 6. Find the resonant encounter (Next code section)

% Q: Should be the same as pre-encounter + full rotation (new turn angle)

% 1. Keyhole point in MTP (PRE-ENCOUNTER COORDINATES)
ik = find( ~isnan(kh_up_xi), 1 ); 

xi0   = kh_up_xi(ik);
zeta0 = kh_up_zeta(ik);

cp = cos(phi);   sp = sin(phi);
ct = cos(theta); st = sin(theta);
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];
r_b_pre_kh = DU*auxR'*[xi0; 0; zeta0]; % Verified at plot
%!! v_b_pre defined in the initial flyby

% 2. Post-encounter coordinates for the keyhole point
% Obtained by formulae
[zeta2,theta1,phi1,xi1,zeta1] = opik_next(U_nd,theta,phi,xi0(1),zeta0(1),t0,h,m);

cp = cos(phi1);   sp = sin(phi1);
ct = cos(theta1); st = sin(theta1);
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];

r_b_post_kh = DU*auxR'*[xi1; 0; zeta1];
v_b_post_kh = U_nd*(DU/TU)*[st*sp; ct; st*cp];

% 3. Transport to perigee
b      = sqrt( xi0^2 + zeta0^2 )*DU;
hgamma = atan( cons.GMe/(b*U*U) );

hv      = cross( r_b_pre_kh, v_b_pre );
DCM     = PRV_2_DCM( -hgamma, hv/norm(hv) );

v_per_post_kh = (V/U)*DCM*v_b_post_kh;

theta_p = acos(v_per_post_kh(2)/V);
phi_p   = atan2( v_per_post_kh(1),v_per_post_kh(3) );
cp = cos(phi_p);   sp = sin(phi_p);
ct = cos(theta_p); st = sin(theta_p);
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];

r_per_post_kh = DU*(U/V)*auxR'*[xi1;0;zeta1]; % Not sure why this equation for rescaling

state_ast_post_kh = [r_per_post_kh; v_per_post_kh];
kep_ast_post_O    = cspice_oscelt( state_ast_post_kh, et+dt_per, cons.GMe );

% 4. Propagate to outgoing assymptote
dt_SOI = 29*3600*24;

state_ast_SOI  = cspice_conics( kep_ast_post_O, et + dt_per - dt_SOI ); % Note the minus only here
state_eat      = cspice_conics( kep_eat,        et + dt_per + dt_SOI );

% 5. Switch to Heliocentric frame
state_ast_post_sun = HillRotInv( state_eat, state_ast_SOI );
kep_ast_post = cspice_oscelt( state_ast_post_sun, et+dt_per+dt_SOI, cons.GMs )

sma_post = kep_ast_post(1)/(1 - kep_ast_post(2))

href = circles(ic,2);
sma_teo  = (kref/href)^(2/3) *kepE_sma(1)


%% Alternative rotation: From pre-encounter
hv      = cross( r_b_pre_kh, v_b_pre );
DCM     = PRV_2_DCM( 2*hgamma, hv/norm(hv) );

v_b_post_kh = DCM*v_b_pre;

theta_b = acos(v_b_post_kh(2)/U);
phi_b   = atan2( v_b_post_kh(1),v_b_post_kh(3) );
cp = cos(phi_b);   sp = sin(phi_b);
ct = cos(theta_b); st = sin(theta_b);
auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];

r_b_post_kh = DU*(U/V)*auxR'*[xi1;0;zeta1];

state_ast_post_kh = [r_b_post_kh; v_b_post_kh];
kep_ast_post_O    = cspice_oscelt( state_ast_post_kh, et+dt_per, cons.GMe );

% 4. Propagate to outgoing assymptote
dt_SOI = 2*3600*24;

state_ast_SOI  = cspice_conics( kep_ast_post_O, et + dt_per - dt_SOI ); % Note the minus only here
state_eat      = cspice_conics( kep_eat,        et + dt_per + dt_SOI );

% 5. Switch to Heliocentric frame
state_ast_post_sun = HillRotInv( state_eat, state_ast_SOI );
kep_ast_post = cspice_oscelt( state_ast_post_sun, et+dt_per+dt_SOI, cons.GMs )

sma_post = kep_ast_post(1)/(1 - kep_ast_post(2))

href = circles(ic,2);
sma_teo  = (kref/href)^(2/3) *kepE_sma(1)


%% Propagations near encounter for validation
% - 2BP with conics
% - 3BP integration
% - Opik

state_pre_kh = [r_b_pre_kh; v_b_pre];
kep_pre_kh   = cspice_oscelt( state_pre_kh, et+dt_per, cons.GMe )

% 2BP conics
tv = ( 0:0.1:15 )*86400;
for i=1:length(tv)
    state_ast_O_t = cspice_conics( kep_pre_kh, et+dt_per - tv(i) );
    state_eat     = cspice_conics( kep_eat,    et+dt_per + tv(i) );

    r_ast_O_2bp(i,:) = state_ast_O_t(1:3);
    
    state_ast_post_sun = HillRotInv( state_eat, state_ast_O_t );
    kep_ast_post_t     = cspice_oscelt( state_ast_post_sun, et+dt_per+tv(i), cons.GMs );
    kep_post_2bp(i,:) = kep_ast_post_t;
    kep_post_2bp(i,1) = kep_ast_post_t(1)/(1-kep_ast_post_t(2));
    
end
figure(11)
plot3( r_ast_O_2bp(:,1)/sc3d, r_ast_O_2bp(:,2)/sc3d, r_ast_O_2bp(:,3)/sc3d,'m')

plot_oe_evolution( 'oe', tv, kep_post_2bp, 15, 'm', tv([1 end])/86400, 86400 )


% 3BP integration
state_eat = cspice_conics( kep_eat, et+dt_per );
X0        = HillRotInv( state_eat, state_pre_kh );

cons_ode.GMcentral = cons.GMs;
cons_ode.GMthird   = cons.GMe;
cons_ode.GMfourth = 0;
cons_ode.OEp   = kepE;
cons_ode.OEp2  = kepJ;
cons_ode.JD0_p = et0;

tol = 1e-13;
options=odeset('RelTol',tol,'AbsTol',ones(1,6)*tol);
[t,X]=ode113(@(t,X) multibodies_3rd4th_oe(t,X,cons_ode),tv,X0,options);
for i=1:length(tv)
    kep0_4bp = cspice_oscelt( X(i,:)', et0+tv(i), cons.GMs );
    kep_post_4bp(i,:) = kep0_4bp;
    kep_post_4bp(i,1) = kep0_4bp(1)/(1-kep0_4bp(2));
    
    % Compute distance to Earth
    xe   = cspice_conics(kepE, et0+tv(i) );
    r_ast_O_4bp(i,:) = HillRot(xe, X(i,:)');
    
    d(i) = norm(xe(1:3) - X(i,1:3)'); % From 4bp integration
end

figure(11)
plot3( r_ast_O_4bp(:,1)/sc3d, r_ast_O_4bp(:,2)/sc3d, r_ast_O_4bp(:,3)/sc3d,'g')

plot_oe_evolution( 'oe', tv, kep_post_4bp, 15, 'g', tv([1 end])/86400, 86400 )


%% Propagation until the next encounter
% Two body distance propagation
tv = (0:.01:50) *cons.yr;
t0 = et+dt_per+dt_SOI ;
for i=1:length(tv)
    
    xa = cspice_conics( kep_ast_post, t0 + tv(i) );
    xe = cspice_conics( kep_eat,      t0 + tv(i) );
    
    d2BP(i) = norm(xe(1:3) - xa(1:3));
end

% Is the next encounter happening?
F = figure(6);
clf

xsc = cons.yr;
ysc = cons.Re;

plot(tv/xsc, d2BP/ysc)
hold on
grid on
xlabel('time (yr)')
ylabel('d (R_\oplus)')

legend('2BP')




