% SOLUTION OF A PLANETARY ENCOUNTER: PDC 2021
% 7th IAA Planetary Defense Conference 2021
clc
clear
close all
format shortG

filename = 'main_sect12_BplaneCircles.m';
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


%% Compute Circles

% xi_nd   = xi/DU;
% zeta_nd = zeta/DU;
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

plot( xi/sc, zeta/sc, '+k' )


%% General Keyholes computation
% Section dependencies: scripts in 'keyholes'
RE_au = cons.Re/DU;

m  = cons.GMe/cons.GMs ;
circles = circ;

F = figure(3);
hold on

sc = cons.Re/DU;
nr = size(circles,1);
kh_good = [];
for i=1:nr
    
    % New circles
    k = circles(i,1);
    h = circles(i,2);
    D = circles(i,3)/cons.Re;    
    R = circles(i,4)/cons.Re;    
    
    [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
        two_keyholes(k, h, D, R, U_nd, theta, phi, m,0,DU);
    
    cc = co(k,:);    
    plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
        'Color',cc);
    plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
         'Color',cc);
     
    % Register keyholes with solutions
%     arcexist = sum(~isnan(kh_up_xi)) + sum(~isnan(kh_down_xi));
    R1 = sqrt(sum(kh_down_xi.^2 + kh_down_zeta.^2,2));
    R2 = sqrt(sum(kh_up_xi.^2 + kh_up_zeta.^2,2));
    arcexist = sum( (R1-RE_au*focus_factor)>0 ) + sum( (R2-RE_au*focus_factor)>0 );
    if arcexist
        kh_good = [kh_good; i];
        fprintf('Keyhole num %g exists\n',i)
    end
    
    
end
colormap(co);
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


%% Keyhole selection
% Pick flyby from the keyholes and generate ICs

% Manually select k
kref = 15;
ic = find( circles(:,1) == kref, 1 ) + 3;
% Max radius from keyholes found
[~,id] = max(abs(circles(kh_good,:)));
ic = kh_good(id(3));
% Finding circles with large radius list
[~,id] = sort(abs(circles(kh_good,3)),'descend');
ic = kh_good(49);
% Final decision
ic = 124;
% ic = 32; % kh_good(7);

k = circles(ic,1);
h = circles(ic,2);
D = circles(ic,3)/cons.Re;
R = circles(ic,4)/cons.Re;

[kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
    two_keyholes(k, h, D, R, U_nd, theta, phi, m,0,DU);

% Plot selected keyhole
figure(2)
cc = [1 0 0];
sc = cons.Re/DU;

plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
    'Color',cc);
plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
    'Color',cc);

% 1. Find a point close to xi=0 (arbitrary choice)
[~,ik] = min( abs(kh_up_xi) );
xi0   = kh_up_xi(ik(1));
zeta0 = kh_up_zeta(ik(1));

% 2. Find the first point of the keyhole arc (arbirary as well)
%--- Select depending on the arch being up or down
% ik = find( ~isnan(kh_up_xi), 1 ); 
% xi0   = kh_up_xi(ik);
% zeta0 = kh_up_zeta(ik);

% ik = find( ~isnan(kh_down_xi), 1 ); 
% xi0   = kh_down_xi(ik);
% zeta0 = kh_down_zeta(ik);

auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];
r0   = auxR'*[xi0; 0; zeta0];

plot(xi0/sc, zeta0/sc,'rd','MarkerSize',8)

%% 6. Solving the encounter with Opik formulae
h = 0; % Number of revolutions until next encounter: only used for zeta2

% Using keyhole point selected: xi0, zeta0
[~,theta1,phi1,xi1,zeta1,r_b_post,v_b_post] = opik_next(U_nd,theta,phi,xi0,zeta0,t0,h,m);


% %% 6.1. TP to MTP
% % Rescaling b-plane coordinates: Angular momentum conservation
% xi1_p   = (U/V)*xi1 ;
% zeta1_p = (U/V)*zeta1 ;
% 
% 
% b    = (V/U)*b_p  ;
% xi   = (V/U)*xi_p ;
% zeta = (V/U)*zeta_p ;
%  
% % Rotation of velocity vector
% hv     = cross( r_ast_O, v_ast_O );
% hv     = cross( r_b_post, v_b_post );
% 
% DCM    = PRV_2_DCM( -hgamma, hv/norm(hv) );
% UVec   = (U/V)*DCM*v_ast_O;
% theta  = acos(UVec(2)/U);
% phi    = atan2(UVec(1),UVec(3));



%% 7. Heliocentric orbit elements from post-encounter coordinates (Opik formulae)
kep_opik_post = opik_bplane_2_oe( theta1,phi1,zeta1,xi1,U_nd,phi,longp,ap )';

kep_opik_post(1) = kep_opik_post(1)*(1-kep_opik_post(2))*DU ; 
% 'opik_bplane_2_oe.m' function returns sma in first element

kep_opik_post(6) = TA_2_MA(kep_opik_post(6),kep_opik_post(2));
% 'opik_bplane_2_oe.m' function returns true anomaly 6th element

kep_opik_post(7:8) = [t0 + dt_per;
                        cons.GMs];
% Include GM of the central body and epoch

% kep_opik_post(1) = a_post_theory*(1-kep_opik_post(2));


%% 8. Plotting: distance over time with heliocentric elements
% Compute distance to Earth
tv = (0.:.001:20) *cons.yr;
et0 = t0 + dt_per;

clear d
for i=1:length(tv)
    xa   = cspice_conics(kep_opik_post, et0 + tv(i) );
    xe   = cspice_conics(kep_eat,       et0 + tv(i) );
    d(i) = norm(xe(1:3) - xa(1:3)); 
end

% Is the next encounter happening?
F = figure(6);
clf
xsc = cons.yr;
ysc = DU; %cons.Re;

plot(tv/xsc, d/ysc)
grid on
xlabel('t (yr)')
ylabel('d (au)')

%% 9. Validation of conversions
% - Does the relationship of the sma's hold?
% - Does the initial position hold? d(0)=b?

k = circles(ic,1);
h = circles(ic,2);

a_post_theory = DU*(k/h)^(2/3)
a_post_opik   = kep_opik_post(1)/(1-kep_opik_post(2))



%% ===== PLOTTING AUX: Trying to create B-plane plot ========
%------- Plot coordinates given by hyperbolic elements
% This part of the section is from older code
kep_ast_O     = cspice_oscelt( state_ast_per, t0 + dt_per, cons.GMe );
r_ast_O = state_ast_per(1:3);
tv = (-24*5:1:24*5) *3600 ;
d  = zeros(length(tv),1);
for i=1:length(tv)
    state_ast_t = cspice_conics(kep_ast_O, t0 + dt_per + tv(i) ); 
    
    d(i) = norm(state_ast_t(1:3)) ;
    dr_ast_hyp(i,:) = state_ast_t(1:3);   
end

F = figure(11); clf;
F.Position = [843 160 944 714];

sc3d = cons.Re;
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

plotBplane( phi, theta, 20, 11, 'b' );

R = 20;
quiver3(r_ast_O(1)/sc3d, r_ast_O(2)/sc3d, r_ast_O(3)/sc3d,...
    R/2*st*sp, R/2*ct, R/2*st*cp,'Color','b',...
        'LineWidth',2,'MaxHeadSize',2)

%------- Plot Bplane
% This part of the section is from older code
plotBplane( phi1, theta1, 20, 11, 'r' );

quiver3(r_ast_O(1)/sc3d, r_ast_O(2)/sc3d, r_ast_O(3)/sc3d,...
    R/2*st*sp, R/2*ct, R/2*st*cp,'Color','b',...
        'LineWidth',2,'MaxHeadSize',2)
