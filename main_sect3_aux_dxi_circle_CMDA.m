% SOLUTION OF A PLANETARY ENCOUNTER: PDC 2021
% 7th IAA Planetary Defense Conference 2021
clc
clear
close all
format shortG

filename = 'main_sect3_aux_dxi_circle_CMDA.m';
filepath = matlab.desktop.editor.getActiveFilename;
currPath = filepath(1:(end-length(filename)));
cd(currPath)
addpath(genpath(pwd))

%% ===== Description ===== 
% Original script: main_sect3_multimethod_dist_circular.m
% Show the variation of the \Delta MOID for every point of the resonant
% circle
% =======================

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
% dir_local_de431 = 'C:\Users\Oscar\Documents\Spice-Kernels\';
% % dir_local_de431 = '/Users/anivid/ExampleMICE/kernels/spk/';
% 
% % addpath(genpath(mice_local_path))
% cspice_furnsh( 'SPICEfiles/naif0012.tls.pc' )
% cspice_furnsh( 'SPICEfiles/gm_de431.tpc' )
% cspice_furnsh( 'SPICEfiles/pck00010.tpc' )
% cspice_furnsh( [dir_local_de431 'de431_part-1.bsp'] )
% cspice_furnsh( [dir_local_de431 'de431_part-2.bsp'] )
% 
% % 2021 PDC
% cspice_furnsh( 'SPICEfiles/2021_PDC-s11-merged-DE431.bsp' )
% ast_id= '-937014';
% epoch = '2021 October 20 TDB';
% et = cspice_str2et( epoch ) + 24*.715*3600;
% 
% % % Bennu
% % ast_id= '-101955';
% % epoch = '2060-Sep-22 00:36';
% % 
% % % 2017 PDC
% % cspice_furnsh( 'SPICEfiles/2017_PDC-merged-DE431.bsp' )
% % ast_id= '-937001';
% % epoch = '2027 July 20 15:00 TDB';
% % et = cspice_str2et( epoch ) + 24*.715*3600;
% % 
% % % Apophis
% % epoch = '2029-Apr-12 21:46:00.0000 TDB';
% % et = cspice_str2et( epoch ) ;
% 
% cons.AU  = cspice_convrt(1, 'AU', 'KM');
% cons.GMs = cspice_bodvrd( 'SUN', 'GM', 10 );
% % cons.GMe = cspice_bodvrd( '399', 'GM', 1 );
% cons.GMe = 398600.43543609593;
% cons.Re  = 6378.140;
% 
% cons.yr  = 365.25 * 24 * 3600 ;
% cons.Day = 3600*24; 
% 
% state_eat = cspice_spkezr( '399', et, 'ECLIPJ2000', 'NONE', '10' );
% kep_eat = cspice_oscelt( state_eat, et, cons.GMs );
% sma_eat = kep_eat(1)/(1-kep_eat(2));
% 
% state_ast = cspice_spkezr( ast_id, et, 'ECLIPJ2000', 'NONE', '10' );
% kep_ast = cspice_oscelt( state_ast, et, cons.GMs );
% sma_ast = kep_ast(1)/(1-kep_ast(2));
% 
% % kep_ast = [1.109226568768932E+08 1.952935182366752E-01 3.415493133874123E+00 2.037874906198153E+02 1.266459636185671E+02 2.518558078117354E+02]';
% % kep_ast(3:6) = kep_ast(3:6)*pi/180;
% % kep_ast(7:8) = [et; cons.GMs];
% 
% 
% 
% %% 2. Modify elements for encounter with simpler Earth model
% % Make Earth Circular-Ecliptic -----
% kep_eat(2:3) = 0;
% cons.yr = 2*pi/sqrt(cons.GMs/kep_eat(1)^3);
% 
% % Modify NEO to try to make dCA smaller ----
% kep_ast(3) = kep_ast(3)+.0;     % Increase inclination
% kep_ast(5) = kep_ast(5)+.06;    % Adjust arg of perihelion to low MOID
% kep_ast(6) = kep_ast(6)-0.0111; % Adjust timing for very close flyby
% 
% % Generate states moments before the flyby
% dt = -4*24 * 3600;
% t0 = et + dt;
% 
% state_ast = cspice_conics(kep_ast, t0);
% state_eat = cspice_conics(kep_eat, t0);


%% 1. and 2. 

%========= 2006 MB14 SETUP ===========
% load('Horizons-2006MB14.mat')
% 
% % Exact encounter date
% epoch  = '1985-Jun-28 00:40 TDB';
% et_enc = cspice_str2et( epoch );
% 
% % Defining a "nominal" time
% id_nom = 74;
% date_cal = char(CalendarDateTDB(id_nom));
% epoch = [date_cal(7:17) ' TDB'];

%========= APOPHIS SETUP ===========
load('horizons_apophis.mat')
% id = 74;
id_nom = 72;
date_cal = char(CalendarDateTDB(id_nom));
epoch = [date_cal(7:17) ' TDB'];

%===================================
% Load constants
dir_local_de431 = 'C:\Users\Oscar\Documents\Spice-Kernels\';
% dir_local_de431 = '/Users/anivid/ExampleMICE/kernels/spk/';

% addpath(genpath(mice_local_path))
cspice_furnsh( 'SPICEfiles/naif0012.tls.pc' )
cspice_furnsh( 'SPICEfiles/gm_de431.tpc' )
cspice_furnsh( 'SPICEfiles/pck00010.tpc' )
cspice_furnsh( [dir_local_de431 'de431_part-1.bsp'] )
cspice_furnsh( [dir_local_de431 'de431_part-2.bsp'] )

cons.AU  = cspice_convrt(1, 'AU', 'KM');
cons.GMs = cspice_bodvrd( 'SUN', 'GM', 10 );
% cons.GMe = cspice_bodvrd( '399', 'GM', 1 );
cons.GMe = 398600.43543609593;
cons.Re  = 6378.140;

cons.yr  = 365.25 * 24 * 3600 ;
cons.Day = 3600*24; 
%===================================

et = cspice_str2et( epoch );
t0 = et;

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
state_ast = cspice_conics(kep_ast, t0);

kep_opik_post = kep_ast;
% kep_opik_post(1) = kep_opik_post(1)*(1-kep_opik_post(2))*DU ; 
kep_opik_post(7:8) = [t0; cons.GMs];




%% 3. Switch to geocentric and find time periapses === [Section for Apophis]
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
subplot(1,2,1)

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

% plot( xi/sc, zeta/sc, '+k' )


%% General Keyholes computation
% Section dependencies: scripts in 'keyholes'
RE_au = cons.Re/DU;

m  = cons.GMe/cons.GMs ;
circles = circ;
% 
% F = figure(3);
% hold on
% 
% sc = cons.Re/DU;
% nr = size(circles,1);
% kh_good = [];
% for i=1:nr
%     
%     % New circles
%     k = circles(i,1);
%     h = circles(i,2);
%     D = circles(i,3)/cons.Re;    
%     R = circles(i,4)/cons.Re;    
%     
%     [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
%         two_keyholes(k, h, D, R, U_nd, theta, phi, m,0,DU);
%     
%     cc = co(k,:);    
%     plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
%         'Color',cc);
%     plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
%          'Color',cc);
%      
%     % Register keyholes with solutions
% %     arcexist = sum(~isnan(kh_up_xi)) + sum(~isnan(kh_down_xi));
%     R1 = sqrt(sum(kh_down_xi.^2 + kh_down_zeta.^2,2));
%     R2 = sqrt(sum(kh_up_xi.^2 + kh_up_zeta.^2,2));
%     arcexist = sum( (R1-3*RE_au*focus_factor)>0 ) + sum( (R2-3*RE_au*focus_factor)>0 );
%     
% %     arcexist = sum( kh_up_zeta > 2*RE_au*focus_factor );
% %     arcexist = sum( kh_down_zeta < -5*RE_au );
% %     arcexist = sum( kh_up_zeta > 3.6*RE_au );
%     
%     if arcexist
%         kh_good = [kh_good; i];
%         fprintf('Keyhole num %g exists\n',i)
%     end
%     
%     
% end
% colormap(co);
% fill(RE_focussed*cos(thv), RE_focussed*sin(thv),'white');
% plot(RE_focussed*cos(thv), RE_focussed*sin(thv),'k');
% plot(cos(thv), sin(thv),'k--');
% 
% grid on
% axis equal
% caxis([1 20])
% cb = colorbar;
% cb.Label.String = 'k';
% axis([-1 1 -1 1]*10)
% xlabel('\xi (R_\oplus)');
% ylabel('\zeta (R_\oplus)');
% 

%% Keyhole selection
ic = 22;

% Pick flyby from the keyholes and generate ICs
F = figure(ic+1); clf;

for i=ic%:nr
    
    k = circ(i,1);
    D = circ(i,3);    
    R = circ(i,4);
    
    xi_circ   = R*cos(thv);
    zeta_circ = D + R*sin(thv);
    
    cc = co(k,:);
    plot( xi_circ/sc, zeta_circ/sc, 'r' )
    hold on
    
end

% Solving for \deltaMOID of all points and plot
subplot(1,2,1)

% Redraw circles ----------------
for i=1:nr
    
    k = circ(i,1);
    D = circ(i,3)/cons.Re;    
    R = circ(i,4)/cons.Re;
    
    xi_circ   = R*cos(thv);
    zeta_circ = D + R*sin(thv);
    
    cc = co(k,:);
    plot( xi_circ, zeta_circ, 'Color', cc )
    hold on
    
end
colormap(co);
fill(RE_focussed*cos(thv), RE_focussed*sin(thv),'white');
plot(RE_focussed*cos(thv), RE_focussed*sin(thv),'k');
plot(cos(thv), sin(thv),'k--');
% plot(3*cos(thv), 3*sin(thv),'k--');
%-----------------------------

% ic = 48;

k = circles(ic,1);
h = circles(ic,2);
D = circles(ic,3)/cons.Re;
R = circles(ic,4)/cons.Re;

kepE_sma = kep_eat';
kepE_sma(1) = kep_eat(1)/(1-kep_eat(2));
et0 = t0 + dt_per;
eti = et0 + 30*86400 ; % Initial ephemeris time for integration

% Secular Propagation
kep_planet = NaN(8,8);
GMvec = cons.GMs;
for i = 2:9
    GMvec(i)          = cspice_bodvrd( num2str(i-1), 'GM', 1);
    state_planet      = cspice_spkezr( num2str(i-1),  eti, 'ECLIPJ2000', 'NONE', '10' );
    kep_planet(:,i-1) = cspice_oscelt( state_planet, eti, cons.GMs );
end
kep_planet(:,3) = kep_eat ;
GMvec(3) = cons.GMe;

% Secular Model: Lagrange-Laplace
kepJ_sma = kep_planet(:,5);
kepJ_sma(1) = kepJ_sma(1)/(1-kepJ_sma(2));

cons_sec.OEp = kepJ_sma';
cons_sec.GMp = GMvec(6);
cons_sec.GMs = GMvec(1);


[kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta,dx_down,dx_up] = ...
    two_keyholes_dxi_sec(k, h, D, R, U_nd, theta, phi, m,0,DU,...
    longp,ap,cons,kepE_sma,cons_sec);

% Plot the secular variation of MOID
axis([-20 20 -20 20]);

nkh = 100;
[xi_up,zeta_up] = res_circle(linspace(0,pi,nkh),D,R);
[xi_down,zeta_down] = res_circle(linspace(-pi,0,nkh),D,R);

% scatter( xi_up, zeta_up, 5, dx_up, 'filled' )

Np = 100;
xv = 0:1/(Np-1):1;
cm = jet(Np);

scc = [min([dx_up; dx_down]) max([dx_up; dx_down])];
scc = [-5 5]*cons.Re/cons.AU;

dx_c = (dx_up - scc(1))/(scc(2)-scc(1));
dx_cm = interp1(xv,cm,dx_c);
scatter( xi_up, zeta_up, 10, dx_cm, 'filled' )
hold on

dx_c = (dx_down - scc(1))/(scc(2)-scc(1));
dx_cm = interp1(xv,cm,dx_c);
scatter( xi_down, zeta_down, 10, dx_cm, 'filled' )

colormap(cm);

caxis(scc*cons.AU/cons.Re)

axis equal
axis([-1 1 -1 1]*20)
cb = colorbar;
cb.Label.String = '\Delta \xi (R_\oplus)';
cb.Position(1) = 0.45;
caxis(scc*cons.AU/cons.Re)

txt = text(-1,12,['k=' num2str(k)]);

F.Units = 'inches';
F.Position(3:4) = [6 3];

title('\Delta\xi - Secular')

%===============================
% Numerical Model Inputs
subplot(1,2,2)

% Redraw circles ----------------
for i=1:nr
    
    k = circ(i,1);
    D = circ(i,3)/cons.Re;    
    R = circ(i,4)/cons.Re;
    
    xi_circ   = R*cos(thv);
    zeta_circ = D + R*sin(thv);
    
    cc = co(k,:);
    plot( xi_circ, zeta_circ, 'Color', cc )
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
%--------------------------------
%

% Solution of the Keyhole
% ic = 29;

k = circles(ic,1);
h = circles(ic,2);
D = circles(ic,3)/cons.Re;
R = circles(ic,4)/cons.Re;

cons_ode.t0        = eti ;
cons_ode.GM_vec    = GMvec ; % 1st element is Sun
cons_ode.IC_planet = kep_planet ; % 1 column per planet

[kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta,dx_down,dx_up] = ...
        two_keyholes_dxi_num(k, h, D, R, U_nd, theta, phi, m,0,DU,longp,ap,cons,kepE_sma,cons_ode);

nkh = 10;
[xi_up,zeta_up] = res_circle(linspace(0,pi,nkh),D,R);
[xi_down,zeta_down] = res_circle(linspace(-pi,0,nkh),D,R);

    
Np = 100;
xv = 0:1/(Np-1):1;
cm = jet(Np);

% figure(3);
scc = [min([dx_up; dx_down]) max([dx_up; dx_down])];
% scc(2) = 3*cons.Re/cons.AU;
% scc = [-20 20]*cons.Re/cons.AU;
scc = [-5 5]*cons.Re/cons.AU;


dx_c = (dx_up - scc(1))/(scc(2)-scc(1));
dx_c(dx_c > 1) = 1;
dx_cm = interp1(xv,cm,dx_c);
scatter( xi_up, zeta_up, 20, dx_cm, 'filled' )
hold on

dx_c = (dx_down - scc(1))/(scc(2)-scc(1));
dx_cm = interp1(xv,cm,dx_c);
scatter( xi_down, zeta_down, 20, dx_cm, 'filled' )

colormap(cm)
caxis(scc*cons.AU/cons.Re)
cb.Label.String = '\Delta \xi (R_\oplus)';
cb.Position(1) = 0.90;
txt = text(-1,-12,['k=' num2str(k)]);
% axis([-10 10 -15 5]);
axis([-1 1 -1 1]*20);


title('\Delta\xi - n3BPert')


%% Auxiliar Functions

% K2S(kepE_sma)
% K2S(kep_nbp(1,:))

% MOID1 = MOID_ORCCA_win( K2S(kepE_sma), K2S(kep_nbp(i,:)) )
% MOID2 = ComputeMOID( K2S(kep_nbp(i,:)), K2S(kepE_sma)  )    
% MOID3 = MOID_SDG_win( kep_nbp(i,[1 2 4 3 5]), kepE_sma([1 2 4 3 5]) )

% Keplerian elements into structure for MOID fxn
function A = K2S(OE,AU) 
A.sma   = OE(1)/AU;
A.e     = OE(2);
A.i     = OE(3);
A.Omega = OE(4);
A.argp  = OE(5);
end
