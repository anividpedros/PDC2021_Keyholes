% MODIFIED KEYHOLES: PDC 2021
% 7th IAA Planetary Defense Conference 2021
clc
clear
close all
format shortG

filename = 'main_sect45_modified_keyholes_circularEarth.m';
filepath = matlab.desktop.editor.getActiveFilename;
currPath = filepath(1:(end-length(filename)));
cd(currPath)
addpath(genpath(pwd))

%% ===== Description ===== 
% Compute modified keyholes by time-varying MOID
% Extension: Generate distance(t) plots given post-encounter coordinates
% =======================

%% 1. Read heliocentric elements
dir_local_de431 = 'C:\Users\Oscar\Documents\Spice-Kernels\';
% dir_local_de431 = '/Users/anivid/ExampleMICE/kernels/spk/';

% addpath(genpath(mice_local_path))
cspice_furnsh( 'SPICEfiles/naif0012.tls.pc' )
cspice_furnsh( 'SPICEfiles/gm_de431.tpc' )
cspice_furnsh( 'SPICEfiles/pck00010.tpc' )
cspice_furnsh( [dir_local_de431 'de431_part-1.bsp'] )
cspice_furnsh( [dir_local_de431 'de431_part-2.bsp'] )

% 2021 PDC
cspice_furnsh( 'SPICEfiles/2021_PDC-s11-merged-DE431.bsp' )
ast_id= '-937014';
epoch = '2021 October 20 TDB';
et = cspice_str2et( epoch ) + 24*.715*3600;

cons.AU  = cspice_convrt(1, 'AU', 'KM');
cons.GMs = cspice_bodvrd( 'SUN', 'GM', 10 );
% cons.GMe = cspice_bodvrd( '399', 'GM', 1 );
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
cons.yr = 2*pi/sqrt(cons.GMs/kep_eat(1)^3);

% Modify NEO to try to make dCA smaller ----
kep_ast(3) = kep_ast(3)+.0;     % Increase inclination
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

b_nd = b/DU;
c_nd = (cons.GMe/cons.GMs)/U_nd^2;

b2 = b_nd*b_nd;
c2 = c_nd*c_nd;
U2 = U_nd*U_nd;


a_pre  = 1/(1-U2-2*U_nd*ct);

kmax = 20;
hmax = 50;

i = 0;
circ = [];

for k=1:kmax
    for h=1:hmax
        
        i = i+1 ;
        a0p(i) = (k/h)^(2/3);
        if sum(find(a0p(i) == a0p(1:i-1))); continue; end
        ct0p = ( 1-U2-ap/a0p(i) )/2/U_nd ;
        
        if abs(ct0p) > 1; continue; end
        st0p = sqrt(1 - ct0p^2);
        
        % Cirle center location and radius
        D = (c_nd*st)/(ct0p - ct);
        R = abs( c_nd*st0p/(ct0p - ct) );
        
        circ = [circ; k h D*DU R*DU];
        
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
subplot(1,2,1)
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
        'Color',cc,'LineStyle','--');
    pl = plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
         'Color',cc,'LineStyle','--');
     
    % Register keyholes with solutions
%     arcexist = sum(~isnan(kh_up_xi)) + sum(~isnan(kh_down_xi));
    R1 = sqrt(sum(kh_down_xi.^2 + kh_down_zeta.^2,2));
    R2 = sqrt(sum(kh_up_xi.^2 + kh_up_zeta.^2,2));
    arcexist = sum( (R1-3*RE_au*focus_factor)>0 ) + sum( (R2-3*RE_au*focus_factor)>0 );
    
%     arcexist = sum( kh_up_zeta > 2*RE_au*focus_factor );
%     arcexist = sum( kh_down_zeta < -5*RE_au );
%     arcexist = sum( kh_up_zeta > 3.6*RE_au );
    
    if arcexist
        kh_good = [kh_good; i];
        fprintf('Keyhole num %g exists\n',i)
    end
    
    if i==1
        h_list(1) = pl(1);
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




%% MOIDIFIED Keyholes Computation
% Section dependencies: scripts in 'keyholes'
RE_au = cons.Re/DU;

m  = cons.GMe/cons.GMs ;
circles = circ;

F = figure(3);
hold on

sc = cons.Re/DU;
nr = size(circles,1);
kh_good = [];

% longitude for MOID computation inside
longp = atan2(state_eat(2),state_eat(1));

% Initial condition
kepE_sma = kep_eat';
kepE_sma(1) = kep_eat(1)/(1-kep_eat(2));

et0 = t0 + dt_per;
eti = et0 + 30*86400 ; % Initial ephemeris time for integration

subplot(1,2,1)

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

% Numerical Model Inputs
cons_ode.t0        = eti ;
cons_ode.GM_vec    = GMvec ; % 1st element is Sun
cons_ode.IC_planet = kep_planet ; % 1 column per planet


for i=1:nr
    
    % New circles
    k = circles(i,1);
    h = circles(i,2);
    D = circles(i,3)/cons.Re;    
    R = circles(i,4)/cons.Re;    
    
%     try
    [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
        two_keyholes_dxi_sec(k, h, D, R, U_nd, theta, phi, m,0,DU,longp,ap,cons,kepE_sma,cons_sec);
%     [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
%           two_keyholes(k, h, D, R, U_nd, theta, phi, m,0,DU);
%     catch        
%     end

%     [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
%         two_keyholes_dxi_num(k, h, D, R, U_nd, theta, phi, m,0,DU,longp,ap,cons,kepE_sma,cons_ode);

    if sum(~isnan(kh_up_xi(:)))
        fprintf('Keyhole %g found!\n',i)
    end
    
    subplot(1,2,1)
    cc = co(k,:);    
%     cc = [0 1 1];
%     cc = [1 0 0];
    plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
        'Color',cc,'LineWidth',1);
    pl = plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
         'Color',cc,'LineWidth',1);
    
%     drawnow 
    % Register keyholes with solutions
%     arcexist = sum(~isnan(kh_up_xi)) + sum(~isnan(kh_down_xi));
    R1 = sqrt(sum(kh_down_xi.^2 + kh_down_zeta.^2,2));
    R2 = sqrt(sum(kh_up_xi.^2 + kh_up_zeta.^2,2));
    arcexist = sum( (R1-RE_au*focus_factor)>0 ) + sum( (R2-RE_au*focus_factor)>0 );
    
    arcexist = sum( kh_up_zeta > 2*RE_au*focus_factor );
%     arcexist = sum( kh_down_zeta < -1.6*RE_au );
    
    if arcexist
        kh_good = [kh_good; i];
        fprintf('Keyhole num %g exists\n',i)
    end
    
    if length(kh_good) == 1
        h_list(2) = pl(1);
    end
    
end
colormap(co);
fill(RE_focussed*cos(thv), RE_focussed*sin(thv),'white');
plot(RE_focussed*cos(thv), RE_focussed*sin(thv),'k');
plot(cos(thv), sin(thv),'k--');

legend(h_list,{'post-enc','sec-LL'})

%%
grid on
axis equal
caxis([1 20])
cb = colorbar;
cb.Label.String = 'k';
axis([-1 1 -1 1]*5)
xlabel('\xi (R_\oplus)');
ylabel('\zeta (R_\oplus)');
