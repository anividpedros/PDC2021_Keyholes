% SOLUTION OF A PLANETARY ENCOUNTER: PDC 2021
% 7th IAA Planetary Defense Conference 2021
clc
clear
% close all
format shortG

filename = 'main_sect12_BplaneCircles_circularEarth_trasteando.m';
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
dir_local_de431 = 'C:\Users\Oscar\Documents\Spice-Kernels\';
% dir_local_de431 = '2nd Author write here your directory';

cspice_furnsh( 'SPICEfiles/naif0012.tls.pc' )
cspice_furnsh( 'SPICEfiles/gm_de431.tpc' )
cspice_furnsh( 'SPICEfiles/pck00010.tpc' )
cspice_furnsh( [dir_local_de431 'de431_part-1.bsp'] )
cspice_furnsh( [dir_local_de431 'de431_part-2.bsp'] )

% 2021 PDC
cspice_furnsh( '2021_PDC-s11-merged-DE431.bsp' )
ast_id= '-937014';
epoch = '2021 October 10 TDB';

et = cspice_str2et( epoch ) + 24*.715*3600;
% et = cspice_str2et( epoch ) - 24* 15;

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
dt = 0;
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
dt_per  = (MA_per - kep_ast_O(6))/n_ast_O; dt_per_day = dt_per/86400


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
b    = (V/U)*b_p ; b_RE = b/cons.Re
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
F  = figure(5);
% subplot(1,2,1);

nr = size(circ,1);
co = winter(22);
thv= 0:(2*pi/99):2*pi ;
sc = cons.Re;

focus_factor = sqrt(1 + 2*cons.GMe/(cons.Re*U^2));
RE_focussed  = focus_factor;

for i=1:nr
    
    k = circ(i,1);
    h = circ(i,2);
    D = circ(i,3);    
    R = circ(i,4);
    
    xi_circ   = R*cos(thv);
    zeta_circ = D + R*sin(thv);
    
    apst = (k/h)^(2/3);
    
    if apst >= 1
        zpl = 1;
    else
        zpl = -1;
    end
    
    cc = co(k,:);
%     plot( xi_circ/sc, zeta_circ/sc, 'Color', cc )
    scatter( xi_circ/sc, zeta_circ/sc, 5, zpl*ones(1,100) ,'filled' ); 
    hold on
    
end
colormap(co);
fill(RE_focussed*cos(thv), RE_focussed*sin(thv),'white');
plot(RE_focussed*cos(thv), RE_focussed*sin(thv),'k');
plot(cos(thv), sin(thv),'k--');
% plot(3*cos(thv), 3*sin(thv),'k--');

grid on
axis equal
% caxis([1 kmax])
caxis([-1 1])
cb = colorbar;
cb.Label.String = 'k';
axis([-1 1 -1 1]*15)
xlabel('\xi (R_\oplus)');
ylabel('\zeta (R_\oplus)');

% plot( xi/sc, zeta/sc, '+k' )



