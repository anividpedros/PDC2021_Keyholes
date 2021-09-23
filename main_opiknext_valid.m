% SOLUTION OF A PLANETARY ENCOUNTER: PDC 2021
% 7th IAA Planetary Defense Conference 2021
clc
clear
close all
format shortG

filename = 'main_opiknext_valid.m';
filepath = matlab.desktop.editor.getActiveFilename;
currPath = filepath(1:(end-length(filename)));
cd(currPath)
addpath(genpath(pwd))

%% ===== Description ===== 
% Original script: main_sect3_multimethod_dist.m
% Extension: Validate post-encounter values
%  =======================


%% Constants declaration
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


%% CASE 1999 AN10
U  = 0.884 ;
xi = 0.000246 ;
theta = 105.3 *pi/180;
phi   = 41.3  *pi/180;
zetav = -.21:0.00005:.21 ;

t0 = 0;
h = 0;
m  = cons.GMe/cons.GMs ;
longp = 0;
ap = 1;

for i=1:length(zetav)
    [~,theta1,phi1,xi1(i),zeta1(i)] = opik_next(U,theta,phi,xi,zetav(i),t0,h,m);
    kep_opik_post(i,:) = opik_bplane_2_oe( theta1,phi1,zeta1(i),xi1(i),U,phi,longp,ap )';
end

F = figure(1);
plot( kep_opik_post(:,1), kep_opik_post(:,2), '.' )
axis([1.35 1.55 0.5 0.6])
grid on
xlabel('a (au)')
ylabel('e')

figure(2)
Per = 2*pi ./sqrt( cons.GMs./(kep_opik_post(:,1)*cons.AU).^3 ) ./cons.yr;
plot( Per,zetav, '.' )
grid on
xlabel('\Per (yr)')
ylabel('\zeta pre')

figure(3)
plot( zetav,xi1 )
grid on
xlabel('\zeta pre')
ylabel('\xi post')


%% Validate keyholes function
c_nd = (cons.GMe/cons.GMs)/U^2;

c2 = c_nd*c_nd;
U2 = U*U;

st = sin(theta);
ct = cos(theta);

h = [7 10 11];
k = [13 17 19];

circ = [];
for i=1:3
        
    a0p(i) = (k(i)/h(i))^(2/3);
    if sum( find(a0p(i) == a0p(1:i-1)) ); continue; end
    
    ct0p = ( 1-U2-ap/a0p(i) )/2/U ;
    if abs(ct0p) > 1; continue; end
    st0p = sqrt(1 - ct0p^2);
    
    % Cirle center location and radius
    D = (c_nd*st)/(ct0p - ct);
    R = abs( c_nd*st0p/(ct0p - ct) );
    
    circ = [circ;
        k(i) h(i) D R];
    
end

% Plot the circles
nr = size(circ,1);
co = winter(22);
thv= 0:(2*pi/99):2*pi ;
sc = cons.Re/cons.AU;

F = figure(3); clf;

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

grid on
axis([-40 40 -40 40])
axis square

plot(xi*[1 1]/sc, zetav([1 end])/sc, 'k')

xlabel('\xi (R_\oplus)');
ylabel('\zeta (R_\oplus)');


% Compute keyholes
DU = cons.AU;
circles = circ;

for ic=1:nr
    
    k = circles(ic,1);
    h = circles(ic,2);
    D = circles(ic,3);
    R = circles(ic,4);
    
    [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
        two_keyholes_Val03(k, h, D/sc, R/sc, U, theta, phi, m,0,DU);

    % cc = co(k,:);
    cc = [0 0 0];
    plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
        'Color',cc);
    plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
        'Color',cc);
end


%% CASE 1997 XF11
U     = 0.459 ;
theta = 84.0 *pi/180 ;
phi   = 99.5 *pi/180 ;
zetav = -.21:0.00005:.21 ;

t0 = 0;
h = 0;
m  = cons.GMe/cons.GMs ;
longp = 0;
ap = 1;


%% Potential Additional validation
% - From heliocentric to Opik elements
% - Check outgoing heliocentric if resonant encounter occurs
%  ~ Final sma, ecc checked in the figure
% - 


%% From heliocentric to Opik elements
TU = sqrt(DU^3/cons.GMs);
% ast_id= '137108';
% epoch = '2021 October 20 TDB';
% epoch = '2027 August 7 TDB';
% tv = (-24*5:.5:0.5) * 3600 ;

load('Horizons-1999AN10.mat')
thv = zeros(length(JDTDB),1);
phv = thv;

for i=1:length(JDTDB)
    
    % et = cspice_str2et( epoch ) + 24*.715*3600;
    % et = cspice_str2et( epoch ) + tv(i);
    et = cspice_str2et( ['JD' num2str(JDTDB(1))] );
    
    state_eat = cspice_spkezr( '399', et, 'ECLIPJ2000', 'NONE', '10' );
    kep_eat = cspice_oscelt( state_eat, et, cons.GMs );
    
    state_ast = [X(i); Y(i); Z(i); VX(i); VY(i); VZ(i)];
%     state_ast(1:3) = state_ast(1:3)*DU;
%     state_ast(4:6) = state_ast(4:6)*(DU/(24*3600));
%     state_ast = cspice_spkezr( ast_id, et, 'ECLIPJ2000', 'NONE', '10' );
%     kep_ast = cspice_oscelt( state_ast, et, cons.GMs );
    state_ast_t = HillRot(state_eat, state_ast);

    rt(i)= norm( state_ast_t(1:3) );
    V(i) = norm( state_ast_t(4:6) )/(DU/TU);
    thv(i) = acos( state_ast_t(5)/V(i) );
    phv(i) = atan2( state_ast_t(4), state_ast_t(6) );
    
end


figure(12);
subplot(2,1,1)
plot( JDTDB/86400, thv*180/pi ); ylabel('\theta')
grid on; hold on
subplot(2,1,2)
plot( JDTDB/86400, phv*180/pi ); ylabel('\phi')
grid on; hold on

figure(13)
subplot(2,1,1)
plot( JDTDB/86400, V )
subplot(2,1,2)
plot( JDTDB/86400, rt/cons.AU )