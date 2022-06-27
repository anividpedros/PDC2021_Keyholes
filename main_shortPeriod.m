% SOLUTION OF A PLANETARY ENCOUNTER: PDC 2021
% 7th IAA Planetary Defense Conference 2021
clc
clear
close all
format shortG

filename = 'main_shortPeriod.m';
filepath = matlab.desktop.editor.getActiveFilename;
currPath = filepath(1:(end-length(filename)));
cd(currPath)
addpath(genpath(pwd))

%% Constants
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

% % Bennu
% ast_id= '-101955';
% epoch = '2060-Sep-22 00:36';
% 
% % 2017 PDC
% cspice_furnsh( 'SPICEfiles/2017_PDC-merged-DE431.bsp' )
% ast_id= '-937001';
% epoch = '2027 July 20 15:00 TDB';
% et = cspice_str2et( epoch ) + 24*.715*3600;
% 
% % Apophis
% epoch = '2029-Apr-12 21:46:00.0000 TDB';
% et = cspice_str2et( epoch ) ;


cons.AU  = cspice_convrt(1, 'AU', 'KM');
cons.GMs = cspice_bodvrd( 'SUN', 'GM', 10 );
% cons.GMe = cspice_bodvrd( '399', 'GM', 1 );
cons.GMe = 398600.43543609593;
cons.Re  = 6378.140;

cons.yr  = 365.25 * 24 * 3600 ;
cons.Day = 3600*24; 

state_eat = cspice_spkezr( '399', et, 'ECLIPJ2000', 'NONE', '10' );
kep_eat = cspice_oscelt( state_eat, et, cons.GMs );
sma_eat = kep_eat(1)/(1-kep_eat(2));

DU = norm( state_eat(1:3) );
t0 = et;

% Initial condition
kepE_sma = kep_eat';
kepE_sma(1) = kep_eat(1)/(1-kep_eat(2));


%% Choose asteroid elements
% ast = 2;
% load('MAT-MCres-25-Oct-2021-multi-1e3.mat');
% id = 1 + 24*(ast-1);
% 
% OE0 = [a_t(id,1) e_t(id,1) i_t(id,1) O_t(id,1) o_t(id,1) 0]';
% 
% kep_opik_post = OE0;
% kep_opik_post(1) = kep_opik_post(1)*(1-kep_opik_post(2))*DU ; 
% kep_opik_post(7:8) = [t0; cons.GMs];

load('Horizons-2006MB14.mat')

% Exact encounter date
epoch  = '1985-Jun-28 00:40 TDB';
et_enc = cspice_str2et( epoch );

% Defining a "nominal" time
id_nom = 74;
date_cal = char(CalendarDateTDB(id_nom));
epoch = [date_cal(7:17) ' TDB'];
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




%% Numerical integration w constant planet elements
% Numerical Integration of point
eti = t0 ; % Initial ephemeris time for integration
X0  = cspice_conics(kep_opik_post, eti );
tv  = ( 0:0.01:100 )*cons.yr;
kep_planet = NaN(8,8);
GMvec = cons.GMs;
for i=2:9
    GMvec(i)          = cspice_bodvrd( [num2str(i-1) ], 'GM', 1);
    state_planet      = cspice_spkezr( [num2str(i-1) ],  eti, 'ECLIPJ2000', 'NONE', '10' );
    kep_planet(:,i-1) = cspice_oscelt( state_planet, eti, cons.GMs );
end
kep_planet(:,3) = kep_eat ;
GMvec(3) = cons.GMe;

cons_ode.t0        = eti ;
cons_ode.GM_vec    = GMvec ; % 1st element is Sun
cons_ode.IC_planet = kep_planet ; % 1 column per planet

tol = 1e-13;
options=odeset('RelTol',tol,'AbsTol',ones(1,6)*tol);
[t,X]=ode113(@(t,X) NBP_propagator(t,X,cons_ode),tv,X0,options);

% Postprocessing
% d_nbp   = zeros(length(tv),1);
kep_nbp = zeros(length(tv),8);
MOIDnbp = zeros(length(tv),1);
for i = 1:length(tv)
    kep0_nbp     = cspice_oscelt( X(i,:)', eti+tv(i), cons.GMs );
    kep_nbp(i,:) = kep0_nbp;
    kep_nbp(i,1) = kep0_nbp(1)/(1-kep0_nbp(2));
   
    MOIDnbp(i) = MOID_ORCCA_win( K2S(kepE_sma,cons.AU), K2S(kep_nbp(i,:),cons.AU) ) *cons.AU;
    
%     % Compute distance to Earth
%     xe   = cspice_conics(kep_eat, eti+tv(i) );
%     d_nbp(i) = norm(xe(1:3) - X(i,1:3)'); % From 4bp integration
end

% Plot orbit elements
plot_oe( kep_nbp, tv/cons.yr, 6, [0 100] )


%% Secular propagation
% Secular Model: Lagrange-Laplace
kepJ_sma = kep_planet(:,5);
kepJ_sma(1) = kepJ_sma(1)/(1-kepJ_sma(2));

cons_sec.OEp = kepJ_sma';
cons_sec.GMp = GMvec(6);
cons_sec.GMs = GMvec(1);

kep0 = kep_opik_post;
kep0_sma = kep_opik_post';
kep0_sma(1) = kep0(1)/(1-kep0(2));

OE0    = kep0_sma ;
% OE0(6) = OE0(6) - sqrt(cons.GMs/OE0(1)^3)*kep0_sma(7) ;
% cons_sec.OEp(7:8) = [];
secular_model_LL = secular_model_10BP_s2(OE0, cons_sec, 1);

kep_LL_t = zeros(length(tv),6);
% kept = kep0_sma;
% d_ll = zeros(length(tv),1);
% MOIDsec = d_ll;
for i = 1:length(tv)
    [~, kep_LL_t(i,:)] = drifted_oe_s2( secular_model_LL, tv(i), OE0, kepJ_sma' );
%     kep_LL_t(1)   = kep0_sma(1)*(1-kep_LL_t(i,2));

    MOIDsec(i) = MOID_ORCCA_win( K2S(kepE_sma,cons.AU), K2S(kep_LL_t(i,:),cons.AU) ) *cons.AU;
    
end

% Plot orbit elements
plot_oe( kep_LL_t, tv/cons.yr, 6, [0 30] )


%% Integrate short-period dynamics
% Hopes that it will be faster!!
tvs  = ( 0:0.01:1 )*cons.yr;

OEJ    = kepJ_sma ;
% OEJ(6) = OEJ(7) - OEJ(6)/sqrt(OEJ(8)/OEJ(1)^3) ;
% This is tau! We need sigma
OEJ(6) = OEJ(6) - sqrt(OEJ(8)/OEJ(1)^3)*OEJ(7) ;


EQ0 = Kep_2_Equinoctial( OE0 );

cons_aux.GMs = cons.GMs;
cons_aux.OEp = OEJ ;

cons_aux.GMcentral = cons.GMs;
cons_aux.GMthird   = GMvec(6);

secular_model_LL = secular_model_10BP_extra(kep0_sma, cons_sec, 1);
Xdot = ShortPeriod_prop(0, EQ0, cons_aux, secular_model_LL);

tol = 1e-13;
options=odeset('RelTol',tol,'AbsTol',ones(1,6)*tol);
[t,EQt]=ode113(@(t,X) ShortPeriod_prop_dev(t,X,cons_aux,secular_model_LL), tvs,EQ0,options);

% To-Dos:
% - Validate all variables defined here
%    Check that we are giving the right inputs
% - Validate new partials with finite differences


for i = 1:length(t)
    kep_sp(i,:) = Equi_2_Keplerian(EQt(i,:));
end

%%
% Plot orbit elements
plot_oe( kep_sp, tvs/cons.yr, 5, [0 30] )


%% Approximating the amplitude of the MOID evolution
% Amplitude in a-e-i

% Variables: kep_nbp, tv/cons.yr
% kep_LL_t, tv/cons.yr, 5, [0 30] )

per = 2*pi/sqrt(cons.GMs./kep_nbp(1).^3)/cons.yr ;
id_end = find( tv/cons.yr >= per*20, 1 );

F = figure(10); clf;

plot( tv/cons.yr, MOIDnbp, tv/cons.yr, MOIDsec )
grid on
xlabel('t (yr)')
ylabel('MOID')
hold on;

% Compute moving average
per = 2*pi/sqrt(cons.GMs./kep_nbp(1).^3) ;

for i=1:length(tv)
    i0 = i;
    i1 = find( tv >= (tv(i) + per*10), 1 );
    if isempty(i1); i1=length(tv); end
    
    mmean(i)= mean( MOIDnbp(i0:i1) );
    mstd(i) = std( MOIDnbp(i0:i1) );
    mstd(i) = max( MOIDnbp(i0:i1) ) - min( MOIDnbp(i0:i1) );
end
xp = tv/cons.yr;
% plot( xp,mmean,'r--', xp,mmean+2*mstd,'r-', xp,mmean-2*mstd,'r-' )
plot( xp,mmean,'r--', xp,mmean+mstd,'r-', xp,mmean-mstd,'r-' )

id_mstd = find( tv/cons.yr > 20,1 );
mean_std = mean( mstd(1:id_mstd) ) ;

plot( tv/cons.yr, MOIDsec-mean_std, 'g',...
    tv/cons.yr, MOIDsec+mean_std, 'g')


%% Repeat plot, but in Earth Radii
F = figure(11); clf;
sc = cons.Re ; 

plot( tv/cons.yr, MOIDnbp/sc, tv/cons.yr, MOIDsec/sc )
grid on
xlabel('t (yr)')
ylabel('MOID (Earth Radii)')
hold on;
% plot( xp,mmean/sc,'r--', xp,mmean/sc+mstd/sc,'r-', xp,mmean/sc-mstd/sc,'r-' )
% 
plot( tv/cons.yr, MOIDsec/sc-mean_std/sc, 'g',...
    tv/cons.yr, MOIDsec/sc+mean_std/sc, 'g')



%% SYSTEMATIC COMPARISON OF DIFFERENT POINTS

xiv = (-25:2:25)*cons.Re ;
zev = (-25:2:25)*cons.Re ;

perc_inside = nan(length(xiv),length(zev));
max_error   = nan(length(xiv),length(zev));
dMOID_sec   = nan(length(xiv),length(zev));
min_sec     = nan(length(xiv),length(zev));
min_nbp     = nan(length(xiv),length(zev));
amplitude   = nan(length(xiv),length(zev));


t0 = et;
plotting = 1;

% profile on

for i=1%:length(xiv)
    for j=1%:length(zev)
        
        xi0   =  -20*cons.Re; %kh_down_xi(ik(1));
        zeta0 =  20*cons.Re; %kh_down_zeta(ik(1));
        
        % Include if arrays are centered at 0
        % if i == 11 && j == 11
        %     continue
        % end
        
        %xi0 = xiv(i);
        %zeta0 = zev(j);
        
        fprintf('[i=%g/%g,j=%g/%g] - b={%g,%g}Re ',i,length(xiv),j,length(zev),xi0/cons.Re,zeta0/cons.Re)
        
        t0 = et;
        
        kep_opik_post = opik_var_post( state_eat, kep_eat, kep_ast, t0, xi0, zeta0, cons );
        if kep_opik_post(1) == 0 || sum(isnan(kep_opik_post))
            continue
        end
        
        t0 = et + 40*3600*24; % New definition to propagate numerically
        
        [kep_nbp, MOIDnbp, tv] = NBP_moid( kep_opik_post, kep_eat, t0, cons );
        id_mstd  = find( tv/cons.yr > 10,1 );
        mean_kep_nbp = mean( kep_nbp(1:id_mstd,:),1 )';
        mean_kep_nbp(1) = mean_kep_nbp(1)*(1-mean_kep_nbp(2));
        
        % [kep_sec, MOIDsec] = LL_moid( kep_opik_post, kep_eat, t0, tv, cons );
        [kep_sec, MOIDsec] = LL_moid( mean_kep_nbp, kep_eat, t0, tv, cons );
        
        % Individual generic postprocessing =========
        
        % Compute moving average---------------------
        % per = 2*pi/sqrt(cons.GMs./kep_nbp(1).^3) ;
        % for i=1:length(tv)
        %     i0 = i;
        %     i1 = find( tv >= (tv(i) + per*10), 1 );
        %     if isempty(i1); i1=length(tv); end
        %
        %     mmean(i)= mean( MOIDnbp(i0:i1) );
        %     mstd(i) = std( MOIDnbp(i0:i1) );
        %     mstd(i) = max( MOIDnbp(i0:i1) ) - min( MOIDnbp(i0:i1) );
        % end
        % xp = tv/cons.yr;
        % % plot( xp,mmean,'r--', xp,mmean+2*mstd,'r-', xp,mmean-2*mstd,'r-' )
        % plot( xp,mmean,'r--', xp,mmean+mstd,'r-', xp,mmean-mstd,'r-' )
        % id_mstd = find( tv/cons.yr > 20,1 );
        % mean_std = mean( mstd(1:id_mstd) ) ;
        % Max variation wrt initial ----------------
        id_mstd = find( tv/cons.yr > 10,1 );
        mean_std = 1*max( abs( MOIDnbp(1:id_mstd)-MOIDnbp(1) ) );
        %-------------------------------------------
        
        if plotting
            F = figure(1); clf;
            F.Position = [540 520 580 253];
            sc = cons.Re ;
            
            plot( tv/cons.yr, MOIDnbp/sc ); hold on;
            plot( tv/cons.yr, MOIDsec/sc, 'c' )
            grid on
            xlabel('t (yr)')
            ylabel('MOID (Earth Radii)')
            hold on;
            % plot( xp,mmean/sc,'r--', xp,mmean/sc+mstd/sc,'r-', xp,mmean/sc-mstd/sc,'r-' )
            
%             plot( tv/cons.yr, MOIDsec/sc-mean_std/sc, 'm--',...
%                 tv/cons.yr, MOIDsec/sc+mean_std/sc, 'm--')
            plot( tv([1 end])/cons.yr, MOIDnbp(1)/sc*[1 1], 'm--')
            axis([0 100 0 50])
            
        end
        
        % Accumulate metrics
        perc_inside(i,j) = 100*sum( abs(MOIDnbp - MOIDsec)< mean_std )/length(tv);
        max_error(i,j)   = max( abs(MOIDnbp - MOIDsec) );
        dMOID_sec(i,j)   = (MOIDsec(2)-MOIDsec(1))/tv(2);
        min_sec(i,j)     = min(MOIDsec);
        min_nbp(i,j)     = min(MOIDnbp);
        amplitude(i,j)   = 1.5*max( abs( MOIDnbp(1:id_mstd)-MOIDnbp(1) ) );
        
        % fprintf('[i=%g/%g,j=%g/%g] - b={%g,%g}Re -  perc=%.1f\n',i,length(xiv),j,length(zev),xi0/cons.Re,zeta0/cons.Re,perc_inside(i,j))
        fprintf(' -  perc=%.1f\n',perc_inside(i,j))
        
    end
end

% profile viewer
% save('ApophisBplanemats-25Re-post.mat', 'xiv','zev','amplitude','perc_inside','max_error','dMOID_sec','min_sec','min_nbp')

%%
sc = cons.Re;
xtck = -200:20:200 ;

%
F = figure(44); clf;

contourf( xiv/sc, zev/sc, dMOID_sec'/sc*cons.yr, 20, 'LineStyle', 'none' );
hold on;
% contour( xiv/sc, zev/sc, dMOID_sec'/sc*cons.yr, [-0.10 0 0.1],'k', 'ShowText', 'on' );
grid on;
hold on;
xlabel('\xi (Re)')
ylabel('\zeta (Re)')
% caxis([-0.15 0.10])
colorbar
% xticks(-200:50:200)
% yticks(-200:50:200)
title('Secular MOID rate (Re/yr)')
axis equal
F.Position([3 4]) = [367 274];
% axis([-180 180 -180 180])

% -------------------------------------------
F = figure(45); clf; 

contourf( xiv/sc, zev/sc, perc_inside', 'LineStyle', 'none' );
grid on;
hold on;
xlabel('\xi (Re)')
ylabel('\zeta (Re)')
colorbar
% xticks(-200:40:200)
% yticks(-200:40:200)
title('Percentage Inside Bounds')
axis equal
% axis([-180 180 -180 180])
F.Position([3 4]) = [367 274];
caxis([70 90])

%-------------------------------------------
F = figure(46); clf; 

contourf( xiv/sc, zev/sc, max_error'/sc, 'LineStyle', 'none' );
grid on;
hold on;
xlabel('\xi (Re)')
ylabel('\zeta (Re)')
colorbar
% xticks(-200:40:200)
% yticks(-200:40:200)
title('Maximum Error (Re)')
axis equal
% axis([-180 180 -180 180])
F.Position([3 4]) = [367 274];

%-------------------------------------------

F = figure(49); clf;
% subplot(1,2,1)
contourf( xiv/sc, zev/sc, min_sec'/sc );
grid on;
hold on;
xlabel('\xi (Re)')
ylabel('\zeta (Re)')
colorbar
% xticks(-200:40:200)
% yticks(-200:40:200)
title('Secular minimum MOID (Re)')


F = figure(47); clf; 

% subplot(1,2,2)
contourf( xiv/sc, zev/sc, (min_nbp')/sc, 20, 'LineStyle', 'none' );
hold on
contour( xiv/sc, zev/sc, (min_nbp')/sc, [0:50:200],'k', 'ShowText', 'on' );
grid on;
hold on;
xlabel('\xi (Re)')
ylabel('\zeta (Re)')
colorbar
% xticks(-200:50:200)
% yticks(-200:50:200)
title('Numerical minimum MOID (Re)')
axis equal
F.Position([3 4]) = [367 274];

%-------------------------------------------
F = figure(48); clf;

contourf( xiv/sc, zev/sc, amplitude'/sc, 20, 'LineStyle', 'none' );
hold on
% contour( xiv/sc, zev/sc, amplitude'/sc, [0:5:50],'k', 'ShowText', 'on' );

grid on;
hold on;
xlabel('\xi (Re)')
ylabel('\zeta (Re)')
c = colorbar;
% caxis([0 8])
% xticks(-200:50:200)
% yticks(-200:50:200)
title('Amplitude (Re)')
axis equal
% axis([-180 180 -180 180])
F.Position([3 4]) = [367 274];


%%

F = figure(50); clf; 

% subplot(1,2,2)
contourf( xiv/sc, zev/sc, (max_nbp')/sc, 20, 'LineStyle', 'none' );
hold on
contour( xiv/sc, zev/sc, (max_nbp')/sc, [0:50:200],'k', 'ShowText', 'on' );
grid on;
hold on;
xlabel('\xi (Re)')
ylabel('\zeta (Re)')
colorbar
% xticks(-200:50:200)
% yticks(-200:50:200)
title('Numerical maximum MOID (Re)')
axis equal
F.Position([3 4]) = [367 274];



%% Variation in semi-major axis
for i=15%:length(xiv)
    for j=1%:length(zev)
%         xi0   =  10*cons.Re; %kh_down_xi(ik(1));
%         zeta0 =  50*cons.Re; %kh_down_zeta(ik(1));
        
        % Include if arrays are centered at 0
        % if i == 11 && j == 11
        %     continue
        % end
        
        xi0 = xiv(i);
        zeta0 = zev(j);
         
        fprintf('[i=%g/%g,j=%g/%g] - b={%g,%g}Re \n',i,length(xiv),j,length(zev),xi0/cons.Re,zeta0/cons.Re)
        
        t0 = et;
        
        kep_opik_post = opik_var_post( state_eat, kep_eat, kep_ast, t0, xi0, zeta0, cons )
        
        sma0 = kep_ast(1)/(1-kep_ast(2));
        sma1 = kep_opik_post(1)/(1-kep_opik_post(2));
        dA(i,j) = sma1 - sma0;
    end
end

%%
F = figure(60);
contourf(xiv/sc,zev/sc,dA')

F = figure(61);
contourf(xiv/sc,zev/sc,100*dA'./sma0)



%% Functions definition
%------------------------------------
% Generate Opik variables and post-encounter keplerian
function kep_opik_post = opik_var_post( state_eat, kep_eat, kep_ast, t0, xi0, zeta0, cons )

state_ast = cspice_conics(kep_ast, t0);

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
dt_per  = (MA_per - kep_ast_O(6))/n_ast_O; dt_per_day = dt_per/86400;

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
% r0   = auxR'*[xi0; 0; zeta0];

% Modify initial B-plane conditions
% To-Do: 
% - Generate circles
% - Propagate different values
% 
% xi0   = 8e5; %kh_down_xi(ik(1));
% zeta0 = 8e5; %kh_down_zeta(ik(1));

kep_opik_post = opik_bplane_2_oe( theta,phi,zeta0/DU,xi0/DU,U_nd,phi,longp,ap )';
% t0 = t0 + 40*3600*24;

if isreal(kep_opik_post)
    % % 7. Heliocentric orbit elements from post-encounter coordinates (Opik formulae)
    % longp = mod( kep_eat(4)+kep_eat(5)+kep_eat(6), 2*pi ) ; % In general sense should be longitude
    % longp = atan2( state_eat(2),state_eat(1) );
    % kep_opik_post = opik_bplane_2_oe( theta1,phi1,zeta1,xi1,U_nd,phi,longp,ap )';
    
    kep_opik_post(1) = kep_opik_post(1)*(1-kep_opik_post(2))*DU ;
    % 'opik_bplane_2_oe.m' function returns sma in first element
    
    kep_opik_post(6) = TA_2_MA(kep_opik_post(6),kep_opik_post(2));
    % 'opik_bplane_2_oe.m' function returns true anomaly 6th element
    
    kep_opik_post(7:8) = [t0 + dt_per;
        cons.GMs];

else
    fprintf('Opik Theory Failed... moving on\n')
    kep_opik_post = zeros(8,1);
end

end

%------------------------------------
% Secular propagation
function [kep_LL_t, MOIDsec] = LL_moid( kep_opik_post, kep_eat, eti, tv, cons )

kep_planet = NaN(8,8);
GMvec = cons.GMs;
for i=2:9
    GMvec(i)          = cspice_bodvrd( [num2str(i-1) ], 'GM', 1);
    state_planet      = cspice_spkezr( [num2str(i-1) ],  eti, 'ECLIPJ2000', 'NONE', '10' );
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

kep0 = kep_opik_post;
kep0_sma = kep_opik_post';
kep0_sma(1) = kep0(1)/(1-kep0(2));

OE0    = kep0_sma ;
% OE0(6) = OE0(6) - sqrt(cons.GMs/OE0(1)^3)*kep0_sma(7) ;
% cons_sec.OEp(7:8) = [];
secular_model_LL = secular_model_10BP_s2(OE0, cons_sec, 1);

kep_LL_t = zeros(length(tv),6);
% kept = kep0_sma;
% d_ll = zeros(length(tv),1);
% MOIDsec = d_ll;

% Initial condition
kepE_sma = kep_eat';
kepE_sma(1) = kep_eat(1)/(1-kep_eat(2));

MOIDsec = zeros(length(tv),1);
for i = 1:length(tv)
    [~, kep_LL_t(i,:)] = drifted_oe_s2( secular_model_LL, tv(i), OE0, kepJ_sma' );
%     kep_LL_t(1)   = kep0_sma(1)*(1-kep_LL_t(i,2));
    MOIDsec(i) = MOID_ORCCA_win( K2S(kepE_sma,cons.AU), K2S(kep_LL_t(i,:),cons.AU) ) *cons.AU;
    
end

end



%------------------------------------
% Numerical integration 
function [kep_nbp, MOIDnbp, tv] = NBP_moid( kep_opik_post, kep_eat, t0, cons )

% Numerical Integration of point
eti = t0 ; % Initial ephemeris time for integration
X0  = cspice_conics(kep_opik_post, eti );
tv  = ( 0:0.01:100 )*cons.yr;
kep_planet = NaN(8,8);
GMvec = cons.GMs;
for i=2:9
    GMvec(i)          = cspice_bodvrd( [num2str(i-1) ], 'GM', 1);
    state_planet      = cspice_spkezr( [num2str(i-1) ],  eti, 'ECLIPJ2000', 'NONE', '10' );
    kep_planet(:,i-1) = cspice_oscelt( state_planet, eti, cons.GMs );
end
kep_planet(:,3) = kep_eat ;
GMvec(3) = cons.GMe;

cons_ode.t0        = eti ;
cons_ode.GM_vec    = GMvec ; % 1st element is Sun
cons_ode.IC_planet = kep_planet ; % 1 column per planet

tol = 1e-13;
options=odeset('RelTol',tol,'AbsTol',ones(1,6)*tol);
[t,X]=ode113(@(t,X) NBP_propagator(t,X,cons_ode),tv,X0,options);

% Postprocessing
% d_nbp   = zeros(length(tv),1);
kep_nbp = zeros(length(tv),8);
MOIDnbp = zeros(length(tv),1);

% Initial condition
kepE_sma = kep_eat';
kepE_sma(1) = kep_eat(1)/(1-kep_eat(2));

for i = 1:length(tv)
    kep0_nbp     = cspice_oscelt( X(i,:)', eti+tv(i), cons.GMs );
    kep_nbp(i,:) = kep0_nbp;
    kep_nbp(i,1) = kep0_nbp(1)/(1-kep0_nbp(2));
   
    MOIDnbp(i) = MOID_ORCCA_win( K2S(kepE_sma,cons.AU), K2S(kep_nbp(i,:),cons.AU) ) *cons.AU;
    
%     % Compute distance to Earth
%     xe   = cspice_conics(kep_eat, eti+tv(i) );
%     d_nbp(i) = norm(xe(1:3) - X(i,1:3)'); % From 4bp integration
end


end

%------------------------------------
% Plotting orbit elements 
function plot_oe( kep_nbp, tv, idfig, xlims )

labels = {'a','e','i','\Omega','\omega','MA'};
F = figure(idfig);

kep_nbp(:,3:6) = wrapTo2Pi(kep_nbp(:,3:6));

for i=1:6

    subplot(2,3,i)
    plot(tv, kep_nbp(:,i))
    hold on; grid on
    xlabel('Time')
    ylabel(labels{i})
    xlim(xlims)
    
end

end











