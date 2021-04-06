function [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = two_keyholes_dxi_sec(k,h,D,R,U,theta,phi,m,t0,DU,longp,ap,cons, kepE_sma,cons_ode)
%% COMPUTE KEYHOLES FOR A GIVEN RESONANCE

% Convert to dimensionless units (au, but here we use the real distance of
% the Earth at the instant of close encounter)
RE_km = 6378.140;
% DU = 152003856.97586098;
RE_au = RE_km/DU;
D_au = D*RE_au;
R_au = R*RE_au;

AU = cons.AU;

tol = 1e-13;
options=odeset('RelTol',tol,'AbsTol',ones(1,6)*tol);

%% The following part should be included in a loop for all -b_Earth < xi < b_Earth

c = m/U^2;
bEarth_au = RE_au*sqrt(1 + 2*c/RE_au);

% Get xi, zeta for the previous ranges of alpha
nkh = 20;
[xi_up,zeta_up] = res_circle(linspace(0,pi,nkh),D_au,R_au);
[xi_down,zeta_down] = res_circle(linspace(-pi,0,nkh),D_au,R_au);

% For each of these values, compute the keyhole width.

% Keyhole - bottom
zeta_edges = nan(nkh,2);
xi_edges = nan(nkh,2);
for i = 1:nkh
    
    xi = xi_down(i);
    zeta = zeta_down(i);
    
%------ New Code ------    
    if sqrt(xi^2+zeta^2) < bEarth_au
        continue
    end
    
    % Check if xi1 <= b_Eart

    [~,theta1,phi1,xi1,zeta1,~] = opik_next(U,theta,phi,xi,zeta,t0,h,m);
    % Post-Encounter Heliocentric Orbit Elements
    kep_opik_post = opik_bplane_2_oe( theta1,phi1,zeta1,xi1,U,phi,longp,ap)';
    
    kep_opik_post(1) = kep_opik_post(1)*(1-kep_opik_post(2))*DU ; % 'opik_bplane_2_oe.m' function returns sma in first element
    kep_opik_post(6) = TA_2_MA(kep_opik_post(6),kep_opik_post(2));% 'opik_bplane_2_oe.m' function returns true anomaly 6th element
    % !!!!!!! INPUT EPHEMERIS TIME IF NEEDED FOR NUM
    kep_opik_post(7:8) = [t0;
                        cons.GMs];

    kep0_sma = kep_opik_post';
    kep0_sma(1) = kep0_sma(1)/(1-kep0_sma(2));

    moid0 = MOID_ORCCA_win(K2S(kepE_sma,AU), K2S(kep0_sma,AU));
    
    At = k * cons.yr;
    
    % NUMERICAL PROPAGATION
    X0  = cspice_conics(kep_opik_post, cons_ode.t0 );
    tv  = [0 At];
    [~,X]=ode113(@(t,X) NBP_propagator(t,X,cons_ode),tv,X0,options);

    kep0_nbp = cspice_oscelt( X(end,:)', cons_ode.t0+tv(end), cons.GMs );
    kep0_nbp(1) = kep0_nbp(1)/(1-kep0_nbp(2));
   
    moid1 = MOID_ORCCA_win( K2S(kepE_sma,AU), K2S(kep0_nbp,AU) );

    dx = moid1 - moid0;
%------------------------
    
    if (abs(xi1 + dx) > bEarth_au)
        zeta_edges(i,:) = NaN;
        xi_edges(i,:)   = NaN;
        continue
    end
    xi2 = xi1 + dx;
    
    % Assign tentative dzeta1, dzeta2
    dz1 = 1E-6; dz2 = 1E-6;
    zetapp1 = zeta - dz1;
    zetapp2 = zeta + dz2;
    
    % Find the value of zeta leading to a direct impact
    try
        [zeta0,~] = fzero(@(zeta) opik_next(U,theta,phi,xi,zeta,t0,h,m),[zetapp1,zetapp2]);
    catch
        %disp(['fzero error at i =',i]);
        zeta_edges(i,:) = NaN;
        xi_edges(i,:)   = NaN;
        continue
    end
    
    % Compute keyhole edges
    zeta2_edges = [sqrt(bEarth_au^2 - xi2^2), -sqrt(bEarth_au^2 - xi2^2)];
    
    [~,theta1,phi1,xi1,zeta1] = opik_next(U,theta,phi,xi,zeta0,t0,h,m);
    dz2dz = dzeta2dzeta(U,theta,phi,xi,zeta0,m,h,theta1,phi1,xi1,zeta1);
    zeta_edges(i,:) = zeta0 + zeta2_edges/dz2dz;
    xi_edges(i,:) = [xi,xi];

end
kh_down_xi = xi_edges;
kh_down_zeta = zeta_edges;

% Keyhole - top
zeta_edges = nan(nkh,2);
xi_edges = nan(nkh,2);
for i = 1:nkh
    
    xi = xi_up(i);
    zeta = zeta_up(i);
    
    if sqrt(xi^2+zeta^2) < bEarth_au
        continue
    end
    
    % Check if xi1 <= b_Eart

    [~,theta1,phi1,xi1,zeta1,~] = opik_next(U,theta,phi,xi,zeta,t0,h,m);
    % Post-Encounter Heliocentric Orbit Elements
    kep_opik_post = opik_bplane_2_oe( theta1,phi1,zeta1,xi1,U,phi,longp,ap)';
    
    kep_opik_post(1) = kep_opik_post(1)*(1-kep_opik_post(2))*DU ; % 'opik_bplane_2_oe.m' function returns sma in first element
    kep_opik_post(6) = TA_2_MA(kep_opik_post(6),kep_opik_post(2));% 'opik_bplane_2_oe.m' function returns true anomaly 6th element
    % !!!!!!! INPUT EPHEMERIS TIME IF NEEDED FOR NUM
    kep_opik_post(7:8) = [t0;
                        cons.GMs];

    kep0_sma = kep_opik_post';
    kep0_sma(1) = kep0_sma(1)/(1-kep0_sma(2));

    moid0 = MOID_ORCCA_win(K2S(kepE_sma,AU), K2S(kep0_sma,AU));
    
    At = k * cons.yr;
    
    % NUMERICAL PROPAGATION
    X0  = cspice_conics(kep_opik_post, cons_ode.t0 );
    tv  = [0 At];
    [~,X]=ode113(@(t,X) NBP_propagator(t,X,cons_ode),tv,X0,options);

    kep0_nbp = cspice_oscelt( X(end,:)', cons_ode.t0+tv(end), cons.GMs );
    kep0_nbp(1) = kep0_nbp(1)/(1-kep0_nbp(2));
   
    moid1 = MOID_ORCCA_win( K2S(kepE_sma,AU), K2S(kep0_nbp,AU) );

    dx = moid1 - moid0;
    
    if (abs(xi1 + dx) > bEarth_au)
        zeta_edges(i,:) = NaN;
        xi_edges(i,:)   = NaN;
        continue
    end
    xi2 = xi1 + dx;
    
    % Assign tentative dzeta1, dzeta2
    dz1 = 1E-6; dz2 = 1E-6;
    zetapp1 = zeta - dz1;
    zetapp2 = zeta + dz2;
    
%     zeta2 = @(zeta) opik_next(U,theta,phi,xi,zeta,t0,h);
%     zv = zetapp1:(zetapp2-zetapp1)/999:zetapp2 ;
%     for iz = 1:length(zv)
%         zv2(iz) = zeta2(zv(iz));
%     end
%     figure;
%     plot(zv-zeta,zv2); grid on;
        
    % Find the value of zeta leading to a direct impact
    try
        [zeta0,~] = fzero(@(zeta) opik_next(U,theta,phi,xi,zeta,t0,h,m),[zetapp1,zetapp2]);
    catch
        %disp(['fzero error at i =',i]);
        zeta_edges(i,:) = NaN;
        xi_edges(i,:)   = NaN;
        continue
    end
    
    % Compute keyhole edges
    zeta2_edges = [sqrt(bEarth_au^2 - xi2^2), -sqrt(bEarth_au^2 - xi2^2)];
    
    [~,theta1,phi1,xi1,zeta1] = opik_next(U,theta,phi,xi,zeta0,t0,h,m);
    dz2dz = dzeta2dzeta(U,theta,phi,xi,zeta0,m,h,theta1,phi1,xi1,zeta1);
    zeta_edges(i,:) = zeta0 + zeta2_edges/dz2dz;
    xi_edges(i,:) = [xi,xi];

end
kh_up_xi = xi_edges;
kh_up_zeta = zeta_edges;

end
