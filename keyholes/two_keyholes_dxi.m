function [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta,dx_down,dx_up] = two_keyholes(k,h,D,R,U,theta,phi,m,t0,DU, type,radius,varargin)
% COMPUTE KEYHOLES FOR A GIVEN RESONANCE
%
% For 3 options of time-varying MOID:
% Defined by input variable "type"
% [1] Constant MOID
% [2] Secular evolution
% [3] Numerical integration
%

if type == 2
    longp = varargin{1};
    ap    = varargin{2};
    cons  = varargin{3};
    kepE_sma = varargin{4};
    cons_sec = varargin{5};
    AU = cons.AU;
    
elseif type == 3
    longp = varargin{1};
    ap    = varargin{2};
    cons  = varargin{3};
    kepE_sma = varargin{4};
    cons_ode = varargin{5};
    AU = cons.AU;
    tol = 1e-13;
    options=odeset('RelTol',tol,'AbsTol',ones(1,6)*tol);
    
end


%%

% Convert to dimensionless units (au, but here we use the real distance of
% the Earth at the instant of close encounter)
RE_km = 6378.140;
% DU = 152003856.97586098;
RE_au = RE_km/DU;
D_au = D*RE_au;
R_au = R*RE_au;

%% The following part should be included in a loop for all -b_Earth < xi < b_Earth

c = m/U^2;
% bEarth_au = RE_au*sqrt(1 + 2*c/RE_au);
if radius == -1
    bEarth_au = RE_au*sqrt(1 + 2*c/RE_au);
else
    bEarth_au = RE_au*radius;
end

% Get xi, zeta for the previous ranges of alpha
nkh = 50;
[xi_up, zeta_up]     = res_circle(linspace(0,pi,nkh), D_au,R_au);
[xi_down, zeta_down] = res_circle(linspace(-pi,0,nkh),D_au,R_au);

% For each of these values, compute the keyhole width.

% Keyhole - bottom
zeta_edges = zeros(nkh*2,2);
xi_edges   = zeros(nkh*2,2);
dx_edges   = zeros(nkh*2);

for i = 1:nkh*2
    
    % Initial conditions of the flyby
    if i <= nkh
        xi   = xi_down(i);
        zeta = zeta_down(i);
    else
        xi   = xi_up(i-nkh);
        zeta = zeta_up(i-nkh);
    end
    
    if sqrt(xi^2+zeta^2) < bEarth_au
        zeta_edges(i,:) = NaN;
        xi_edges(i,:)   = NaN;
        continue
    end
    
    % COMPUTING THE MOID AT SUBSEQUENT ENCOUNTER: 3 METHODS
    [~,theta1,phi1,xi1,zeta1,~] = opik_next(U,theta,phi,xi,zeta,t0,h,m);
    
    switch type
        
        case 1  %------- Case constant MOID ---------
            dx = 0;
            dx_edges(i) = dx;
            
        case 2  %------- Case secular model ---------
            
            % Post-Encounter Heliocentric Orbit Elements
            kep_opik_post = opik_bplane_2_oe( theta1,phi1,zeta1,xi1,U,phi,longp,ap)';
            if isnan(kep_opik_post) | ~isreal(kep_opik_post) 
                dx = nan ;
                dx_edges(i) = nan;
                
                zeta_edges(i,:) = NaN;
                xi_edges(i,:)   = NaN;
                continue
            end
            
            kep_opik_post(1) = kep_opik_post(1)*(1-kep_opik_post(2))*DU ; % 'opik_bplane_2_oe.m' function returns sma in first element
            kep_opik_post(6) = TA_2_MA(kep_opik_post(6),kep_opik_post(2));% 'opik_bplane_2_oe.m' function returns true anomaly 6th element
            kep_opik_post(7:8) = [t0;
                cons.GMs];
            
            kep0_sma = kep_opik_post';
            kep0_sma(1) = kep0_sma(1)/(1-kep0_sma(2));
            
            moid0 = MOID_ORCCA_win(K2S(kepE_sma,AU), K2S(kep0_sma,AU));
%             moid0 = ComputeMOID_mex_MAC(K2S(kepE_sma,AU), K2S(kep0_sma,AU));
            
            At = k * cons.yr;
            % SECULAR PROPAGATION
            secular_model_LL = secular_model_10BP_s2(kep0_sma, cons_sec, 1);
            [~, kep0_LL_t] = drifted_oe_s2( secular_model_LL, At, kep0_sma, cons_sec.OEp);
            
            moid1 = MOID_ORCCA_win( K2S(kepE_sma,AU), K2S(kep0_LL_t,AU));
%             moid1 = ComputeMOID_mex_MAC(K2S(kepE_sma,AU), K2S(kep0_LL_t,AU));
            dx = moid1 - moid0;
            dx_edges(i) = dx;
            
        case 3  %------- Case NBP integration model ---------
            
            % Post-Encounter Heliocentric Orbit Elements
            kep_opik_post = opik_bplane_2_oe( theta1,phi1,zeta1,xi1,U,phi,longp,ap)';            
            if isnan(kep_opik_post) | ~isreal(kep_opik_post) 
                dx = nan ;
                dx_edges(i) = nan;
                
                zeta_edges(i,:) = NaN;
                xi_edges(i,:)   = NaN;
                continue
            end
            
            kep_opik_post(1) = kep_opik_post(1)*(1-kep_opik_post(2))*DU ; % 'opik_bplane_2_oe.m' function returns sma in first element
            kep_opik_post(6) = TA_2_MA(kep_opik_post(6),kep_opik_post(2));% 'opik_bplane_2_oe.m' function returns true anomaly 6th element
            kep_opik_post(7:8) = [t0;
                cons.GMs];
            
            kep0_sma = kep_opik_post';
            kep0_sma(1) = kep0_sma(1)/(1-kep0_sma(2));
            
            moid0 = MOID_ORCCA_win(K2S(kepE_sma,AU), K2S(kep0_sma,AU));
%             moid0 = ComputeMOID_mex_MAC(K2S(kepE_sma,AU), K2S(kep0_sma,AU));
            
            At = k * cons.yr;
            
            % NUMERICAL PROPAGATION
            X0  = cspice_conics(kep_opik_post, cons_ode.t0 );
            tv  = [0 At];
            [~,X]=ode113(@(t,X) NBP_propagator(t,X,cons_ode),tv,X0,options);
            
            kep0_nbp = cspice_oscelt( X(end,:)', cons_ode.t0+tv(end), cons.GMs );
            kep0_nbp(1) = kep0_nbp(1)/(1-kep0_nbp(2));
            
            moid1 = MOID_ORCCA_win( K2S(kepE_sma,AU), K2S(kep0_nbp,AU) );
%             moid1 = ComputeMOID_mex_MAC(K2S(kepE_sma,AU), K2S(kep0_nbp,AU));
            
            dx = moid1 - moid0;
            dx_edges(i) = dx;
            
    end
    

    % Check if xi1 <= b_Earth
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
    xi_edges(i,:)   = [xi,xi];

end

kh_down_xi   = xi_edges(1:nkh,:);
kh_down_zeta = zeta_edges(1:nkh,:);

kh_up_xi     = xi_edges(nkh+1:end,:);
kh_up_zeta   = zeta_edges(nkh+1:end,:);

dx_down      = dx_edges(1:nkh);
dx_up        = dx_edges(nkh+1:end);

end
