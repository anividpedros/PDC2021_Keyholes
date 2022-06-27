function secular_model = secular_model_10BP(OE0, cons, pl_list)

% Adapted from main_oe_driftsPlanets_PlanetsModels
% List of Outputs:

% ev: Matrix of eigenvectors of H and K
% iv: Matrix of eigenvectors of P and Q
% g:  Frequencies of the H-K eigenvalue problem
% f:  Frequencies of the P-Q eigenvalue problem
% phaseE: Integration constants set to match ICs at t0
% phaseI: Integration constants set to match ICs at t0
% p: Position of the particle in the secular model 
%    (ordered by semi-major axis!)
% Ejk: Drifts in sigma
% s0:  Initial sigmas

% Evaluate secular solution
%     H = eji*sin( gi*t + Betai );
%     K = eji*cos( gi*t + Betai );
%     P = Iji*sin( fi*t + Gammai );
%     Q = Iji*cos( fi*t + Gammai );
%     S = Ejk*t + s0 ;

%% Reading inputs
OEp     = cons.OEp;
mu      = cons.GMs;
GMp_GMs = cons.GMp/mu ;

% Select Planets
% pl_list = 1:length(cons.GMp);
% pl_list = 5;

np  = length(pl_list);
OEp = OEp(pl_list,:);
GMp_GMs = GMp_GMs(pl_list);


%% Adding extra entry for the particle!
OEp = [OEp; OE0];
i0  = find( OEp(:,1) > OE0(1) , 1 );

OEp = sortrows( OEp, 1 );

GMp_GMs = GMp_GMs([1:i0 i0:np]);
GMp_GMs(i0) = 0;

%% N-body Lagrange-Laplace

dLaplace_coeff = @(p,s,j,al) cos(j*p)./(1-2*al*cos(p)+al^2).^(s) /pi ;
np = length( GMp_GMs );

b_1_32_jk = zeros(np);
b_2_32_jk = zeros(np);
D_b_0_12_jk = zeros(np);

al_jk     = zeros(np);
al_bar_jk = ones(np);
mk_mcj    = zeros(np);
nj        = zeros(np);

EQp = zeros(np,6);
for j = 1:np
    for k = 1:(j-1) % In all cases j > k
        
        al_jk(j,k) = OEp(k,1) / OEp(j,1);
        al = al_jk(j,k);
        
        b_1_32 = integral2( @(p) dLaplace_coeff(p, 1.5, 1, al), 0, 2*pi );
        b_2_32 = integral2( @(p) dLaplace_coeff(p, 1.5, 2, al), 0, 2*pi );
        
        b_1_32_jk(j,k) = b_1_32;
        b_2_32_jk(j,k) = b_2_32;
        
        b_0_32 = integral2( @(p) dLaplace_coeff(p,  1.5, 0, al), 0, 2*pi );
        D_b_0_12 = (b_1_32 - al*b_0_32);
        D_b_0_12_jk(j,k) = D_b_0_12 ;
        
        mk_mcj(j,k) = GMp_GMs(k)/(1 + GMp_GMs(j));
        mk_mcj(k,j) = GMp_GMs(j)/(1 + GMp_GMs(k));
        
    end
    
    al_jk(1:(j-1),j)     = al_jk(j,1:(j-1)); 
    al_bar_jk(1:(j-1),j) = al_jk(j,1:(j-1));
    b_1_32_jk(1:(j-1),j) = b_1_32_jk(j,1:(j-1));
    b_2_32_jk(1:(j-1),j) = b_2_32_jk(j,1:(j-1));
    D_b_0_12_jk(1:(j-1),j) = D_b_0_12_jk(j,1:(j-1));
    
    nj(j,:) = sqrt(mu/OEp(j,1)^3);
    EQp(j,:) = Kep_2_Equinoctial( OEp(j,:) );
    
end

% Assemble matrices
Bjk = .25 *mk_mcj .* nj.* al_jk.* al_bar_jk.* b_1_32_jk ;
Ajk =-.25 *mk_mcj .* nj.* al_jk.* al_bar_jk.* b_2_32_jk ;
ids = logical(eye(np));
Bjk(ids) =-sum(Bjk,2);
Ajk(ids) =-Bjk(ids);

[eji_b, gi] = eig(Ajk,'vector');
[Iji_b, fi] = eig(Bjk,'vector');

hj = EQp(:,2);
kj = EQp(:,3);
pj = EQp(:,4);
qj = EQp(:,5);

SisBi = eji_b\hj;
SicBi = eji_b\kj;

Tisyi = Iji_b\pj;
Ticyi = Iji_b\qj;

Si = sqrt( SisBi.^2 + SicBi.^2 );
Ti = sqrt( Tisyi.^2 + Ticyi.^2 );

Betai = atan2( SisBi, SicBi );
Gammai= atan2( Tisyi, Ticyi );

eji = (ones(np,1)*Si').*eji_b ;
Iji = (ones(np,1)*Ti').*Iji_b ;

%% Assemble matrices for sigma(t) computation
Ejk = zeros(np);
for j=1:np
    for k=1:np
        
        if k < j % Internal perturber
            Ejk(j,k) = GMp_GMs(k)*mu *D_b_0_12_jk(j,k) /OEp(j,1)^3/nj(j,k) ;            
        elseif k > j % External perturber
            Ejk(j,k) = -GMp_GMs(k)*mu *D_b_0_12_jk(j,k) /OEp(j,1)/nj(j,k)/OEp(k,1)^2 ;
        end
        
    end
end
Ejkv = sum(Ejk,2);


%% Setup Output Structure
secular_model.ev = eji;
secular_model.iv = Iji;
secular_model.f  = fi;
secular_model.g  = gi;
secular_model.phaseE = Betai;
secular_model.phaseI = Gammai;
secular_model.p  = i0;

secular_model.ev_n = eji_b;
secular_model.iv_n = Iji_b;

secular_model.Ejk = Ejkv;
secular_model.s0  = OEp(:,6);

secular_model.pl = pl_list;
