function Xdot = ShortPeriod_prop(t,semieq,cons,secular_model)
    

OE0 = Equi_2_Keplerian(semieq);
OEe = cons.OEp;

% Full perturbing potential component    
dR_d3BP = dR_3BP_dr(OE0,OEe,t,cons)';
dX_dOE  = dCart_dKep(OE0,t,cons) ;
dOE_dseq = dKep_dSEqui(semieq,OE0, cons);

% Partial of potential wrt. to semieq
% dR_dseq
Xd = (dR_d3BP*dX_dOE)*dOE_dseq ;

sc = 1/sqrt(cons.GMs*OE0(1)); % 1/L  or  1/na^2
Xdot_F = [ Xd(6) ; 
    Xd(2)*sc; -Xd(3)*sc;
    Xd(4)*sc; -Xd(4)*sc;
    -Xd(1)];

% Secular perturbing potential component
% Xdot_S = zeros(6,1);
Xdot_S = dX_sec_LL(secular_model, t );
% Fxn to dev!

Xdot = Xdot_F + Xdot_S ;

% Fxns to dev list:
% X- Equi to Keplerian (Recycled)
% X- dKep_dSemieq
% X- dsecular_dt

end

function X = dKep_dSEqui(x,oe, cons)
% Partial of Keplerian with respect to Semi-equinoctial elements
% Inputs: x: Semi-equinoctial elements [L,h,k,p,q,sigma]
%         oe: Keplerian elements [a,e,i,Omega,omega,sigma]
mu = cons.GMs;
X  = zeros(6);

X(1,1) = 2*x(1)/mu ;
X(2,2:3) = [x(2)/oe(2) x(3)/oe(2)];
X(3,4:5) = [x(4)/oe(3) x(5)/oe(3)];

i2 = oe(3)*oe(3);
e2 = oe(2)*oe(2);
X(4,4:5) = [x(5)/i2 -x(4)/i2];
X(5,2:5) = [x(3)/e2 -x(2)/e2 -x(5)/i2 +x(4)/i2];

X(6,6) = 1;

end


function X = dX_sec_LL(secular_model, dt )

    Ajk = secular_model.Ajk;
    Bjk = secular_model.Bjk;
    Ej  = secular_model.Ejk;
    
    EQd = drifted_eq_s2( secular_model, dt );
    
    X = [0;0;0; 0;0;0];
    X(1) = 0;
    X(2) =  Ajk(1,:)*EQd(:,3);
    X(3) = -Ajk(1,:)*EQd(:,2);
    X(4) =  Bjk(1,:)*EQd(:,5);
    X(5) = -Bjk(1,:)*EQd(:,4);
    X(6) = Ej(1) ;

end