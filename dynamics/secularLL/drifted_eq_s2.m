function EQd = drifted_oe_at_secss( secular_model, dt )

% Using secular model generated by secular_model_10BP

% Structure set-up
% ev: Matrix of eigenvectors of H and K
% iv: Matrix of eigenvectors of P and Q
% g:  Frequencies of the H-K eigenvalue problem
% f:  Frequencies of the P-Q eigenvalue problem
% phaseE: Integration constants set to match ICs at t0
% phaseI: Integration constants set to match ICs at t0
% p: Position of the particle in the secular model 
%    (ordered by semi-major axis!)

eji = secular_model.ev;
Iji = secular_model.iv;
gi  = secular_model.g;
fi  = secular_model.f;
Betai = secular_model.phaseE;
Gammai= secular_model.phaseI;
Ev  = secular_model.Ejk ;
s0  = secular_model.s0  ;

i0    = secular_model.p;
np    = length( gi );
pl    = secular_model.pl;

% Evaluate secular solution
H = eji*sin( gi*dt + Betai );
K = eji*cos( gi*dt + Betai );
P = Iji*sin( fi*dt + Gammai );
Q = Iji*cos( fi*dt + Gammai );

S = s0 + Ev*dt ;

% Common for both
zcol = zeros(np,1);
EQd = [zcol H K P Q S];