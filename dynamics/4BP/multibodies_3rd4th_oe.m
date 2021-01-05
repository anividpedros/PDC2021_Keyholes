function xdot = multibodies_3rd4th_oe(t,X,c)

% Acceleration on particle given by:
% central gravity + third body + fourth body
% the perturbers are modelled by a constant set of elements
% uses function cspice_conics.

% Constants
mu  = c.GMcentral; 
mus = c.GMthird; 
muj = c.GMfourth;

OEp = c.OEp;
OEp2= c.OEp2;

% Coordinates
r = norm(X(1:3));
r2 = r*r;
r3 = r2*r;

%======== TWO BODY EQUATIONS ==========
xdot = [X(4:6);
    -X(1:3)*mu/r3];


%======== OTHER PERTURBATIONS ==========
% 3rd Body
X_s = cspice_conics(OEp, c.JD0_p+t );

R_s = X_s(1:3);
R_i_sc  = R_s - X(1:3); % R_Sun/Spacecraft
R_ri1 = sqrt( sum( R_i_sc.^2 ) );
R_ri2 = R_ri1*R_ri1;
R_ri3 = R_ri2*R_ri1;
R = norm(R_s);
R3 = R*R*R;

a_3b_s = mus *( R_i_sc/R_ri3 - R_s/R3 );

% 4th Body
X_s = cspice_conics(OEp2, c.JD0_p+t );

R_s = X_s(1:3);
R_i_sc  = R_s - X(1:3); % R_Sun/Spacecraft
R_ri1 = sqrt( sum( R_i_sc.^2 ) );
R_ri2 = R_ri1*R_ri1;
R_ri3 = R_ri2*R_ri1;
R = norm(R_s);
R3 = R*R*R;

a_4b_s = muj *( R_i_sc/R_ri3 - R_s/R3 );

%======================================
% Sum of forces
xdot(4:6) = xdot(4:6) + a_3b_s + a_4b_s ;
%======================================


end