function Vrot = ax_angle(V,theta,e)
% Rotate the vector "V" by the angle "theta" counter-clockwise around the
% axis "e". The latter is a unit vector.

sth = sin(theta); cth = cos(theta);
Vdote = dot(V,e);
ecrossV = cross(e,V);

Vrot = V*cth + ecrossV*sth + (1 - cth)*Vdote*e;

end