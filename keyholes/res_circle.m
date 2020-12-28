function [xi,zeta] = res_circle(alpha,D,R)

xi = R.*cos(alpha);
zeta = D + R.*sin(alpha);

end