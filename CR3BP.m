function Xdot = CR3BP(t,X,mu)

    x = X(1);
    y = X(2);
    z = X(3);
    dx = X(4);
    dy = X(5);
    dz = X(6);

    partialx = x - (1-mu).*((x+mu).^2+y.^2+z.^2)^(-3/2).*(x+mu) - mu.*((x-1+mu).^2+y.^2+z.^2)^(-3/2).*(x-1+mu);
    partialy = y - (1-mu).*((x+mu).^2+y.^2+z.^2)^(-3/2).*y - mu.*((x-1+mu).^2+y.^2+z.^2)^(-3/2).*y;

    ax = 2.*dy + partialx;
    ay = -2.*dx + partialy;
    az = -(1-mu).*z.*((x+mu).^2+y.^2+z.^2)^(-3/2) - mu.*z.*((x-1+mu).^2+y.^2+z.^2)^(-3/2);

    Xdot = [dx; dy; dz; ax; ay; az];

end