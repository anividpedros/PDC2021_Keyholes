function F = plotBplane( phi, theta, R, idfig, c )

N = 10;

cp = cos(phi);   sp = sin(phi);
ct = cos(theta); st = sin(theta);

Rphi = rotationM(-phi,2);
Rth  = rotationM(-theta,1);
ROT  = Rth*Rphi ;
ROTt = ROT';

Rphi = rotationM(phi,2);
Rth  = rotationM(theta,1);
ROTt = Rphi*Rth ;

auxR = [cp 0 -sp;st*sp ct st*cp;ct*sp -st ct*cp];
ROTt = auxR';

x = (0:1:N)*2*R/N - R;
z = x;

for i=1:(N+1)
    for j=1:(N+1)
        
        rb = [x(i); 0; z(j)];
        ro = ROTt*rb;
        
        Xo(i,j) = ro(1);
        Yo(i,j) = ro(2);
        Zo(i,j) = ro(3);
        
    end
end

F = figure(idfig);

hold on
surf( Xo, Yo, Zo, 'FaceAlpha', .5, ...
    'EdgeColor',c,'FaceColor', c );
quiver3(0,0,0, R/2*st*sp, R/2*ct, R/2*st*cp,'Color',c,...
    'LineWidth',2,'MaxHeadSize',2)