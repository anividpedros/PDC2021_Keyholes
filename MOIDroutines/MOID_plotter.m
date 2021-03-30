function [F, D, V, v] = MOID_plotter( OE1, OE2, MA, inp )

% MOIDi = MOID_SDG( OE1([1 2 4 3 5]), OE2([1 2 4 3 5]) );

v = -180:2:180 ; v = v*pi/180 ;
V = -180:2:180 ; V = V*pi/180 ;

nv = length(v);
nV = length(V);

if MA % Map in MA space
    
    D = zeros(nV,nv);
    for i=1:nV
        ri = coordinates_AnAstro(0, [OE1(1:5) V(i)], 0,0,1,1);
        for j=1:nv
            rj = coordinates_AnAstro(0, [OE2(1:5) v(j)], 0,0,1,1);
            D(i,j) = norm( ri(1:3) - rj(1:3) );
        end
    end
    
    
    else % Map in TA space
    
    D = zeros(nV,nv);
    for i=1:nV
        ri = coordinates_AnAstro(0, [OE1(1:5) 0 V(i)], 0,0,1,1);
        for j=1:nv
            rj = coordinates_AnAstro(0, [OE2(1:5) 0 v(j)], 0,0,1,1);
            D(i,j) = norm( ri(1:3) - rj(1:3) );
        end
    end
    
end

%     Dmin = min(min(D)) ;
%     Dmax = max(max(D)) ;

Dmin = 1e-4;
Dmax = 1.5;

% lv = Dmin:(Dmax-Dmin)/100:Dmax;
lv = logspace( log10(Dmin),log10(Dmax),50);

if length(inp) == 1
    F = figure(inp);
    %         clf
    
    % F.Position = [100.2 238.6 726.4 476];
    % F.Position = [-560 120 501.6 281.6];
    %         F.Position = [-483 104 424.6 297.6];
    
else
    F = figure(inp(1));
    subplot(inp(2),inp(3),inp(4))
end

contour(V*180/pi,v*180/pi,D',lv)
%     colorbar
hold on

% xlabel('M (Mean anomaly)')
% ylabel('m (mean anomaly)')
xlabel('V (True anomaly)')
ylabel('v (true anomaly)')
xticks(-180:45:180);
yticks(-180:45:180);

axis equal
axis([-180 180 -180 180])

% title(['MOID = ' num2str(MOIDi) ' (au)'])
% title(['MOID = ' num2str(MOIDi) ' (au)'])
    
end