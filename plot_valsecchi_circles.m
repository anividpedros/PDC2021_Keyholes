function plot_valsecchi_circles(n,circ,cons,U,kmax,xi,zeta)

    % Plot Valsecchi circles!
    % n = number of figure

    F  = figure(n);

    nr = size(circ,1);
    co = winter(22);
    thv= 0:(2*pi/99):2*pi ;
    sc = cons.Re;

    focus_factor = sqrt(1 + 2*cons.GMe/(cons.Re*U^2));
    RE_focussed  = focus_factor;

    for i=1:nr

        k = circ(i,1);
        D = circ(i,3);    
        R = circ(i,4);

        xi_circ   = R*cos(thv);
        zeta_circ = D + R*sin(thv);

        cc = co(k,:);
        plot( xi_circ/sc, zeta_circ/sc, 'Color', cc )
        hold on

    end
    colormap(co);
    fill(RE_focussed*cos(thv), RE_focussed*sin(thv),'white');
    plot(RE_focussed*cos(thv), RE_focussed*sin(thv),'k');
    plot(cos(thv), sin(thv),'k--');
    % plot(3*cos(thv), 3*sin(thv),'k--');

    grid on
    axis equal
    caxis([1 kmax])
    cb = colorbar;
    cb.Label.String = 'k';
    axis([-1 1 -1 1]*15)
    xlabel('\xi (R_\oplus)');
    ylabel('\zeta (R_\oplus)');

    plot( xi/sc, zeta/sc, '+k' )
end