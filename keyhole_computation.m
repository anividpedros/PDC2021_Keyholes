function [circles, kh_good] = keyhole_computation(n,circ,cons,U,DU,U_nd,theta,phi)
    RE_au = cons.Re/DU;
    
    co = winter(22);
    thv= 0:(2*pi/99):2*pi ;
    focus_factor = sqrt(1 + 2*cons.GMe/(cons.Re*U^2));

    RE_focussed  = focus_factor;
    m  = cons.GMe/cons.GMs ;
    circles = circ;

    F = figure(n);
    hold on

    sc = cons.Re/DU;
    nr = size(circles,1);
    kh_good = [];
    
    for i=1:nr

        % New circles
        k = circles(i,1);
        h = circles(i,2);
        D = circles(i,3)/cons.Re;    
        R = circles(i,4)/cons.Re;    

        [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
            two_keyholes(k, h, D, R, U_nd, theta, phi, m,0,DU);

        cc = co(k,:);    
        plot(kh_down_xi(:,1)/sc,kh_down_zeta(:,1)/sc,kh_down_xi(:,2)/sc,kh_down_zeta(:,2)/sc,...
            'Color',cc);
        plot(kh_up_xi(:,1)/sc,kh_up_zeta(:,1)/sc,kh_up_xi(:,2)/sc,kh_up_zeta(:,2)/sc,...
             'Color',cc);

        % Register keyholes with solutions
    %     arcexist = sum(~isnan(kh_up_xi)) + sum(~isnan(kh_down_xi));
        R1 = sqrt(sum(kh_down_xi.^2 + kh_down_zeta.^2,2));
        R2 = sqrt(sum(kh_up_xi.^2 + kh_up_zeta.^2,2));
        arcexist = sum( (R1-3*RE_au*focus_factor)>0 ) + sum( (R2-3*RE_au*focus_factor)>0 );

    %     arcexist = sum( kh_up_zeta > 2*RE_au*focus_factor );
    %     arcexist = sum( kh_down_zeta < -5*RE_au );
    %     arcexist = sum( kh_up_zeta > 3.6*RE_au );

        if arcexist
            kh_good = [kh_good; i];
            fprintf('Keyhole num %g exists\n',i)
        end


    end
    colormap(co);
    fill(RE_focussed*cos(thv), RE_focussed*sin(thv),'white');
    plot(RE_focussed*cos(thv), RE_focussed*sin(thv),'k');
    plot(cos(thv), sin(thv),'k--');

    grid on
    axis equal
    caxis([1 20])
    cb = colorbar;
    cb.Label.String = 'k';
    axis([-1 1 -1 1]*10)
    xlabel('\xi (R_\oplus)');
    ylabel('\zeta (R_\oplus)');
end