function plot_distancevstime(n,cons,xt,tv1,d_pe,d_nbp,d_ll)
    F = figure(n);
%     clf
%     F.Position = [-1203 239 560 189];

    xsc = cons.yr;
    ysc = cons.AU; %cons.Re;

    plot(tv1/xsc, d_pe/ysc, 'b')
    grid on
    hold on
    xlabel('t (yr)')
    ylabel('d (au)')


    plot( xt, d_nbp/ysc, 'r' )

    plot( xt, d_ll/ysc, 'g')
    xlim(tv1([1 end])/xsc)
    legend('post-enc','n3BPert','sec-LL')

end