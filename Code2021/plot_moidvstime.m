function plot_moidvstime(n,cons,tv1,xt,MOID0,MOIDnbp,MOIDsec)


    F = figure(n);
%     clf
%     F.Position = [-1202 -54 560 189];

    xsc = cons.yr;
    ysc = cons.Re;

    plot( tv1([1 end])/xsc, MOID0*[1 1]/ysc, 'b--' )
    hold on
    plot( xt, MOIDnbp/ysc, 'r' )

    plot( xt, MOIDsec/ysc, 'g' )

    grid on
    xlabel('time (yr)')
    ylabel('MOID R_\oplus)')%(DU)');

    legend('post-enc','n3BPert','sec-LL')%,'4BP')
    xlim(tv1([1 end])/xsc)

end