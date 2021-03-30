function plot_oe_evolution(set,t,OE,id,spec,tb,t_scale)

t = t/t_scale ;

F = figure(id);
%set(F, 'Position',[-1239.8 98.6 1199.2 634.4])
% F.Position = [-1239.8 98.6 1199.2 634.4];

% Fullscreen almost in center screen
% F.Position = [46.6 133.8 1260 628];

% Powerpoint format
% F.Position = [317 426.6 989.6 335.2];

% % Bigger format
% F.Position = [189.8 189.8 1168 526.4]; % Centered
% F.Position = [71.4 164.2 991.2 526.4]; % West

if set == 'oe'
    
    legends_oe = {'a (au)','e','i (deg)','\Omega (deg)','\omega (deg)','\sigma (deg)'};
    ylims      = [ OE(1,1)*max( (abs(OE(1,1)-OE(end,1))/OE(1,1))*[-1 1],[0.95 1.05]); -0.01 0.1; 0 pi; 0 2*pi; -2*pi 2*pi; -2*pi 2*pi ];
    for i=1:6
        subplot(2,3,i)
        
        if i>=3
            plot( t, wrapTo360(OE(:,i)),spec )
            hold on
        else
            plot( t, OE(:,i),spec )
            hold on
        end
        grid on
        ylabel( legends_oe{i} )
%         ylim( ylims(i,: ))
        xlim(tb)    
    end
    
    
elseif set == 'xyz'
    
    legends_xyz = {'x (km)','y (km)','z (km)','v_x (km/s)','v_y (km/s)','v_z (km/s)',};
    for i=1:6
        subplot(2,3,i)
        plot( t, OE(:,i) )
        grid on
        hold on
        ylabel( legends_xyz{i} )
        
    end

end

% suptitle(titola)
% suptitle([titola ' - [\Omega_0 = ' num2str(OE(1,4)*180/pi) ...
%     ', inc_0 = ' num2str(OE(1,3)*180/pi) ']'])