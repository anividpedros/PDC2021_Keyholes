%% "Water Procedure"
% To avoid missing minima .
% 4 positions of the meridional plane, evenly distributed along the 
% inclined orbit. We move these points along as water droplets.

%Special care when only 1 minima found

function [vtrueB, vL, vdis, N] = WaterProcedure(incliB, omegaB, argpB, N, vtrueB, vL, vdis)
if N<2 
    vtrueB = zeros (4,1);
    vL = zeros (4,1);
    vdis = zeros (4,1);
    N = 4;
    for j = 1:4
        vtrueB(j) = (0.25+0.5*j)*pi; %evenly distributed 0.75 - 1.25 - 1.75 - 2.25
        a = (cos(omegaB)*cos(argpB+vtrueB(j)) - sin(omegaB)*sin(argpB+vtrueB(j))*cos(incliB));
        b = (sin(omegaB)*cos(argpB+vtrueB(j)) + cos(omegaB)*sin(argpB+vtrueB(j))*cos(incliB));
        vL(j) = atan2(b,a);
        vdis(j) = 1e6; %set to something large
    end
end
end