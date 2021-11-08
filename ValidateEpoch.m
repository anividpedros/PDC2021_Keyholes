% SOLUTION OF A PLANETARY ENCOUNTER: PDC 2021
% 7th IAA Planetary Defense Conference 2021
clc
clear
close all
format shortG

filename = 'ValidateEpoch.m';
filepath = matlab.desktop.editor.getActiveFilename;
currPath = filepath(1:(end-length(filename)));
cd(currPath)
addpath(genpath(pwd))

% addpath(genpath/('/Users/anividpedrosfaura/PDC2021_Keyholes/MOIDroutines/MOID_MATLAB'));

%% ===== Description ===== 
% Original script: main_1flyby_pdc2021.m
% Extension: Validates convergence in heliocentric
% orbit elements for a chosen epoch

%% Script steps
% 1. Read heliocentric elements (spice file) - date close to the encounter
% 2. Modify for circular Earth and adjust NEO orbit for encounter still happening
% 3. Switch to geocentric Keplerian elements and find time of perigeetp
% 4. Coordinates of the modified target plane(ξp,ζp,vp)
% 5. Compute the corresponding target plane pre-encounter coordinates (input for Öpik).(ξ,η,ζ)
% 6. Use Öpik theory to compute post-encounter coordinates
% 7. Convert to Heliocentric: Use Öpik theory formulae

%% 1. Read heliocentric elements
cspice_furnsh( 'SPICEfiles/naif0012.tls.pc' )
cspice_furnsh( 'SPICEfiles/gm_de431.tpc' )
cspice_furnsh( 'SPICEfiles/pck00010.tpc' )
cspice_furnsh( 'SPICEfiles/de431_part-1.bsp' )
cspice_furnsh( 'SPICEfiles/de431_part-2.bsp' )

% 2021 PDC
cspice_furnsh( 'SPICEfiles/2021_PDC-s11-merged-DE431.bsp' )
ast_id= '-937014';

cons.AU  = cspice_convrt(1, 'AU', 'KM');
cons.GMs = cspice_bodvrd( 'SUN', 'GM', 10 );
cons.GMe = 398600.43543609593;
cons.Re  = 6378.140;

cons.yr  = 365.25 * 24 * 3600 ;
cons.Day = 3600*24; 

epoch0 = '2021 October 20 TDB';
et0 = cspice_str2et( epoch0 ) + 24*.715*3600;

N = 20;
kep_ast_all = zeros(8,N);

for i = 1:N
    
    et = et0 - 24*3600*(i-1);
    
    state_ast = cspice_spkezr( ast_id, et, 'ECLIPJ2000', 'NONE', '10' );
    kep_ast = cspice_oscelt( state_ast, et, cons.GMs );
    sma_ast = kep_ast(1)/(1-kep_ast(2));
    kep_ast(1) = sma_ast;
    
    kep_ast_all(:,i) = kep_ast;
       
end

%%
figure
subplot(5,1,1)
plot(20:-1:(21-N),kep_ast_all(1,:),'LineWidth',1.5)
grid on
xlabel('Date in October')
ylabel('sma [km]')
title('Asteroid Heliocentric Orbit Elements')
subplot(5,1,2)
plot(20:-1:(21-N),kep_ast_all(2,:),'LineWidth',1.5)
grid on
xlabel('Date in October')
ylabel('e')
subplot(5,1,3)
plot(20:-1:(21-N),kep_ast_all(3,:),'LineWidth',1.5)
grid on
xlabel('Date in October')
ylabel('inc [rad]')
subplot(5,1,4)
plot(20:-1:(21-N),kep_ast_all(4,:),'LineWidth',1.5)
grid on
xlabel('Date in October')
ylabel('RAAN [rad]')
subplot(5,1,5)
plot(20:-1:(21-N),kep_ast_all(5,:),'LineWidth',1.5)
grid on
xlabel('Date in October')
ylabel('\omega [rad]')



