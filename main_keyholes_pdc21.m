% GENERATION OF MODIFIED KEYHOLES ON THE B-PLANE
% 7th IAA Planetary Defense Conference 2021
clc
clear

%% Algorithm steps
% 1. Extract Earth and Asteroid trajectories
% 2. States of Earth and Asteroid
% 3. Coordinates in Earth-Centered Orbital Frame
% 4. Opik Variables at CA
% 5. Scan for resonances (Compute Valsecchi circles)

%% 0. Constants
cspice_furnsh('mice-kernels\naif0012.tls.pc')
cspice_furnsh('mice-kernels\gm_de431.tpc')
cspice_furnsh('mice-kernels\pck00010.tpc')
cspice_furnsh('C:\Users\Oscar\Documents\Spice-Kernels\de431_part-1.bsp')
cspice_furnsh('C:\Users\Oscar\Documents\Spice-Kernels\de431_part-2.bsp')



%% 1. 
