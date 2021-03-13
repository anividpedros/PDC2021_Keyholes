%% MOID test - SCRIPT
clc
clear
close all
format shortG

filename = 'MOIDtest.m';
filepath = matlab.desktop.editor.getActiveFilename;
currPath = filepath(1:(end-length(filename)));
cd(currPath)
addpath(genpath(pwd))


A.sma = 1.00000011;
A.e = 0.01671022;
A.argp = deg2rad(102.94719 );
A.Omega = deg2rad(-11.26064);
A.i = deg2rad(0.00005);

B.sma = 1.5;
B.e = 0.3;
B.argp = deg2rad(140);
B.Omega = deg2rad(-35);
B.i = deg2rad(10);

[moid] = ComputeMOID(A,B)