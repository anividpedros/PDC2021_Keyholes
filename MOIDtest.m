%% MOID test - SCRIPT
clc
clear
close all
format shortG

filename = 'main_1flyby_pdc2021.m';
filepath = matlab.desktop.editor.getActiveFilename;
currPath = filepath(1:(end-length(filename)));
cd(currPath)
addpath(genpath(pwd))
