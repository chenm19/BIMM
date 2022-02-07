%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script of Genetic Algorithm optimizations for Multiscale modeling of 
% tumorigenesis regulatory network. 
% 27 parameters need to be optimized
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

close all; clc; clear all; format short;
nvars = 27; %variable

%Upper and lower bounds 
Ub = ones(1,nvars) * 10;
Lb = ones(1,nvars) * 0.01;

opts = optimoptions('ga','PlotFcn', 'gaplotbestf', 'MaxGenerations', 50, ...
    'PopulationSize', 300);

%imply ga algorithm
[value, error] = ga(@fitness,nvars,[],[],[],[],Lb,Ub,[],[],opts); 

save("cancer_ga_outputs.mat");


