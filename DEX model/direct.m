%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script of direct optimizations for Multiscale modeling of tumorigenesis
% regulatory network. 
% 27 parameters need to be optimized
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
close all; clc; clear all; format short;
nvars = 27; %variable


%optimization
 Inits = rand(1,nvars)*5;  

Ub = ones(1,nvars) * 10;
Lb = ones(1,nvars) * 0.01;

options = optimoptions('patternsearch','TolFun',1e-10,'TolX',1e-8,...
    'TolMesh',1e-8,'TolCon',1e-10,'Vectorized','off',...
    'MaxIter', 800,'Display','iter','PlotFcn',{@psplotbestf,@psplotbestx});


%imply direct optimization algorithm
[value, error] = patternsearch(@fitness,Inits,[],[],[],[],Lb,Ub,[],options); 
save("cancer_drct.mat");



