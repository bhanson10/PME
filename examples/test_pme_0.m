clear all; close all; clc; format shortG;

%% parameters 
d = 2; 
mu = zeros(d, 1); S = eye(d); 
[Z, W] = CUT6(mu, S);
[X, P, l] = pme_cut2pdf(Z, W); 

l,