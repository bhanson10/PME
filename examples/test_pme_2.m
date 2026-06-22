clear all; close all; clc; format long;

%% parameters
tic,
load("colors.mat"); M = 2; n = 25; d = 5; 
mu = randn(d, 1); A = randn(d, d); S = A*A';     
[Z, W] = cut6(mu, S);
[X, P, l] = pme_cut2pdf(Z, W, 'M', M, 'n', n); 

l,

figure("Position", [362,158,753,625]); hold on; axis equal; 
set(gca, 'FontName', 'times', 'FontSize', 18, "LineWidth", 2); box on; 
p.color = hangreen; p.fill = 1; p.type = "grid"; p.hist = 1; 
plot_corner_pdf(X, 'P', P, 'p', p);
p.color = hanred; p.fill = 0; p.type = "cukf"; p.marker = "d"; p.ms = 50; p.points_too = 1; 
plot_corner_pdf(Z, 'P', W, 'p', p); 

LH(1) = scatter(nan, nan, 100, 'd', 'MarkerFaceColor', hanred, 'MarkerEdgeColor', 'none');
L{1} = "CUT6";
LH(2) = plot(nan, nan,  "Color", hanred, "LineWidth", 2);
L{2} = "$\mathcal{N}(\mu, \Sigma)$";
LH(3) = fill(nan, nan, hangreen, 'FaceAlpha', 0.5, "EdgeColor", "none");
L{3} = "PME ($M=2$)";
leg = legend(LH, L, "FontSize", 18, "FontName", "times", "Interpreter", "latex",...
             "LineWidth", 1, "Orientation", "horizontal");
leg.Layout.Tile = 'south'; 
toc,