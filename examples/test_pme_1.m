clear all; close all; clc; format shortG;

%% parameters
load("colors.mat"); 
mu = [1; 2]; S = [1 0.4; 0.4 1]; M = 2; 
[Z, W] = CUT6(mu, S); 
[X, P, l] = pme_cut2pdf(Z, W, 'M', M); 

l,

figure("Position", [362,158,753,625]); 
tiledlayout(1, 1, "TileSpacing", "compact"); nexttile(1); hold on; axis equal; 
set(gca, 'FontName', 'times', 'FontSize', 24, "LineWidth", 2); box on; 
p.color = hangreen; p.fill = 1; p.name = "PME $(M=2)$"; 
plot_nongaussian_surface(X, 'P', P, 'p', p); 
scatter(Z(:, 1), Z(:, 2), 250, 'filled', 'd', 'filled', 'MarkerFaceColor', hanred, "DisplayName", "CUT6 $(N=13)$");
p.color = hanred; p.fill = 0; p.name = "$\mathcal{N}(\mu, \Sigma)$";  
plot_gaussian_ellipsoid(mu, S, 'sd', [1, 2, 3], 'p', p);
xlim([-3, 5]); 
ylim([-2, 6]);  
leg = legend("FontSize", 22, "FontName", "times", "Interpreter", "latex",...
             "LineWidth", 1, "Orientation", "horizontal");
leg.Layout.Tile = 'south'; 

drawnow;