clear all; close all; clc; format long;

%% parameters
load("colors.mat"); 
d = 2; 
t = 0:0.005:0.2; Nt = numel(t); 
ts = reshape(t,[1, 1, length(t)]); t = t'; 
mu0 = zeros(d, 1); S0 = eye(d); 

%% CUT
[Z0, W] = cut6(mu0, S0); N = size(Z0, 1);
Z0 = reshape(Z0,[N,d,1]);
Z = Z0 ./ (1 - Z0 .* ts);

%% reconstruction 
f = figure; clf; f.Position = [250, 250, 1500, 550];
tiledlayout(1, 3, "Padding", "compact"); 
nexttile(1); hold on; box on; axis equal; 
set(gca, 'FontName', 'Times','FontSize', 24, "LineWidth", 2); 
leg = legend(gca, 'Orientation', 'Horizontal', 'FontSize', 24, ...
             'FontName', 'times', 'Interpreter', 'latex', "LineWidth", 2);
leg.Layout.Tile = 'north'; 
xlabel("$x$", "Interpreter", "latex", "FontSize", 32);
ylabel("$y$", "Interpreter", "latex", "FontSize", 32);
leg = legend(gca, 'Orientation', 'Horizontal', 'FontSize', 24, ...
             'FontName', 'times', 'Interpreter', 'latex', "LineWidth", 2);
leg.Layout.Tile = 'south';

% parameters
mu0 = [0; 0]; S0 = eye(2); M = 4; 
[xgrid,ygrid] = meshgrid(linspace(-8,10,1000)); XY = [xgrid(:), ygrid(:)]; 

% truth, initial
P0 = quadGrowthPDF(xgrid, ygrid, 0, mu0, S0);
p.color = hanblue; p.fill = 1; p.name = "Truth"; 
plot_nongaussian_surface(XY, 'P', P0(:), 'p', p); clear p; 
title("$t=0$", 'FontSize', 28, 'FontName', 'times', 'Interpreter', 'latex');

% reconstruct, initial
mu_0_r = reconstruct_cut(squeeze(Z(:, :, 1)), W, 1, 6); % \mathcal{N}(\mu, \Sigma)
S_0_r = reconstruct_cut(squeeze(Z(:, :, 1)), W, 2, 6); % \mathcal{N}(\mu, \Sigma)
[X_0_r, P_0_r] = pme_cut2pdf(squeeze(Z(:, :, 1)), W, 'M', M); % PME (order M)
scatter(squeeze(Z(:, 1, 1)), squeeze(Z(:, 2, 1)), 250, "d", "filled", "MarkerFaceColor", hanred, "HandleVisibility", "off");  
p.color = hangreen; p.ls = "-"; p.lw = 3; 
plot_nongaussian_surface(X_0_r, 'P', P_0_r, 'p', p);
p.color = hanred; p.ls = "--"; 
plot_gaussian_ellipsoid(mu_0_r, S_0_r, 'p', p); 
scatter(nan, nan, 250, "d", "filled", "MarkerFaceColor", hanred, "DisplayName", "CUT6 $(N=$ " + num2str(N) + "$)$"); 
plot(nan, nan, "Color", hanred, "LineWidth", 3, "LineStyle", "--", "DisplayName", "$\mathcal{N}(\hat{x}, P)$"); 
plot(nan, nan, "Color", hangreen, "LineWidth", 3, "LineStyle", "-", "DisplayName", "PME $(M=$ " + num2str(M) + "$)$"); 
xlim([-4, 10]); ylim([-4, 10]);  

nexttile(2); hold on; box on; axis equal; 
set(gca, 'FontName', 'Times','FontSize', 24, "LineWidth", 2); 
xlabel("$x$", "Interpreter", "latex", "FontSize", 32);
ylabel("$y$", "Interpreter", "latex", "FontSize", 32);

% truth, middle
i = round(size(Z, 3) / 2); 
t_i = t(i); 
Pi = quadGrowthPDF(xgrid, ygrid, t_i, mu0, S0);
if any(isnan(Pi(:)))
    disp('Vector contains at least one NaN');
end
p.color = hanblue; p.fill = 1; 
plot_nongaussian_surface(XY, 'P', Pi(:), 'p', p); clear p;  
title("$t=$ " + num2str(t_i), 'FontSize', 28, 'FontName', 'times', 'Interpreter', 'latex');

% reconstruct, middle
[X_i_r, P_i_r] = pme_cut2pdf(squeeze(Z(:, :, i)), W, 'M', M); 
mu_i_r = reconstruct_cut(squeeze(Z(:, :, i)), W, 1, 6); 
S_i_r = reconstruct_cut(squeeze(Z(:, :, i)), W, 2, 6); 
scatter(squeeze(Z(:, 1, i)), squeeze(Z(:, 2, i)), 250, "d", "filled", "MarkerFaceColor", hanred, "HandleVisibility", "off");  
p.color = hangreen; p.ls = "-"; p.lw = 3; 
plot_nongaussian_surface(X_i_r, 'P', P_i_r, 'p', p);
p.color = hanred; p.ls = "--"; 
plot_gaussian_ellipsoid(mu_i_r, S_i_r, 'p', p); 
xlim([-4, 10]); ylim([-4, 10]);  

nexttile(3); hold on; box on; axis equal; 
set(gca, 'FontName', 'Times','FontSize', 24, "LineWidth", 2); 
xlabel("$x$", "Interpreter", "latex", "FontSize", 32);
ylabel("$y$", "Interpreter", "latex", "FontSize", 32);

% truth, final
Pf = quadGrowthPDF(xgrid, ygrid, t(end), mu0, S0);
p.color = hanblue; p.fill = 1; 
plot_nongaussian_surface(XY, 'P', Pf(:), 'p', p); clear p; 
title("$t=$ " + num2str(t(end)), 'FontSize', 28, 'FontName', 'times', 'Interpreter', 'latex');

% reconstruct, final
[X_f_r, P_f_r] = pme_cut2pdf(squeeze(Z(:, :, end)), W, 'M', M);
mu_f_r = reconstruct_cut(squeeze(Z(:, :, end)), W, 1, 6); 
S_f_r = reconstruct_cut(squeeze(Z(:, :, end)), W, 2, 6); 
scatter(squeeze(Z(:, 1, end)), squeeze(Z(:, 2, end)), 250, "d", "filled", "MarkerFaceColor", hanred, "HandleVisibility", "off");  
p.color = hangreen; p.ls = "-"; p.lw = 3; 
plot_nongaussian_surface(X_f_r, 'P', P_f_r, 'p', p);
p.color = hanred; p.ls = "--"; 
plot_gaussian_ellipsoid(mu_i_r, S_i_r, 'p', p); 
xlim([-4, 10]); ylim([-4, 10]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = quadGrowthPDF(x, y, t, mu, P)

    % Preimage
    denom_x = 1 + t*x;
    denom_y = 1 + t*y;

    x0 = x ./ denom_x;
    y0 = y ./ denom_y;

    % Jacobian
    J = 1 ./ (denom_x.^2 .* denom_y.^2);

    % Gaussian terms
    Pinv = inv(P);
    norm_const = 1 / (2*pi*sqrt(det(P)));

    dx = x0 - mu(1);
    dy = y0 - mu(2);

    exponent = -0.5 * ( ...
        Pinv(1,1)*dx.^2 + ...
        2*Pinv(1,2)*dx.*dy + ...
        Pinv(2,2)*dy.^2 );

    p0 = norm_const .* exp(exponent);

    % ---- FIX 1: remove Inf/NaN from exponent path
    badMap = ~isfinite(x0) | ~isfinite(y0) | ~isfinite(J);

    % enforce safe values so no Inf/NaN propagate
    p0(badMap) = 0;
    J(badMap)  = 0;

    % ---- final density
    p = p0 .* J;

    % ---- FIX 2 (extra safety)
    p(~isfinite(p)) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%