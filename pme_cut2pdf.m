function varargout = pme_cut2pdf(Z, W, varargin)
% pme_cut2pdf.m
% Benjamin Hanson, 2026
% 
%   Given the CUT sigma points Z and weights W, this function constructs a
%   a PDF based on the principle of maximum entropy (PME)
% 
% Inputs:
%          Z -- N-by-D sigma points
%          W -- N-by-1 weights
%   varargin -- optional arguments
%               * M -- order of PME reconstruction
%               * X -- support of PDF
%               * opts -- fminunc options
%
%   [X, P] = pme_cut2pdf(...) returns just the PDF P over the space X
%
%   [X, P, l] = pme_cut2pdf(...) returns the PDF P over the space X and
%   the p-Lagrange multipliers l
% 
% Example:
%   load("colors.mat"); 
%   mu = [1; 2]; S = [1 0.4; 0.4 1]; M = 2; 
%   [Z, W] = CUT6(mu, S);
%   [X, P] = pme_cut2pdf(Z, W, 'M', M); 
%   figure; hold on; axis equal; 
%   p.color = hangreen; p.fill = 1; p.name = "Maximum Entropy"; 
%   plot_nongaussian_surface(X, 'P', P, 'p', p); 
%   scatter(Z(:, 1), Z(:, 2), 100, 'filled', "DisplayName", "CUT");
%   p.color = hanred; p.fill = 0; p.name = "$\mathcal{N}(\mu, \Sigma)$";  
%   plot_gaussian_ellipsoid(mu, S, 'sd', 3, 'p', p);
%   legend("Location", "north", "Orientation", "horizontal", "FontName", "Times", ...
%       "FontSize", 12, "Interpreter", "latex", "LineWidth", 1); 
%   drawnow;
%
% Copyright 2026 by Benjamin L. Hanson, published under BSD 2-Clause License.

    for i=1:2:length(varargin)
        if strcmp('M',varargin{i})
            M = varargin{i+1};
        elseif strcmp('X',varargin{i})
            X = varargin{i+1};
        elseif strcmp('n',varargin{i})
            n = varargin{i+1};
        elseif strcmp('opts',varargin{i})
            opts = varargin{i+1};
        else
            error(append("Unspecified argument: ", varargin{i}));
        end
    end
    
    if ~exist("M", "var")
        M = 2; 
    end
    if ~exist("n", "var")
        n = 100; 
    end
        d = size(Z, 2); 
    if ~exist("X", "var")
        Omega = zeros(d, 2);
        for i = 1:d
            Omega(i, 1) = min(Z(:, i)) - 0.1 * (max(Z(:, i)) - min(Z(:, i))); 
            Omega(i, 2) = max(Z(:, i)) + 0.1 * (max(Z(:, i)) - min(Z(:, i))); 
        end  
        axes = cell(d,1);
        for j = 1:d
            axes{j} = linspace(Omega(j,1), Omega(j,2), n);
        end
        [grids{1:d}] = ndgrid(axes{:});
        X = zeros(numel(grids{1}), d);
        for j = 1:d
            X(:,j) = grids{j}(:);
        end
    else 
        Omega = zeros(d, 2);
        for i = 1:d
            Omega(i, 1) = min(X(:, i)) - 0.1 * (max(X(:, i)) - min(X(:, i))); 
            Omega(i, 2) = max(X(:, i)) + 0.1 * (max(X(:, i)) - min(X(:, i))); 
        end
    end
    if ~exist("opts", "var")
        opts = optimoptions('fminunc', ...
            'Algorithm','trust-region', ...
            'SpecifyObjectiveGradient',true, ...
            'HessianFcn','objective', ...
            'Display','iter', ...
            'OptimalityTolerance',1e-8, ...
            'StepTolerance',1e-10, ...
            'MaxFunctionEvaluations',1e7);
    end
    
    [Mi, E] = cut2mom(Z, W, M),
    l = solve_pme_eqns(Mi, E, X, Omega, opts);
    pdf = make_pme_pdf(l, E);

    P = pdf(X); P = P ./ sum(P);

    varargout{1} = X;
    varargout{2} = P;
    if nargout > 2
        varargout{3} = l;
    end
end