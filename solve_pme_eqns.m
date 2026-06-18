function lambda = solve_pme_eqns(M, E, X, Omega, opts)
% SOLVE_PME Maximum-entropy moment matching
%
% lambda = solve_pme(M, E, X, Omega)
%
% Inputs
% -------
% M      : p x 1 moment vector
% E      : p x d exponent matrix
%          E(i,:) = exponents for basis function i
% X      : N x d quadrature points
% Omega  : d x 2 domain bounds
%
% Output
% -------
% lambda : p x 1 MaxEnt coefficients
%
% The basis functions are assumed to be monomials
%
%   g_i(x) = prod_j x_j^(E(i,j))
%
% and the density is
%
%   p(x) = exp(lambda' g(x)).
%
% The dual objective is
%
%   Psi(lambda)
%      = w * sum(exp(Phi*lambda))
%      - lambda'*M
%
% which is strictly convex.

    if nargin < 5 || isempty(opts)
        opts = optimoptions('fminunc', ...
            'Algorithm','trust-region', ...
            'SpecifyObjectiveGradient',true, ...
            'HessianFcn','objective', ...
            'Display','iter', ...
            'OptimalityTolerance',1e-8, ...
            'StepTolerance',1e-10, ...
            'MaxFunctionEvaluations',1e7);
    end
    
    [N,d] = size(X);
    p     = size(E,1);
    
    vol = prod(Omega(:,2)-Omega(:,1));
    w   = vol/N;
    
    %% ------------------------------------------------------------
    % Build basis matrix
    % Phi(n,i) = prod_j X(n,j)^E(i,j)
    Phi = ones(N,p);
    
    for i = 1:p
        for j = 1:d
    
            e = E(i,j);
    
            if e ~= 0
                Phi(:,i) = Phi(:,i).*X(:,j).^e;
            end
    
        end
    end
    
    %% ------------------------------------------------------------
    % Initial guess
    
    lambda0 = zeros(p,1);
    
    %% ------------------------------------------------------------
    % Solve
    
    lambda = fminunc( ...
        @(l) dual_objective(l,Phi,w,M), ...
        lambda0, ...
        opts);

end

function [f,g,H] = dual_objective(lambda,Phi,w,M)
    s = Phi*lambda;
    
    % Log-sum-exp stabilization
    smax = max(s);
    e = exp(s-smax);
    Z = sum(e);
    f = w*exp(smax)*Z - lambda'*M;
    
    if nargout > 1
        e = exp(s);
        g = w*(Phi'*e) - M;
    end
    if nargout > 2
        e = exp(s);
        H = Phi'*(w*e.*Phi);
    end
end
