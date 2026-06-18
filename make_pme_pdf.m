function pdf = make_pme_pdf(l, E)

    pdf = @(X) exp(phi(X, l, E));

end


function y = phi(X, l, E)

    p = size(E,1);
    G = zeros(size(X,1), p);

    for k = 1:p
        G(:,k) = prod(X.^E(k,:), 2);
    end

    y = G * l;
end