function [Mi, E] = cut2mom(Z, W, M)
% CUT2MOM Conjugate Unscented Transform (CUT) to finite scalar moments (M) 
%
%   M = cut2mom(Z, W, k) returns the k-by-1 vector of expected scalar moment
%   values defined by the CUT sigma points Z and weights W, with p being
%   the number of scalar moments fit
% 
%   This function assumes that the monomial basis functions are MATLAB syms
%
% Example:
%   mu = [0; 0]; S = [1 0; 0 1]; 
%   [Z, W] = CUT4(mu, S);
%   M = 2; 
%   M = cut2mom(Z, W, M);
%
% Copyright 2026 by Benjamin L. Hanson, published under BSD 2-Clause License.
    
    [N, d] = size(Z); 
    p = nchoosek(d + M, d); 
    Mi = zeros(p, 1); 
    E = zeros(p, d);

    for i = 1:p
        E(i, :) = gen_E(d, M, i); 
        for j = 1:N
            Mi(i) = Mi(i) + W(j) * prod(Z(j, :) .^ E(i, :));
        end
    end
end


function E = gen_E(d, M, idx)
    E = zeros(d, 1); 
    remaining_idx = idx; 
    for k = 0:M
        block_size = comp_count(d, k);
        if(remaining_idx > block_size)
            remaining_idx = remaining_idx - block_size;
            continue;
        end

        rem_dim = d; 
        rem_sum = k; 
        for i = 1:d-1
            for val = rem_sum:-1:0
                count = comp_count(rem_dim - 1, rem_sum - val);
                if(remaining_idx > count)
                    remaining_idx = remaining_idx - count; 
                else
                    E(i) = val; 
                    rem_sum = rem_sum - val; 
                    rem_dim = rem_dim - 1; 
                    break;
                end
            end
        end
        E(d) = rem_sum; 
        return; 
    end
end

function c = comp_count(r,k)
    if r == 0
        c = double(k == 0);
        return
    end
    c = nchoosek(r + k - 1, k);
end