function e = getEntropy(w, cov_matrix)
% This function computes the entropy in asset volatility of a vector of
% weights.
%
% INPUTS:
% w:            Vector of weights of the portfolio
% cov_matrix:   Covariance matrix of the returns
%
% OUTPUT:
% e:            Entropy of the portfolio

% Initialize the Entropy
e = 0;
% Extract the variances vector
variances = diag(cov_matrix);

% Compute the Entropy
for i=1:length(w)
    if w(i)>0
        e = e + sum( w(i)^2 * variances(i) / sum(w.^2 .* variances ) .* ...
                        log( w(i)^2 * variances(i) / sum(w.^2 .* variances) ) );
    end
end
% Invert the sign
e = -e;
end