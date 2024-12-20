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

% e = 0;
% for i=1:length(w)
%     if(w(i)<1e-6)
% 
%     else
%         e = e + w(i)*log(w(i));
%     end
% end
% e = -e;

e = 0;

variances = diag(cov_matrix);

for i=1:length(w)
    if w(i)>0
        e = e + sum( w(i)^2 * variances(i) / sum(w.^2 .* variances ) .* ...
                        log( w(i)^2 * variances(i) / sum(w.^2 .* variances) ) );
    end
end
%e = sum(w.^2 .* diag(cov_matrix)/ sum(w.^2 .* diag(cov_matrix)) .* ...
%                        log(w.^2 .* diag(cov_matrix)/ sum(w.^2 .* diag(cov_matrix))));

e = -e;
end