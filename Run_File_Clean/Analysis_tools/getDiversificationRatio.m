function DR = getDiversificationRatio(x, Ret)
% This function computes the diversification ratio of a portfolio
% 
% INPUTS:
% x:        Vector of portfolio's weights
% Ret:      Log-Returns
%
% OUTPUT:
% DR        Diversification Ratio

% Volatilities of the returns
vola = std(Ret); 
% Covariance matrix of the returns
V = cov(Ret);
% Portfolio's volatility
volaPtf =sqrt(x'*V*x);

% Compute the Diversification Ratio
DR = (x'*vola')/volaPtf;

end