function [relRC, RC, mVol] = getRiskContributions(x, Ret)
% This function computes the relative risk contributions of a portfolio
%
% INPUTS:
% x:        Portfolio's weights
% Ret:      Log-Returns
%
% OUTPUTS:
% relRC:    Relative risk contributions vector
% RC        Absolute risk contributions vector
% mVol      Marginal risk conrtibutions vector

% Covariance matrix of the returns
V = cov(Ret);
% Portfolio's volatility
VolaPtf = sqrt(x'*V*x);

% Marginal risk conrtibutions
mVol = V*x/VolaPtf;
% Absolute risk contributions
RC = mVol.*x;
% Relative risk contributions
relRC = RC/sum(RC);

end