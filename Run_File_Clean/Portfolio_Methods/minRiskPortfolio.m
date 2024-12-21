function Output_struct = minRiskPortfolio(Ptf, pwgt, pf_risk_Ptf, name_ptf)
% Compute the Minimum Variance Portfolio of the frontier.
% 
% INPUTS
% Ptf:              Portfolio object
% pwgt:             The weights of the portfolios on the frontier
% pf_risk_Ptf:      The volatility of the portfolios on the frontier
% name_ptf:         Name of the portfolio
%
% OUTPUTS
% Output_struct:    A struct containing the volatility, weights, return, 
%                   Sharpe ratio, name and object protfolio of the
%                   minimum variance portfolio

% find minimum risk portfolio (Minimum Variance Portfolio)
minRisk_Ptf = min(pf_risk_Ptf);                         % Volatility
minRiskWgt_Ptf = pwgt(:, pf_risk_Ptf == minRisk_Ptf);   % Weights
minRiskRet_Ptf = Ptf.AssetMean' * minRiskWgt_Ptf;       % Return
minRiskSR_Ptf = (minRiskRet_Ptf - Ptf.RiskFreeRate) / minRisk_Ptf; % Sharpe Ratio

% Display the minimum risk portfolio
print_portfolio(minRiskWgt_Ptf, Ptf.AssetList, minRiskRet_Ptf, minRisk_Ptf, minRiskSR_Ptf, name_ptf);

% Build a struct for the output
Output_struct = struct('Volatility', minRisk_Ptf, 'Weights', minRiskWgt_Ptf,...
    'Return', minRiskRet_Ptf, 'Sharpe_Ratio', minRiskSR_Ptf);
Output_struct.Name = name_ptf;
Output_struct.Ptf = Ptf;

end