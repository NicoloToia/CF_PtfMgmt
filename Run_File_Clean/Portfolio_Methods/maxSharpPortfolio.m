function Output_struct = maxSharpPortfolio(Ptf, name_ptf)
% Compute the Maximum Sharpe Ratio Portfolio
%
% INPUTS
% Ptf:          Portfolio object
% name_ptf:     Name of the portfolio
%
% OUTPUTS
% Output_struct: A struct containing the volatility, weights, return, 
%                Sharpe ratio, name and object protfolio of the maximum
%                Sharpe ratio portfolio

% Find maximum Sharpe ratio portfolio (Maximum Sharpe Ratio Portfolio)
maxSharpeWgt_Ptf = estimateMaxSharpeRatio(Ptf);         % Weights
maxSharpeRet_Ptf = Ptf.AssetMean' * maxSharpeWgt_Ptf;   % Return
maxSharpeRisk_Ptf = sqrt(maxSharpeWgt_Ptf' * Ptf.AssetCovar * maxSharpeWgt_Ptf);    % Volatility
maxSharpeSR_Ptf = (maxSharpeRet_Ptf - Ptf.RiskFreeRate) / maxSharpeRisk_Ptf;        % Sharpe Ratio

% Display the maximum Sharpe ratio portfolio
print_portfolio(maxSharpeWgt_Ptf, Ptf.AssetList, maxSharpeRet_Ptf, maxSharpeRisk_Ptf, maxSharpeSR_Ptf, name_ptf)

% Build a struct for the output
Output_struct = struct('Volatility', maxSharpeRisk_Ptf, 'Weights', maxSharpeWgt_Ptf,...
    'Return', maxSharpeRet_Ptf, 'Sharpe_Ratio', maxSharpeSR_Ptf);
Output_struct.Name = name_ptf;
Output_struct.Ptf = Ptf;

end