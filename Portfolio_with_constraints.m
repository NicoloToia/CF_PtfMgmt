function [P2, minRisk_P2, minRiskWgt_P2, minRiskRet_P2,...
     minRiskSR_P2, maxSharpeWgt_P2, maxSharpeRet_P2, maxSharpeRisk_P2, maxSharpeSR_P2] =...
      Portfolio_with_constraints(mean_returns, cov_matrix, names, risk_free_rate, mkt)


%  Compute the efficient frontier under the following constraints (to be
%  considered all at once):
%  • Standard constraints,
%  • The total exposure to sensible sectors has to be greater than 10%
%  and the total exposure on cyclical sectors has to be lower than
%  30%,
%  • The weights of the sector Consumer Staples and the factor
%  Low Volatility have to be set equal to 0,
%  • The maximum exposure on sectors has to be lower than 80%.
%  Compute the Minimum Variance Portfolio, named Portfolio C, and the
%  Maximum Sharpe Ratio Portfolio, named Portfolio D, of the frontier.


sectors = [mkt.sector.sensible, mkt.sector.cyclical, mkt.sector.defensive];

% Create a portfolio object
P2 = Portfolio('AssetList', names);
P2 = setAssetMoments(P2, mean_returns, cov_matrix);
% Define the additional constraints
% 1. Standard constraints, all weights sum to 1, no shorting, and 100% investment in risky assets
% 2. Total exposure to sensible sectors (Energy, Information Technology, Communication Services) 
%    should be greater than 10% and total exposure to cyclical sectors should be lower than 30%
% 3. Consumer Staples and Low Volatility should not be included in the portfolio (weights set to 0)
% 4. Maximum exposure to sectors (Sensible, Cyclical, Defensive) should be less than 80%

sensibleIdx = ismember(P2.AssetList, mkt.sector.sensible);
cyclicalIdx = ismember(P2.AssetList, mkt.sector.cyclical);
excludeIdx_CS = ismember(P2.AssetList, "ConsumerStaples");
excludeIdx_LV = ismember(P2.AssetList, "LowVolatility");
sectorIdx = ismember(P2.AssetList, sectors);

Aineq = [
    -sensibleIdx;
    cyclicalIdx;
    excludeIdx_CS;
    excludeIdx_LV;
    sectorIdx;
];

bineq = [-0.1; 0.3; 0; 0; 0.8];
% set the constraints
P2 = setDefaultConstraints(P2);
P2.AInequality = Aineq;
P2.bInequality = bineq;

% Estimate efficient frontier with constraints
pwgt2 = estimateFrontier(P2, 100);
% Estimate portfolio moments
[pf_risk_P2, pf_Retn_P2] = estimatePortMoments(P2, pwgt2);
% Plot efficient frontier with constraints
% plot(pf_risk_P2, pf_Retn_P2)


% Find minimum risk portfolio (Portfolio C - Minimum Variance Portfolio with constraints)
minRisk_P2 = min(pf_risk_P2);
minRiskWgt_P2 = pwgt2(:, pf_risk_P2 == minRisk_P2);
minRiskRet_P2 = mean_returns' * minRiskWgt_P2;
minRiskSR_P2 = (minRiskRet_P2 - risk_free_rate) / minRisk_P2;



% Find maximum Sharpe ratio portfolio (Portfolio D - Maximum Sharpe Ratio Portfolio with constraints)
maxSharpeWgt_P2 = estimateMaxSharpeRatio(P2);
maxSharpeRet_P2 = mean_returns' * maxSharpeWgt_P2;
maxSharpeRisk_P2 = sqrt(maxSharpeWgt_P2' * cov_matrix * maxSharpeWgt_P2);
maxSharpeSR_P2 = (maxSharpeRet_P2 - risk_free_rate) / maxSharpeRisk_P2;

% Display Portfolio C - Minimum Variance Portfolio with constraints
disp('===========================================================================')
disp('Minimum Risk Portfolio with Constraints (Portfolio C)')
disp('===========================================================================')
print_portfolio(minRiskWgt_P2, names, minRiskRet_P2, minRisk_P2, minRiskSR_P2)

% Display Portfolio D - Maximum Sharpe Ratio Portfolio with constraints
disp('===========================================================================')
disp('Maximum Sharpe Ratio Portfolio with Constraints (Portfolio D)')
disp('===========================================================================')
print_portfolio(maxSharpeWgt_P2, names, maxSharpeRet_P2, maxSharpeRisk_P2, maxSharpeSR_P2)

end