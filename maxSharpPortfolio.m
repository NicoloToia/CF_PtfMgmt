function [P1,maxSharpeRisk_P1, maxSharpeWgt_P1, maxSharpeRet_P1, maxSharpeSR_P1] = maxSharpPortfolio(names, mean_returns, cov_matrix, risk_free_rate)
%   Compute the efficient frontier under the standard constraints.
%   Compute the Minimum Variance Portfolio, named Portfolio A, and the Maximum Sharpe 
%   Ratio Portfolio, named Portfolio B, of the frontier.

% Create a portfolio object
P1 = Portfolio('AssetList', names);
P1 = setDefaultConstraints(P1); % all weights sum to 1, no shorting, and 100% investment in risky assets
P1 = setAssetMoments(P1, mean_returns, cov_matrix); % set mean returns and covariance matrix
% estimate efficient frontier, via object method and 100 P1oints
pwgt1 = estimateFrontier(P1, 100);
% estimate portfolio moments
[pf_risk_P1, pf_Retn_P1] = estimatePortMoments(P1, pwgt1);
% plot efficient frontier
% plot(pf_risk_P1, pf_Retn_P1)

% find maximum Sharpe ratio portfolio (Portfolio B - Maximum Sharpe Ratio Portfolio)
maxSharpeWgt_P1 = estimateMaxSharpeRatio(P1);
maxSharpeRet_P1 = mean_returns' * maxSharpeWgt_P1;
maxSharpeRisk_P1 = sqrt(maxSharpeWgt_P1' * cov_matrix * maxSharpeWgt_P1);
maxSharpeSR_P1 = (maxSharpeRet_P1 - risk_free_rate) / maxSharpeRisk_P1;

% Display Portfolio B - Maximum Sharpe Ratio Portfolio
disp('===========================================================================')
disp(' Maximum Sharpe Ratio Portfolio (Portfolio B)')
disp('===========================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(maxSharpeWgt_P1)
    fprintf('%-25s %.4f\n', names{i}, maxSharpeWgt_P1(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', maxSharpeRet_P1);
fprintf('%-25s %.4f\n', 'Volatility', maxSharpeRisk_P1);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', maxSharpeSR_P1);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(maxSharpeWgt_P1));
disp('  ')