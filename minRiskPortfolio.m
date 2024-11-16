function [P1,minRisk_P1, minRiskWgt_P1, minRiskRet_P1, minRiskSR_P1] = minRiskPortfolio(names, mean_returns, cov_matrix, risk_free_rate)
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

% find minimum risk portfolio (Portfolio A - Minimum Variance Portfolio)
minRisk_P1 = min(pf_risk_P1);
minRiskWgt_P1 = pwgt1(:, pf_risk_P1 == minRisk_P1);
minRiskRet_P1 = mean_returns' * minRiskWgt_P1;
minRiskSR_P1 = (minRiskRet_P1 - risk_free_rate) / minRisk_P1;

% Display Portfolio A - Minimum Variance Portfolio
disp('===========================================================================')
disp('    Minimum Risk Portfolio (Portfolio A)   ')
disp('===========================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(minRiskWgt_P1)
    fprintf('%-25s %.4f\n', names{i}, minRiskWgt_P1(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', minRiskRet_P1);
fprintf('%-25s %.4f\n', 'Volatility', minRisk_P1);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', minRiskSR_P1);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(minRiskWgt_P1));
disp('  ')