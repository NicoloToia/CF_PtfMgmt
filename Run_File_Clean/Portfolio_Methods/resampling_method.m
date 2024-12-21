function Output_struct = resampling_method(Ptf, name_ptf, flag)
% This function estimates the minimum risk portfolio or the maximum Sharpe ratio portfolio
% using the resampling method.
%
% INPUTS:
% Ptf:          Portfolio object
% name_ptf:     The name of the portfolio
% flag:         0 for minimum risk portfolio
%               1 for maximum Sharpe ratio portfolio
%
% OUTPUT:
% Output_struct: A struct containing the portfolio weights,
%                expected return, volatility, and Sharpe ratio

% fix the seed for reproducibility between the portfolios
rng(42); % the answer to the ultimate question of life, the universe, and everything

% Get the number of assets in the portfolio
num_assets = Ptf.NumAssets;

% number of resampling simulations
N = 500;
% # of points along the efficient frontier
num_frontier_points = 100;

% Initialize variables
Ret_sim = zeros(num_frontier_points, N);
Risk_sim = zeros(num_frontier_points, N);
Weights_sim = zeros(num_assets, num_frontier_points, N);

for n = 1:N
    % Resample the returns and6 estimate the frontier
    resampledReturns = mvnrnd(Ptf.AssetMean, Ptf.AssetCovar, 252);
    New_mean_returns = mean(resampledReturns)';
    NewCov = cov(resampledReturns);

    Ptf_sim = setAssetMoments(Ptf, New_mean_returns, NewCov);

    w_sim = estimateFrontier(Ptf_sim, num_frontier_points);

    [pf_risk_sim, pf_Retn_sim] = estimatePortMoments(Ptf_sim, w_sim);

    Ret_sim(:,n) = pf_Retn_sim;
    Risk_sim(:, n) = pf_risk_sim;
    Weights_sim(:,:,n) = w_sim;
  
end

% Compute the average of the resampled portfolios
pwgt = mean(Weights_sim, 3);
pf_risk = mean(Risk_sim, 2);
pf_Retn = mean(Ret_sim, 2);

if flag == 0    % Minimum risk portfolio

    % Find the minimum risk portfolio (Minimum Variance Portfolio with resampling)
    [minRisk_Rsim, idx_minRisk] = min(pf_risk);
    minRiskWgt_Rsim = pwgt(:, idx_minRisk);
    minRiskRet_Rsim = pf_Retn(idx_minRisk);
    minRiskSR_Rsim = (minRiskRet_Rsim - Ptf.RiskFreeRate) / minRisk_Rsim;

    % Display the minimum risk portfolio
    print_portfolio(minRiskWgt_Rsim, Ptf.AssetList, minRiskRet_Rsim, minRisk_Rsim, minRiskSR_Rsim, name_ptf);

    % Build a struct for the output
    Output_struct = struct('Volatility', minRisk_Rsim, 'Weights', minRiskWgt_Rsim,...
        'Return', minRiskRet_Rsim, 'Sharpe_Ratio', minRiskSR_Rsim);
    Output_struct.Name = name_ptf;
    Output_struct.Ptf = Ptf;
    
elseif flag == 1    % Maximum Sharpe ratio portfolio

    % Find the maximum Sharpe ratio portfolio (Maximum Sharpe Ratio Portfolio with resampling)
    sharpeRatio_Rsim = (pf_Retn - Ptf.RiskFreeRate) ./ pf_risk;
    [maxSharpeSR_Rsim, idx_maxSharpe] = max(sharpeRatio_Rsim);
    maxSharpeWgt_Rsim = pwgt(:, idx_maxSharpe);
    maxSharpeRet_Rsim = pf_Retn(idx_maxSharpe);
    maxSharpeRisk_Rsim = pf_risk(idx_maxSharpe);

    % Display the maximum Sharpe ratio portfolio
    print_portfolio(maxSharpeWgt_Rsim, Ptf.AssetList, maxSharpeRet_Rsim, maxSharpeRisk_Rsim, maxSharpeSR_Rsim, name_ptf);

    % Build a struct for the output
    Output_struct = struct('Volatility', maxSharpeRisk_Rsim, 'Weights', maxSharpeWgt_Rsim,...
        'Return', maxSharpeRet_Rsim, 'Sharpe_Ratio', maxSharpeSR_Rsim);
    Output_struct.Name = name_ptf;
    Output_struct.Ptf = Ptf;

else
    % Error message, not valid flag
    error('Not valid flag, please use 0 for minimum risk portfolio and 1 for maximum Sharpe ratio portfolio')
end

end