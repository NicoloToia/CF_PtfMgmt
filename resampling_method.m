function [minRisk_P1_Rsim, minRiskWgt_P1_Rsim, minRiskRet_P1_Rsim, minRiskSR_P1_Rsim, ...
    minRisk_P2_Rsim, minRiskWgt_P2_Rsim, minRiskRet_P2_Rsim, minRiskSR_P2_Rsim, ...
    maxSharpeSR_P1_Rsim, maxSharpeWgt_P1_Rsim, maxSharpeRet_P1_Rsim, maxSharpeRisk_P1_Rsim, ...
    maxSharpeSR_P2_Rsim, maxSharpeWgt_P2_Rsim, maxSharpeRet_P2_Rsim, maxSharpeRisk_P2_Rsim] =...
     resampling_method(mean_returns, cov_matrix, P1, P2, risk_free_rate,num_assets,names)

    %  Compute the frontiers in step 1 and 2 using the resampling method in
    %  order to obtain 2 robust frontier. For each frontier save the Minimum
    %  Variance Portfolios, named Portfolios E and F, and the Maximum
    %  Sharpe Ratio Portfolios, named Portfolios G and H, of the frontiers.

    % INPUTS
    % mean_returns: expected returns of the assets
    % cov_matrix: covariance matrix of the assets
    % P1: Portfolio object with no constraints
    % P2: Portfolio object with constraints
    % risk_free_rate: risk free rate
    % num_assets: number of assets
    % names: names of the assets

    % number of resampling simulations
    N = 100;
    % # of points along the efficient frontier
    num_frontier_points = 100;

    % Initialize variables
    Ret_P1_sim = zeros(num_frontier_points, N);
    Risk_P1_sim = zeros(num_frontier_points, N);
    Weights_P1_sim = zeros(num_assets, num_frontier_points, N);

    Ret_P2_sim = zeros(num_frontier_points, N);
    Risk_P2_sim = zeros(num_frontier_points, N);
    Weights_P2_sim = zeros(num_assets, num_frontier_points, N);

    for n = 1:N
        
        resampledReturns = mvnrnd(mean_returns, cov_matrix, 252); % 252 mi suggerisce chat per resaplare i rendimenti giornalieri i.e.  252 giorni
        New_mean_returns = mean(resampledReturns)';               % altrimenti mi dice che non ho un resampling completo ma non sono convinto
        NewCov = cov(resampledReturns);

        P1_sim = setAssetMoments(P1, New_mean_returns, NewCov);
        P2_sim = setAssetMoments(P2, New_mean_returns, NewCov);

        w1_sim = estimateFrontier(P1_sim, num_frontier_points);
        w2_sim = estimateFrontier(P2_sim, num_frontier_points);

        [pf_risk_P1_sim, pf_Retn_P1_sim] = estimatePortMoments(P1_sim, w1_sim);
        [pf_risk_P2_sim, pf_Retn_P2_sim] = estimatePortMoments(P2_sim, w2_sim);

        Ret_P1_sim(:,n) = pf_Retn_P1_sim;
        Risk_P1_sim(:, n) = pf_risk_P1_sim;
        Weights_P1_sim(:,:,n) = w1_sim;

        Ret_P2_sim(:,n) = pf_Retn_P2_sim;
        Risk_P2_sim(:, n) = pf_risk_P2_sim;
        Weights_P2_sim(:,:,n) = w2_sim;
        
    end

    pwgt_1 = mean(Weights_P1_sim, 3);
    pf_risk_1 = mean(Risk_P1_sim, 2);
    pf_Retn_1 = mean(Ret_P1_sim, 2);

    pwgt_2 = mean(Weights_P2_sim, 3);
    pf_risk_2 = mean(Risk_P2_sim, 2);
    pf_Retn_2 = mean(Ret_P2_sim, 2);

    % Find minimum risk portfolio (Portfolio E - Minimum Variance Portfolio with resampling)
    [minRisk_P1_Rsim, idx_minRisk_P1] = min(pf_risk_1);
    minRiskWgt_P1_Rsim = pwgt_1(:, idx_minRisk_P1);
    minRiskRet_P1_Rsim = pf_Retn_1(idx_minRisk_P1);
    minRiskSR_P1_Rsim = (minRiskRet_P1_Rsim - risk_free_rate) / minRisk_P1_Rsim;

    % Find minimum risk portfolio (Portfolio F - Minimum Variance Portfolio with resampling and constraints)
    [minRisk_P2_Rsim, idx_minRisk_P2] = min(pf_risk_2);
    minRiskWgt_P2_Rsim = pwgt_2(:, idx_minRisk_P2);
    minRiskRet_P2_Rsim = pf_Retn_2(idx_minRisk_P2);
    minRiskSR_P2_Rsim = (minRiskRet_P2_Rsim - risk_free_rate) / minRisk_P2_Rsim;

    % Find maximum Sharpe ratio portfolio (Portfolio G - Maximum Sharpe Ratio Portfolio with resampling)
    sharpeRatio_P1_Rsim = (pf_Retn_1 - risk_free_rate) ./ pf_risk_1;
    [maxSharpeSR_P1_Rsim, idx_maxSharpe_P1] = max(sharpeRatio_P1_Rsim);
    maxSharpeWgt_P1_Rsim = pwgt_1(:, idx_maxSharpe_P1);
    maxSharpeRet_P1_Rsim = pf_Retn_1(idx_maxSharpe_P1);
    maxSharpeRisk_P1_Rsim = pf_risk_1(idx_maxSharpe_P1);

    % Find maximum Sharpe ratio portfolio (Portfolio H - Maximum Sharpe Ratio Portfolio with resampling and constraints)
    sharpeRatio_P2_Rsim = (pf_Retn_2 - risk_free_rate) ./ pf_risk_2;
    [maxSharpeSR_P2_Rsim, idx_maxSharpe_P2] = max(sharpeRatio_P2_Rsim);
    maxSharpeWgt_P2_Rsim = pwgt_2(:, idx_maxSharpe_P2);
    maxSharpeRet_P2_Rsim = pf_Retn_2(idx_maxSharpe_P2);
    maxSharpeRisk_P2_Rsim = pf_risk_2(idx_maxSharpe_P2);

    % Display Portfolio E - Minimum Variance Portfolio with resampling
    print_portfolio(minRiskWgt_P1_Rsim, names, minRiskRet_P1_Rsim, minRisk_P1_Rsim, minRiskSR_P1_Rsim,'E')
    
    % Display Portfolio F - Minimum Variance Portfolio with resampling and constraints
    print_portfolio(minRiskWgt_P2_Rsim, names, minRiskRet_P2_Rsim, minRisk_P2_Rsim, minRiskSR_P2_Rsim,'F')

    % Display Portfolio G - Maximum Sharpe Ratio Portfolio with resampling
    print_portfolio(maxSharpeWgt_P1_Rsim, names, maxSharpeRet_P1_Rsim, maxSharpeRisk_P1_Rsim, maxSharpeSR_P1_Rsim,'G')

    % Display Portfolio H - Maximum Sharpe Ratio Portfolio with resampling and constraints
    print_portfolio(maxSharpeWgt_P2_Rsim, names, maxSharpeRet_P2_Rsim, maxSharpeRisk_P2_Rsim, maxSharpeSR_P2_Rsim,'H')
end