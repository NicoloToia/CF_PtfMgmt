function [P1,minRisk_P1, minRiskWgt_P1, minRiskRet_P1, minRiskSR_P1] = ...
    minRiskPortfolio(names, mean_returns, cov_matrix, risk_free_rate , const,flag )

    %   Compute the Minimum Variance Portfolio of the frontier.
    % INPUTS
    % names: cell array of asset names
    % mean_returns: vector of mean returns
    % cov_matrix: covariance matrix
    % risk_free_rate: risk free rate
    % const: constraints struct 

    % Create a portfolio object
    P1 = Portfolio('AssetList', names);
    P1 = setDefaultConstraints(P1); % all weights sum to 1, no shorting, and 100% investment in risky assets
    P1 = setAssetMoments(P1, mean_returns, cov_matrix); % set mean returns and covariance matrix

    P1.AInequality = const.Aineq;
    P1.bInequality = const.bineq;
    % estimate efficient frontier, via object method and 100 P1oints
    pwgt1 = estimateFrontier(P1, 100);
    % estimate portfolio moments
    [pf_risk_P1, pf_Retn_P1] = estimatePortMoments(P1, pwgt1);
    % plot efficient frontier
    % plot(pf_risk_P1, pf_Retn_P1)

    % find minimum risk portfolio (Portfolio A or C - Minimum Variance Portfolio)
    minRisk_P1 = min(pf_risk_P1);
    minRiskWgt_P1 = pwgt1(:, pf_risk_P1 == minRisk_P1);
    minRiskRet_P1 = mean_returns' * minRiskWgt_P1;
    minRiskSR_P1 = (minRiskRet_P1 - risk_free_rate) / minRisk_P1;



    % Display Portfolio A or C - Minimum Variance Portfolio
    if flag == 0
        disp('===========================================================================')
        disp(' Maximum Sharpe Ratio Portfolio (Portfolio A)')
        disp('===========================================================================')
    elseif flag == 1
        disp('===========================================================================')
        disp(' Maximum Sharpe Ratio Portfolio (Portfolio C)')
        disp('===========================================================================')
    end
    print_portfolio(minRiskWgt_P1, names, minRiskRet_P1, minRisk_P1, minRiskSR_P1)
end